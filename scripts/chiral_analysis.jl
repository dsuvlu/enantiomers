using CliffordNumbers
using Chemfiles
using LinearAlgebra
using DelimitedFiles
using BenchmarkTools
using StaticArrays
using ArgParse
using Serialization

include("/home/dsuvlu/git/enantiomers/src/julia/main.jl")

args = parse_command_line()
directory = args["directory"]
traj_name = args["traj"]
pdb_name = args["pdb"]
output = args["output"]
solute_flag = args["solute"]
#solute_flag = true
#solute_flag = false

# set parameters
params = Parameters(3.5, 30.0, SVector{1,Float64}(0.87243), 0.05, 3.5)

# initialize histograms and per frame data structures
trivec_hist_total, trivec_hist_frame, 
mean_trivec_hist_total, mean_trivec_hist_frame, chain_hist = initialize_histograms(params, solute_flag)

trivec_sum_frame, trivec_sum_total, 
mean_trivec_sum_frame, mean_trivec_sum_total = initialize_statistics(params, solute_flag)

# clifford algebra basis vectors
@basis_vars(VGA3D, Int)

# load trajectory, determine number of frames
#directory = "/home/dsuvlu/git/enantiomers/data/alanine/"
#directory = "/home/dsuvlu/git/enantiomers/data/water/"
#output = "/home/dsuvlu/git/enantiomers/results/water/"
#output = "/home/dsuvlu/git/enantiomers/results/alanine/r/"
#traj_name = "trajectory_10ns.dcd"
#traj_name = "r_trajectory_10ns.dcd"
#pdb_name = "water.pdb"
#pdb_name = "r_alanine_water.pdb"
trajfile = joinpath(directory, traj_name)
pdbfile = joinpath(directory, pdb_name)
trajectory = Trajectory(trajfile)
set_topology!(trajectory, pdbfile)
n_frames = Int(length(trajectory))

# grab frame, number of atoms
f = 0
frame = read_step(trajectory, f)
n_atoms = Int(length(frame))

# select all water oxygens
sel_oxy = Selection("name O")
oxy = Int.(evaluate(sel_oxy, frame))
n_o = length(oxy)

# initialize dictionaries
nhbs, donors, chains, trivecs, min_dist = initialize_dictionaries(oxy)

# select hydrogen bonding solute atoms, initialize solute dictionary
if solute_flag
    sel_sol = Selection("(type O or type N) and not resname HOH")
    solute = Int.(evaluate(sel_sol, frame))
    sol_dist = Dict{Int, Float64}()
    for i in solute
        sol_dist[i] = 0.0
    end
end

# loop over all frames
for f in 0:n_frames-1

    # grab current frame
    println("Frame: ", f)
    frame = read_step(trajectory, f)
    
    # grab current unit cell
    unitcell = lengths(UnitCell(frame))
    ucell = SVector{3, Float64}(unitcell)

    # grab current positions, convert to static arrays
    pos = positions(frame)
    pos_static = Array{SVector{3,Float64}}(undef, n_atoms)
    for i in 1:n_atoms
        pos_static[i] = SVector{3, Float64}(pos[:,i])
    end

    # clear dictionaries
    clear_values!(oxy, nhbs, donors, chains, trivecs, min_dist)

    # compute minimum distance between solute and water oxygens
    if solute_flag
        minimum_distance!(min_dist, pos_static, ucell, oxy, solute, sol_dist)
    end
    
    # find neighbors within dOO distance
    #@time find_neighbors!(nhbs, oxy, n_o, pos_static, params, ucell)
    find_neighbors_fast!(nhbs, oxy, n_o, pos_static, params, ucell)

    # find hydrogen bond donors
    find_donors!(donors, nhbs, n_o, oxy, frame, params)

    # find hydrogen bond chains
    find_hbond_chains!(donors, n_o, oxy, chains)

    # compute trivectors
    compute_trivectors!(pos_static, params, chains, trivecs, e1, e2, e3)

    # compute per atom mean trivectors at pi
    p1_mean = compute_mean_pi(trivecs.p1)
    p2_mean = compute_mean_pi(trivecs.p2)
    p3_mean = compute_mean_pi(trivecs.p3)
    p4_mean = compute_mean_pi(trivecs.p4)
    all_mean = compute_mean_all(trivecs, oxy)

    # reset per frame values 
    reset_statistics!(trivec_sum_frame, solute_flag)
    reset_statistics!(trivec_hist_frame, solute_flag)
    reset_statistics!(mean_trivec_sum_frame, solute_flag)    
    reset_statistics!(mean_trivec_hist_frame, solute_flag)

    # per atom mean flag
    mean_flag = true
    # compute histograms for per atom mean trivectors
    histogram_trivecs!(p1_mean, mean_trivec_sum_frame.p1, mean_trivec_hist_frame.p1, 
        mean_trivec_hist_frame.bins1, mean_trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(p2_mean, mean_trivec_sum_frame.p2, mean_trivec_hist_frame.p2, 
        mean_trivec_hist_frame.bins1, mean_trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(p3_mean, mean_trivec_sum_frame.p3, mean_trivec_hist_frame.p3, 
        mean_trivec_hist_frame.bins1, mean_trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(p4_mean, mean_trivec_sum_frame.p4, mean_trivec_hist_frame.p4, 
        mean_trivec_hist_frame.bins1, mean_trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(all_mean, mean_trivec_sum_frame.all, mean_trivec_hist_frame.all, 
        mean_trivec_hist_frame.bins1, mean_trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    # per atom mean flag
    mean_flag = false
    # compute histograms for trivectors
    histogram_trivecs!(trivecs.p1, trivec_sum_frame.p1, trivec_hist_frame.p1, 
        trivec_hist_frame.bins1, trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(trivecs.p2, trivec_sum_frame.p2, trivec_hist_frame.p2, 
        trivec_hist_frame.bins1, trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(trivecs.p3, trivec_sum_frame.p3, trivec_hist_frame.p3, 
        trivec_hist_frame.bins1, trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    histogram_trivecs!(trivecs.p4, trivec_sum_frame.p4, trivec_hist_frame.p4, 
        trivec_hist_frame.bins1, trivec_hist_frame.bins2, min_dist, 
        solute_flag, mean_flag)

    # update total statistics 
    update_totals!(mean_trivec_sum_total, mean_trivec_sum_frame)
    update_totals!(mean_trivec_hist_total, mean_trivec_hist_frame)
    update_totals!(trivec_sum_total, trivec_sum_frame)
    update_totals!(trivec_hist_total, trivec_hist_frame)

    # compute histograms of chain participation
    histogram_chain_participation!(trivecs.p1, chain_hist.p1, 
        chain_hist.bins1, chain_hist.bins2, min_dist, solute_flag)
    histogram_chain_participation!(trivecs.p2, chain_hist.p2, 
        chain_hist.bins1, chain_hist.bins2, min_dist, solute_flag)
    histogram_chain_participation!(trivecs.p3, chain_hist.p3, 
        chain_hist.bins1, chain_hist.bins2, min_dist, solute_flag)
    histogram_chain_participation!(trivecs.p4, chain_hist.p4, 
        chain_hist.bins1, chain_hist.bins2, min_dist, solute_flag)

    # save chain info for chains close to the enantiomer to a binary file
    if solute_flag
        chain_info = ChainList(Dict{Int, Vector{Chains}}(), Dict{Int, Vector{Chains}}(),
        Dict{Int, Vector{Chains}}(), Dict{Int, Vector{Chains}}())
        write_chains!(chain_info, f, trivecs, min_dist, chains, params, output)
    end

    # save the average trivector and per atom mean trivector for each frame
    write_trivec_frame!(trivec_sum_frame, trivec_hist_frame, 
        mean_trivec_sum_frame, mean_trivec_hist_frame, output)
end

# write histograms to file
if !solute_flag
    
    writedlm(joinpath(output, "chain_hist.dat"), 
        hcat(chain_hist.p1, chain_hist.p2, chain_hist.p3, chain_hist.p4))

    writedlm(joinpath(output, "trivec_time_avg.dat"), 
        hcat(trivec_sum_total.p1./vec(sum(trivec_hist_total.p1, dims=1)), 
        trivec_sum_total.p2./vec(sum(trivec_hist_total.p2, dims=1)), 
        trivec_sum_total.p3./vec(sum(trivec_hist_total.p3, dims=1)), 
        trivec_sum_total.p4./vec(sum(trivec_hist_total.p4, dims=1))))

    writedlm(joinpath(output, "trivec_hist.dat"), 
        hcat(trivec_hist_total.p1, trivec_hist_total.p2, 
        trivec_hist_total.p3, trivec_hist_total.p4))

    writedlm(joinpath(output, "mean_trivec_time_avg.dat"),
        hcat(mean_trivec_sum_total.p1./vec(sum(mean_trivec_hist_total.p1, dims=1)), 
        mean_trivec_sum_total.p2./vec(sum(mean_trivec_hist_total.p2, dims=1)), 
        mean_trivec_sum_total.p3./vec(sum(mean_trivec_hist_total.p3, dims=1)), 
        mean_trivec_sum_total.p4./vec(sum(mean_trivec_hist_total.p4, dims=1)), 
        mean_trivec_sum_total.all./vec(sum(mean_trivec_hist_total.all, dims=1))))

    writedlm(joinpath(output, "mean_trivec_hist.dat"), 
        hcat(mean_trivec_hist_total.p1, mean_trivec_hist_total.p2, 
        mean_trivec_hist_total.p3, mean_trivec_hist_total.p4, mean_trivec_hist_total.all))

elseif solute_flag
    
    writedlm(joinpath(output, "chain_hist_p1.dat"), chain_hist.p1)
    writedlm(joinpath(output, "chain_hist_p2.dat"), chain_hist.p2)
    writedlm(joinpath(output, "chain_hist_p3.dat"), chain_hist.p3)
    writedlm(joinpath(output, "chain_hist_p4.dat"), chain_hist.p4)

    writedlm(joinpath(output, "trivec_time_avg_p1.dat"), 
        trivec_sum_total.p1./vec(sum(trivec_hist_total.p1, dims=1)))
    writedlm(joinpath(output, "trivec_time_avg_p2.dat"), 
        trivec_sum_total.p2./vec(sum(trivec_hist_total.p2, dims=1)))
    writedlm(joinpath(output, "trivec_time_avg_p3.dat"), 
        trivec_sum_total.p3./vec(sum(trivec_hist_total.p3, dims=1)))
    writedlm(joinpath(output, "trivec_time_avg_p4.dat"), 
        trivec_sum_total.p4./vec(sum(trivec_hist_total.p4, dims=1)))

    writedlm(joinpath(output, "trivec_hist_p1.dat"), trivec_hist_total.p1)
    writedlm(joinpath(output, "trivec_hist_p2.dat"), trivec_hist_total.p2)
    writedlm(joinpath(output, "trivec_hist_p3.dat"), trivec_hist_total.p3)
    writedlm(joinpath(output, "trivec_hist_p4.dat"), trivec_hist_total.p4)

    writedlm(joinpath(output, "mean_trivec_time_avg_p1.dat"), 
        mean_trivec_sum_total.p1./vec(sum(mean_trivec_hist_total.p1, dims=1)))
    writedlm(joinpath(output, "mean_trivec_time_avg_p2.dat"),
        mean_trivec_sum_total.p2./vec(sum(mean_trivec_hist_total.p2, dims=1)))
    writedlm(joinpath(output, "mean_trivec_time_avg_p3.dat"),
        mean_trivec_sum_total.p3./vec(sum(mean_trivec_hist_total.p3, dims=1)))
    writedlm(joinpath(output, "mean_trivec_time_avg_p4.dat"),
        mean_trivec_sum_total.p4./vec(sum(mean_trivec_hist_total.p4, dims=1)))
    writedlm(joinpath(output, "mean_trivec_time_avg_all.dat"),
        mean_trivec_sum_total.all./vec(sum(mean_trivec_hist_total.all, dims=1)))

    writedlm(joinpath(output, "mean_trivec_hist_p1.dat"), mean_trivec_hist_total.p1)
    writedlm(joinpath(output, "mean_trivec_hist_p2.dat"), mean_trivec_hist_total.p2)
    writedlm(joinpath(output, "mean_trivec_hist_p3.dat"), mean_trivec_hist_total.p3)
    writedlm(joinpath(output, "mean_trivec_hist_p4.dat"), mean_trivec_hist_total.p4)
    writedlm(joinpath(output, "mean_trivec_hist_all.dat"), mean_trivec_hist_total.all)

end
