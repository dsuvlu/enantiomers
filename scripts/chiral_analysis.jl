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
#solute_flag = "yes"

# set parameters
params = Parameters(3.5, 30.0, SVector{1,Float64}(0.87243), 0.05, 3.5)

# initialize histograms
trivec_hist, mean_trivec_hist, chain_hist = initialize_histograms(params, solute_flag)

# clifford algebra basis vectors
@basis_vars(VGA3D, Int)

# load trajectory, determine number of frames
#directory = "/home/dsuvlu/git/enantiomers/data/alanine/"
#output = "/home/dsuvlu/git/enantiomers/results/alanine/r/"
#traj_name = "r_trajectory_10ns.dcd"
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
nhbs, donors, chains, trivecs = initialize_dictionaries(oxy)

# select hydrogen bonding solute atoms
if solute_flag == "yes"
    sel_sol = Selection("(type O or type N) and not resname HOH")
    solute = Int.(evaluate(sel_sol, frame))
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

    # compute minimum distance between solute and water oxygens
    if solute_flag == "yes"
        min_dist = @time minimum_distance(pos_static, ucell, oxy, solute)
    else
        min_dist = Dict{Int, Vector{Float64}}()
    end

    # clear dictionaries
    clear_values!(oxy, nhbs, donors, chains, trivecs)
    
    # find neighbors within dOO distance
    #@time find_neighbors!(nhbs, oxy, n_o, pos_static, params, ucell)
    @time find_neighbors_fast!(nhbs, oxy, n_o, pos_static, params, ucell)

    # find hydrogen bond donors
    @time find_donors!(donors, nhbs, n_o, oxy, frame, params)

    # find hydrogen bond chains
    @time find_hbond_chains!(donors, n_o, oxy, chains)

    # compute trivectors
    @time compute_trivectors!(pos_static, params, chains, trivecs)

    # compute mean trivectors
    p1_mean = compute_mean_pi(trivecs.p1)
    p2_mean = compute_mean_pi(trivecs.p2)
    p3_mean = compute_mean_pi(trivecs.p3)
    p4_mean = compute_mean_pi(trivecs.p4)
    all_mean = compute_mean_all(trivecs, oxy)

    # compute histograms for mean trivectors
    histogram_trivecs!(mean_trivec_hist.p1, p1_mean, mean_trivec_hist.bins1, mean_trivec_hist.bins2, min_dist)
    histogram_trivecs!(mean_trivec_hist.p2, p2_mean, mean_trivec_hist.bins1, mean_trivec_hist.bins2, min_dist)
    histogram_trivecs!(mean_trivec_hist.p3, p3_mean, mean_trivec_hist.bins1, mean_trivec_hist.bins2, min_dist)
    histogram_trivecs!(mean_trivec_hist.p4, p4_mean, mean_trivec_hist.bins1, mean_trivec_hist.bins2, min_dist)
    histogram_trivecs!(mean_trivec_hist.all, all_mean, mean_trivec_hist.bins1, mean_trivec_hist.bins2, min_dist)

    # compute histograms for trivectors
    histogram_trivecs!(trivec_hist.p1, trivecs.p1, trivec_hist.bins1, trivec_hist.bins2, min_dist)
    histogram_trivecs!(trivec_hist.p2, trivecs.p2, trivec_hist.bins1, trivec_hist.bins2, min_dist)
    histogram_trivecs!(trivec_hist.p3, trivecs.p3, trivec_hist.bins1, trivec_hist.bins2, min_dist)
    histogram_trivecs!(trivec_hist.p4, trivecs.p4, trivec_hist.bins1, trivec_hist.bins2, min_dist)

    # compute histograms of chain participation
    histogram_chain_participation!(chain_hist.p1, trivecs.p1, chain_hist.bins1, chain_hist.bins2, min_dist)
    histogram_chain_participation!(chain_hist.p2, trivecs.p2, chain_hist.bins1, chain_hist.bins2, min_dist)
    histogram_chain_participation!(chain_hist.p3, trivecs.p3, chain_hist.bins1, chain_hist.bins2, min_dist)
    histogram_chain_participation!(chain_hist.p4, trivecs.p4, chain_hist.bins1, chain_hist.bins2, min_dist)

    # save chain info for chains close to the enantiomer to a binary file
    if solute_flag == "yes"
        chain_info = Chains(Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}())
        write_chains!(chain_info, f, trivecs, min_dist, chains, params, output)
    end
end

# write histograms to file
if solute_flag == "no"
    
    writedlm(joinpath(output, "chain_hist.dat"), 
        hcat(chain_hist.p1, chain_hist.p2, chain_hist.p3, chain_hist.p4))

    writedlm(joinpath(output, "trivec_hist.dat"), 
        hcat(trivec_hist.p1, trivec_hist.p2, trivec_hist.p3, trivec_hist.p4))

    writedlm(joinpath(output, "mean_trivec_hist.dat"), 
        hcat(mean_trivec_hist.p1, mean_trivec_hist.p2, mean_trivec_hist.p3, 
        mean_trivec_hist.p4, mean_trivec_hist.all))

elseif solute_flag == "yes"
    
    writedlm(joinpath(output, "chain_hist_p1.dat"), chain_hist.p1)
    writedlm(joinpath(output, "chain_hist_p2.dat"), chain_hist.p2)
    writedlm(joinpath(output, "chain_hist_p3.dat"), chain_hist.p3)
    writedlm(joinpath(output, "chain_hist_p4.dat"), chain_hist.p4)

    writedlm(joinpath(output, "trivec_hist_p1.dat"), trivec_hist.p1)
    writedlm(joinpath(output, "trivec_hist_p2.dat"), trivec_hist.p2)
    writedlm(joinpath(output, "trivec_hist_p3.dat"), trivec_hist.p3)
    writedlm(joinpath(output, "trivec_hist_p4.dat"), trivec_hist.p4)

    writedlm(joinpath(output, "mean_trivec_hist_p1.dat"), mean_trivec_hist.p1)
    writedlm(joinpath(output, "mean_trivec_hist_p2.dat"), mean_trivec_hist.p2)
    writedlm(joinpath(output, "mean_trivec_hist_p3.dat"), mean_trivec_hist.p3)
    writedlm(joinpath(output, "mean_trivec_hist_p4.dat"), mean_trivec_hist.p4)
    writedlm(joinpath(output, "mean_trivec_hist_all.dat"), mean_trivec_hist.all)

end
