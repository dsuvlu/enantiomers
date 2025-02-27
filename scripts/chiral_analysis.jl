using CliffordNumbers
using Chemfiles
using LinearAlgebra
using DelimitedFiles
using BenchmarkTools
using StaticArrays
using ArgParse

include("/home/dsuvlu/git/enantiomers/src/julia/main.jl")

args = parse_command_line()
directory = args["directory"]
traj_name = args["traj"]
pdb_name = args["pdb"]
output = args["output"]
solute_flag = args["solute"]

deltabins = 0.05

# initialize histograms
if solute_flag == "no"

    bins = collect(-1:deltabins:1)
    n_bins = length(bins) - 1
    bins_mid = (bins[1:end-1] + bins[2:end]) / 2

    hist = Hist(zeros(n_bins), zeros(n_bins), zeros(n_bins),
        zeros(n_bins), zeros(n_bins), bins, bins_mid, zeros(length(bins)))

    hist_mean = Hist(zeros(n_bins), zeros(n_bins), zeros(n_bins),
        zeros(n_bins), zeros(n_bins), bins, bins_mid, zeros(length(bins)))

elseif solute_flag == "yes"
    
    dbins = collect(0.0:3.0:33.0)
    n_dbins = length(dbins) - 1

    bins = collect(-1:deltabins:1)
    n_bins = length(bins) - 1
    bins_mid = (bins[1:end-1] + bins[2:end]) / 2

    hist = Hist(zeros((n_bins, n_dbins)), zeros((n_bins, n_dbins)), zeros((n_bins, n_dbins)),
        zeros((n_bins, n_dbins)), zeros((n_bins, n_dbins)), bins, bins_mid, dbins)
        
    hist_mean = Hist(zeros((n_bins, n_dbins)), zeros((n_bins, n_dbins)), zeros((n_bins, n_dbins)),
        zeros((n_bins, n_dbins)), zeros((n_bins, n_dbins)), bins, bins_mid, dbins)

end

# clifford algebra basis vectors
@basis_vars(VGA3D, Int)

# set parameters
params = Parameters(3.5, 30.0, 0.87243)

# load trajectory, determine number of frames
#directory = "/home/dsuvlu/git/enantiomers/data/alanine/"
#traj_name = "r_trajectory.dcd"
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

# select all solute atoms
if solute_flag == "yes"
    sel_sol = Selection("not resname HOH")
    solute = Int.(evaluate(sel_sol, frame))
end

# loop over all frames
for f in 0:n_frames-1

    # grab current frame
    println("Frame: ", f)
    frame = read_step(trajectory, f)
    
    # grab current unit cell
    ucell = lengths(UnitCell(frame))

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

    # find neighbors within dOO distance
    nhbs = Dict{Int, Vector{Neighbors}}()
    @time find_neighbors!(nhbs, oxy, pos_static, params, ucell)

    # find hydrogen bond donors
    donors = Dict{Int, Vector{HBond}}()
    @time find_donors!(donors, nhbs, n_o, oxy, frame, params)
    #@time find_donors_chemfiles!(donors, n_o, oxy, frame, params)

    # find hydrogen bond chains
    chains = Chains(Dict{Int, Array{Vector{HBond}}}(), Dict{Int, Array{Vector{HBond}}}(),
        Dict{Int, Array{Vector{HBond}}}(), Dict{Int, Array{Vector{HBond}}}())
    @time find_hbond_chains!(donors, n_o, oxy, chains)

    # compute trivectors
    trivecs = Trivecs(Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}())
    @time compute_trivectors!(pos, params, chains, trivecs)

    # compute mean trivectors
    p1_mean = compute_mean_pi(trivecs.p1)
    p2_mean = compute_mean_pi(trivecs.p2)
    p3_mean = compute_mean_pi(trivecs.p3)
    p4_mean = compute_mean_pi(trivecs.p4)
    all_mean = compute_mean_all(trivecs, oxy)

    # compute histograms for mean trivectors
    histogram_trivecs!(hist_mean.p1, p1_mean, hist_mean.bins, hist_mean.dbins, min_dist)
    histogram_trivecs!(hist_mean.p2, p2_mean, hist_mean.bins, hist_mean.dbins, min_dist)
    histogram_trivecs!(hist_mean.p3, p3_mean, hist_mean.bins, hist_mean.dbins, min_dist)
    histogram_trivecs!(hist_mean.p4, p4_mean, hist_mean.bins, hist_mean.dbins, min_dist)
    histogram_trivecs!(hist_mean.all, all_mean, hist_mean.bins, hist_mean.dbins, min_dist)

    # compute histograms for trivectors
    histogram_trivecs!(hist.p1, trivecs.p1, hist.bins, hist.dbins, min_dist)
    histogram_trivecs!(hist.p2, trivecs.p2, hist.bins, hist.dbins, min_dist)
    histogram_trivecs!(hist.p3, trivecs.p3, hist.bins, hist.dbins, min_dist)
    histogram_trivecs!(hist.p4, trivecs.p4, hist.bins, hist.dbins, min_dist)
end

if solute_flag == "no"
    # write histograms to file
    histfile = joinpath(output, "hist.dat")
    histmeanfile = joinpath(output, "hist_mean.dat")
    writedlm(histmeanfile, hcat(hist_mean.p1, hist_mean.p2, 
        hist_mean.p3, hist_mean.p4, hist_mean.all))
    writedlm(histfile, hcat(hist.p1, hist.p2, hist.p3, hist.p4))
elseif solute_flag == "yes"
    writedlm(joinpath(output, "hist_mean_p1.dat"), hist_mean.p1)
    writedlm(joinpath(output, "hist_mean_p2.dat"), hist_mean.p2)
    writedlm(joinpath(output, "hist_mean_p3.dat"), hist_mean.p3)
    writedlm(joinpath(output, "hist_mean_p4.dat"), hist_mean.p4)
    writedlm(joinpath(output, "hist_mean_all.dat"), hist_mean.all)
    writedlm(joinpath(output, "hist_p1.dat"), hist.p1)
    writedlm(joinpath(output, "hist_p2.dat"), hist.p2)
    writedlm(joinpath(output, "hist_p3.dat"), hist.p3)
    writedlm(joinpath(output, "hist_p4.dat"), hist.p4)
end

#=
# Plots
using Plots
plot(hist_mean.mbins, hist_mean.p1[:,6])
plot!(hist_mean.mbins, hist_mean.p2)
plot!(hist_mean.mbins, hist_mean.p3)
plot!(hist_mean.mbins, hist_mean.p4)
plot!(hist_mean.mbins, hist_mean.all)

plot(hist.mbins, hist.p1)
plot!(hist.mbins, hist.p2)
plot!(hist.mbins, hist.p3)
plot!(hist.mbins, hist.p4)

histogram(p1_mean)
histogram!(p2_mean)
histogram!(p3_mean)
histogram!(p4_mean)

p1_tmp = Float64[]
for (key, value) in p1_trivecs
    for v in value
        push!(p1_tmp, v)
    end
end
histogram(p1_tmp)
=#