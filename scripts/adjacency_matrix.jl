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
solute_flag = "no"

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
directory = "/home/dsuvlu/git/enantiomers/data/water/"
traj_name = "trajectory.dcd"
pdb_name = "water.pdb"
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

# initialize adjacency Matrix
A = zeros(n_o, n_o)

# select all solute atoms
if solute_flag == "yes"
    sel_sol = Selection("not resname HOH")
    solute = Int.(evaluate(sel_sol, frame))
end

# loop over all frames
#for f in 0:n_frames-1

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

    # find neighbors withing dOO distance
    nhbs = Dict{Int, Vector{Neighbors}}()
    @time find_neighbors!(nhbs, oxy, pos_static, params, ucell)

    # find hydrogen bond donors
    donors = Dict{Int, Vector{HBond}}()
    @time find_donors!(donors, nhbs, n_o, oxy, frame, params)
    #@time find_donors_chemfiles!(donors, n_o, oxy, frame, params)

    @time adjacency_matrix!(A, donors, n_o, oxy)
    
#end
