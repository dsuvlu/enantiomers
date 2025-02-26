using CliffordNumbers
using Chemfiles
using LinearAlgebra
using DelimitedFiles
using StaticArrays
using ArgParse

include("/home/dsuvlu/git/enantiomers/src/julia/main.jl")

# set parameters
params = Parameters(3.5, 30.0, 0.87243)

dbins = 0.05

# clifford algebra basis vectors
@basis_vars(VGA3D, Int)

v1 = 1.0 * e1 + 0.0 * e2 + 0.0 * e3
v2 = 0.0 * e1 + 1.0 * e2 + 0.0 * e3
v3 = 0.0 * e1 + 0.0 * e2 + 1.0 * e3

V = v1 ∧ v2 ∧ v3


traj_name = "minimized.pdb"
trajfile = joinpath(directory, traj_name)
#pdbfile = joinpath(directory, pdb_name)
trajectory = Trajectory(trajfile)
#set_topology!(trajectory, pdbfile)
n_frames = Int(length(trajectory))

# grab frame, number of atoms
f = 0
frame = read_step(trajectory, f)
n_atoms = Int(length(frame))

# select all water oxygens
sel_oxy = Selection("name O")
oxy = Int.(evaluate(sel_oxy, frame))
n_o = length(oxy)

# grab current unit cell
ucell = lengths(UnitCell(frame))

# grab current positions, convert to static arrays
pos = positions(frame)
pos_static = Array{SVector{3,Float64}}(undef, n_atoms)
for i in 1:n_atoms
    pos_static[i] = SVector{3, Float64}(pos[:,i])
end

# find neighbors withing dOO distance
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