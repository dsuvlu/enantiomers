using CliffordNumbers
using Chemfiles
using LinearAlgebra
using DelimitedFiles
using StaticArrays
using ArgParse

# load trajectory, determine number of frames
directory = "/home/dsuvlu/git/enantiomers/tests/"
traj_name = "four_waters.pdb"
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

delta = 2.5
center = [20.0, 20.0, 20.0]
origin = [0.0, 0.0, 0.0]

pos = positions(frame)

shift_position!(pos, center, oxy[1]+1)
shift_position!(pos, center.+[delta,0,0], oxy[2]+1)
shift_position!(pos, center.+[delta,delta,0], oxy[3]+1)
shift_position!(pos, center.+[delta,delta,delta], oxy[4]+1)

Trajectory("shifted.pdb", 'w') do traj
    #set_topology!(traj, topology(frame))
    write(traj, frame)
end

rotate_to_target!(pos, center, [1,0,0], oxy[1]+1)
rotate_to_target!(pos, center.+[delta,0,0], [0,1,0], oxy[2]+1)
rotate_to_target!(pos, center.+[delta,delta,0], [0,0,1], oxy[3]+1)
rotate_to_target!(pos, center.+[delta,delta,delta], [1,1,1], oxy[4]+1)

Trajectory("rotated.pdb", 'w') do traj
    #set_topology!(traj, topology(frame))
    write(traj, frame)
end
