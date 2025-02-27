using CliffordNumbers
using Chemfiles
using LinearAlgebra
using DelimitedFiles
using StaticArrays
using ArgParse

include("/home/dsuvlu/git/enantiomers/src/julia/main.jl")

delta = 2.5
center = [20.0, 20.0, 20.0]
origin = [0.0, 0.0, 0.0]

# load trajectory, determine number of frames
directory = "/home/dsuvlu/git/enantiomers/tests/"
traj_name = "four_waters.pdb"
trajfile = joinpath(directory, traj_name)
trajectory = Trajectory(trajfile)

shift = [1.0 1.0 1.0;
         0.0 -1.0 -1.0;
         0.0 0.0 1.0]*delta

target = [1.0 0.0 0.0 1.0;
          0.0 -1.0 0.0 -1.0;
          0.0 0.0 1.0 1.0]

shift_and_rotate!(trajectory, shift, center, target)