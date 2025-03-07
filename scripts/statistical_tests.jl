using Statistics
using Printf
using Plots
using BlockingMethod
using HypothesisTests
using Distances
using Serialization
using NPZ


# Reference: Flyvbjerg and Petersen, J. Chem. Phys., 91, 461, 1989

"""
    estimate(x::Vector{Float64})::NTuple{2, Float64}

Performs time series analysis on the `x` column of data outputting
the mean and standard error (variance).

# Flyvbjerb and Petersen 1989 - Abstract
We describe how the true statistical error on an average of correlated data
can be obtained with ease and efficiency by renormalization group method [...]
Reference article https://doi.org/10.1063/1.457480 for more info.
"""
function estimate(x::Vector{Float64})::NTuple{2, Float64}

    n = length(x)
    fn = Float64(n)

    xm = sum(x) / fn
    x2m = sum(y -> y^2, x) / fn
    σ = 0.0

    # if too small exit out
    if (n < 2) @goto fin end

    # number of bins
    nrbin = log(fn) / log(2) - 1.
    @printf("number of blocks = %d\n", nrbin)

    c_max = (x2m - xm^2) / (fn - 1.)
    fac = 1. / sqrt(2. * fn - 2.)

    dc = fac * c_max

    c_old = 0.0

    # calculating sigma
    for i in 1:nrbin
        
        c = 0.0
        n = floor(Int64, n / 2)
        fn = Float64(n)

        for j in 1:n
            if(n == 1)
                continue
            else
                fac = 1. / sqrt(2. * fn - 2.)
                x[j] = (x[2*j] + x[2*j-1]) / 2.
                c += (x[j] - xm)^2 / (fn^2 - fn)
            end
        end

        dc = fac*sqrt(c)
        diff = sqrt(c) - sqrt(c_old)
        @printf("i = %d, xm = %12.5e, cmax = %12.5e,\n", i, xm, sqrt(c_max))

        if (abs(diff) < dc && c > c_max) c_max = c end

        c_old = c

    end

    σ = sqrt(c_max)

    @label fin
    return (xm, σ)
end

using NPZ
n_frames = 1000000
directory = "/home/dsuvlu/git/enantiomers/results/alanine/s/1000ns/"
sdata = zeros(Float64, (10, n_frames))
open(joinpath(directory, "p1_trivec.bin"), "r") do IO
    for f in 1:n_frames
        sdata[:,f] = deserialize(IO)
    end
end
NPZ.npzwrite(joinpath(directory, "p1_trivec.npy"), transpose(sdata))

directory = "/home/dsuvlu/git/enantiomers/results/alanine/r/1000ns/"
rdata = zeros(Float64, (10, n_frames))
open(joinpath(directory, "p1_trivec.bin"), "r") do IO
    for f in 1:n_frames
        rdata[:,f] = deserialize(IO)
    end
end
NPZ.npzwrite(joinpath(directory, "p1_trivec.npy"), transpose(rdata))

estimate(filter(!isnan, sdata[1,:]))
mean(filter(!isnan, sdata[1,:]))
estimate(filter(!isnan, rdata[1,:]))

enddata = 996000

# Example usage:
# Here we generate some synthetic data. In practice, x would be your time series (e.g., measurements of x(t)).
block_vars, block_errs, block_sizes = blocking_analysis(filter(!isnan, sdata[1,:])[1:enddata])

@printf("Blocking analysis results:\n")
for i in 1:length(block_vars)
    @printf("Level %2d: N = %4d, Var(mean) = %12.5e, Error = %12.5e\n",
            i, block_sizes[i], block_vars[i], block_errs[i])
end

# A simple heuristic for the final error estimate is to take the value at the last blocking level,
# where the block size is largest (and the blocks are assumed uncorrelated).
final_error_estimate = block_errs[end]
@printf("\nFinal error estimate on the mean: %12.5e\n", final_error_estimate)

histogram(sdata[1,:], bins=300, label="", xlims=(-0.25, 0.25))
histogram!(rdata[1,:], bins=300, label="")

ks_test = ApproximateTwoSampleKSTest(filter(!isnan, sdata[1,:]), filter(!isnan, rdata[1,:]))


n_frames = 1000000
directory = "/home/dsuvlu/git/enantiomers/results/alanine/s/1000ns/"
sdata = zeros(Float64, (10, n_frames))
sdata = Dict{Int64, Vector{Vector{Union{Float64, HBond}}}}()
open(joinpath(directory, "p1_info.bin"), "r") do IO
    for f in 0:n_frames
        #sdata[:,f] = deserialize(IO)
        tmp = deserialize(IO)
        sdata[f] = tmp[f]
        #println(tmp)
    end
end
NPZ.npzwrite(joinpath(directory, "p1_trivec.npy"), transpose(sdata))