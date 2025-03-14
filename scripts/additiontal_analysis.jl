using Chemfiles # for reading in the trajectory
using LinearAlgebra # for computing the norm of vectors
using StaticArrays # for static arrays which are faster for certain operations (e.g. distance calculations)
using Serialization # for saving data structures to binary files
using Plots
using ArgParse
using NPZ
using Statistics
using BlockingMethod

# include functions
include("/home/dsuvlu/git/enantiomers/src/julia/main.jl")

directory = "/home/dsuvlu/git/enantiomers/data/2-butanol/"
r_output = "/home/dsuvlu/git/enantiomers/results/2-butanol/r/1000ns/"
s_output = "/home/dsuvlu/git/enantiomers/results/2-butanol/s/1000ns/"

r_chain_info = ChainList(Dict{Int, Vector{Chains}}(), Dict{Int, Vector{Chains}}(),
        Dict{Int, Vector{Chains}}(), Dict{Int, Vector{Chains}}())

r_sol_hbonds = Dict{Int, Vector{HBond}}()

s_chain_info = ChainList(Dict{Int, Vector{Chains}}(), Dict{Int, Vector{Chains}}(),
        Dict{Int, Vector{Chains}}(), Dict{Int, Vector{Chains}}())

s_sol_hbonds = Dict{Int, Vector{HBond}}()

n_frames = 1000000

for f in 0:n_frames
    r_chain_info.p1[f] = Vector{Chains}()
    r_chain_info.p2[f] = Vector{Chains}()
    r_chain_info.p3[f] = Vector{Chains}()
    r_chain_info.p4[f] = Vector{Chains}()

    s_chain_info.p1[f] = Vector{Chains}()
    s_chain_info.p2[f] = Vector{Chains}()
    s_chain_info.p3[f] = Vector{Chains}()
    s_chain_info.p4[f] = Vector{Chains}()
end

file = "p1_info.bin"
npyfile = "p1_trivec.npy"
r_chi_p1 = Float64[]
open(joinpath(r_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
            
            f = first(keys(tmp))
            #println(f)
            r_chain_info.p1[f] = tmp[f]
            for chain in r_chain_info.p1[f]
                push!(r_chi_p1, chain.trivec)
            end

        end
    end
end
mean(r_chi_p1)
NPZ.npzwrite(joinpath(r_output, npyfile), r_chi_p1)

s_chi_p1 = Float64[]
open(joinpath(s_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
           
            f = first(keys(tmp))
            #println(f)
            s_chain_info.p1[f] = tmp[f]
            for chain in s_chain_info.p1[f]
                push!(s_chi_p1, chain.trivec)
            end

        end
    end
end
mean(s_chi_p1)
NPZ.npzwrite(joinpath(s_output, npyfile), s_chi_p1)


file = "p2_info.bin"
npyfile = "p2_trivec.npy"
r_chi_p2 = Float64[]
open(joinpath(r_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
            
            f = first(keys(tmp))
            #println(f)
            r_chain_info.p2[f] = tmp[f]
            for chain in r_chain_info.p2[f]
                push!(r_chi_p2, chain.trivec)
            end

        end
    end
end
mean(r_chi_p2)
NPZ.npzwrite(joinpath(r_output, npyfile), r_chi_p2)

s_chi_p2 = Float64[]
open(joinpath(s_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
           
            f = first(keys(tmp))
            #println(f)
            s_chain_info.p2[f] = tmp[f]
            for chain in s_chain_info.p2[f]
                push!(s_chi_p2, chain.trivec)
            end

        end
    end
end
mean(s_chi_p2)
NPZ.npzwrite(joinpath(s_output, npyfile), s_chi_p2)



file = "p3_info.bin"
npyfile = "p3_trivec.npy"
r_chi_p3 = Float64[]
open(joinpath(r_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
            
            f = first(keys(tmp))
            #println(f)
            r_chain_info.p3[f] = tmp[f]
            for chain in r_chain_info.p3[f]
                push!(r_chi_p3, chain.trivec)
            end

        end
    end
end
mean(r_chi_p3)
NPZ.npzwrite(joinpath(r_output, npyfile), r_chi_p3)

s_chi_p3 = Float64[]
open(joinpath(s_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
           
            f = first(keys(tmp))
            #println(f)
            s_chain_info.p3[f] = tmp[f]
            for chain in s_chain_info.p3[f]
                push!(s_chi_p3, chain.trivec)
            end

        end
    end
end
mean(s_chi_p3)
NPZ.npzwrite(joinpath(s_output, npyfile), s_chi_p3)



file = "p4_info.bin"
npyfile = "p4_trivec.npy"
r_chi_p4 = Float64[]
open(joinpath(r_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
            
            f = first(keys(tmp))
            #println(f)
            r_chain_info.p4[f] = tmp[f]
            for chain in r_chain_info.p4[f]
                push!(r_chi_p4, chain.trivec)
            end

        end
    end
end
mean(r_chi_p4)
NPZ.npzwrite(joinpath(r_output, npyfile), r_chi_p4)

s_chi_p4 = Float64[]
open(joinpath(s_output, file), "r") do io
    while !eof(io)
        tmp = deserialize(io)
        if !isempty(tmp)
           
            f = first(keys(tmp))
            #println(f)
            s_chain_info.p4[f] = tmp[f]
            for chain in s_chain_info.p4[f]
                push!(s_chi_p4, chain.trivec)
            end

        end
    end
end
mean(s_chi_p4)
NPZ.npzwrite(joinpath(s_output, npyfile), s_chi_p4)

histogram(r_chi_p1, bins=101, normalize=:pdf, label="r-p1", xlabel="χ", ylabel="Frequency", title="Chi distribution")
histogram!(s_chi_p1, bins=101, normalize=:pdf, label="s-p1", xlabel="χ", ylabel="Frequency", title="Chi distribution")

histogram!(r_chi_p2, bins=101, normalize=:pdf, label="r-p2", xlabel="χ", ylabel="Frequency", title="Chi distribution")
histogram!(s_chi_p2, bins=101, normalize=:pdf, label="s-p2", xlabel="χ", ylabel="Frequency", title="Chi distribution")

histogram!(r_chi_p3, bins=101, normalize=:pdf, label="r-p3", xlabel="χ", ylabel="Frequency", title="Chi distribution")
histogram!(s_chi_p3, bins=101, normalize=:pdf, label="s-p3", xlabel="χ", ylabel="Frequency", title="Chi distribution")

histogram!(r_chi_p4, bins=101, normalize=:pdf, label="r-p4", xlabel="χ", ylabel="Frequency", title="Chi distribution")
histogram!(s_chi_p4, bins=101, normalize=:pdf, label="s-p4", xlabel="χ", ylabel="Frequency", title="Chi distribution")

savefig("2-butanol_chi.pdf")

file = "sol_hbonds.bin"
r_hbonds = []
open(joinpath(r_output, file), "r") do io
    while !eof(io)
        push!(r_hbonds, deserialize(io))
    end
end

r_nhbonds = Int[]
for hbonds in r_hbonds
    push!(r_nhbonds, length(hbonds[2][3]))
end
mean(r_nhbonds)
estimate(Float64.(r_nhbonds))

s_hbonds = []
open(joinpath(s_output, file), "r") do io
    while !eof(io)
        push!(s_hbonds, deserialize(io))
    end
end

s_nhbonds = Int[]
for hbonds in s_hbonds
    push!(s_nhbonds, length(hbonds[2][3]))
end
mean(s_nhbonds)
estimate(Float64.(s_nhbonds))

histogram(r_nhbonds, bins=10, normalize=:pdf, label="r", xlabel="Number of solute-solvent hydrogen bonds", ylabel="Frequency", title="Number of solute-solvent hydrogen bonds")
histogram!(s_nhbonds, bins=10, normalize=:pdf, label="s", xlabel="Number of solute-solvent hydrogen bonds", ylabel="Frequency", title="Number of solute-solvent hydrogen bonds")

