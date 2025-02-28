### check for existence of hbond between i, j pairs of water molecules ###

struct Parameters 
    dOO::Float64 # maximum oxygen-oxygen distance
    angle::Float64 # maximum HOO angle
    dOH::SVector{1, Float64} # OH bond length
    delta_tbins::Float64 # bin width for trivector histograms
    delta_dbins::Float64 # bin width for distance histograms
end

mutable struct Hist1D
    p1::Vector{Float64}
    p2::Vector{Float64}
    p3::Vector{Float64}
    p4::Vector{Float64}
    all::Vector{Float64}
    bins1::Vector{Float64}
    bins2::Vector{Float64}
end

mutable struct Hist2D
    p1::Matrix{Float64}
    p2::Matrix{Float64}
    p3::Matrix{Float64}
    p4::Matrix{Float64}
    all::Matrix{Float64}
    bins1::Vector{Float64}
    bins2::Vector{Float64}
end

struct Neighbors 
    index::Int
    dist::Float64
end

struct HBond
    donor::Int
    acceptor::Int
    hydrogen::Int
    angle::Float64
    dist::Float64
end

const hbond_list = Vector{Union{HBond, Float64}}
mutable struct Chains
    p1::Dict{Int, Vector{hbond_list}}
    p2::Dict{Int, Vector{hbond_list}}
    p3::Dict{Int, Vector{hbond_list}}
    p4::Dict{Int, Vector{hbond_list}}
end


mutable struct Trivecs
    p1::Dict{Int, Vector{Float64}}
    p2::Dict{Int, Vector{Float64}}
    p3::Dict{Int, Vector{Float64}}
    p4::Dict{Int, Vector{Float64}}
end

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--directory"
            help = "Directory path"
            arg_type = String
            required = true
        "--traj"
            help = "File name"
            arg_type = String
            required = true
        "--pdb"
            help = "PDB file name"
            arg_type = String
            required = true
        "--output"
            help = "Output directory"
            arg_type = String
            required = true
        "--solute"
            help = "Solute atoms"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

function initialize_dictionaries(oxy)
    
    nhbs = Dict{Int, Vector{Neighbors}}()
    
    donors = Dict{Int, Vector{HBond}}()
    
    chains = Chains(Dict{Int, Vector{HBond}}(), Dict{Int, Vector{HBond}}(),
        Dict{Int, Vector{HBond}}(), Dict{Int, Vector{HBond}}())
    
    trivecs = Trivecs(Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}())
    
    for o in oxy
        
        nhbs[o] = Vector{Neighbors}()
        donors[o] = Vector{HBond}()
        chains.p1[o] = Vector{HBond}()
        chains.p2[o] = Vector{HBond}()
        chains.p3[o] = Vector{HBond}()
        chains.p4[o] = Vector{HBond}()
        trivecs.p1[o] = Vector{Float64}()
        trivecs.p2[o] = Vector{Float64}()
        trivecs.p3[o] = Vector{Float64}()
        trivecs.p4[o] = Vector{Float64}()
        
    end

    return nhbs, donors, chains, trivecs

end

function clear_values!(oxy, nhbs, donors, chains, trivecs)

    for o in oxy
        
        nhbs[o] = Vector{Neighbors}()
        donors[o] = Vector{HBond}()
        chains.p1[o] = Vector{HBond}()
        chains.p2[o] = Vector{HBond}()
        chains.p3[o] = Vector{HBond}()
        chains.p4[o] = Vector{HBond}()
        trivecs.p1[o] = Vector{Float64}()
        trivecs.p2[o] = Vector{Float64}()
        trivecs.p3[o] = Vector{Float64}()
        trivecs.p4[o] = Vector{Float64}()
        
    end

end

function initialize_histograms(params, solute_flag)

    if solute_flag == "no"

        tbins = collect(-1:params.delta_tbins:1)
        n_tbins = length(tbins) - 1

        cbins = collect(0:1:20)
        n_cbins = length(cbins) - 1

        trivec_hist = Hist1D(zeros(n_tbins), zeros(n_tbins), zeros(n_tbins),
            zeros(n_tbins), zeros(n_tbins), tbins, tbins)

        mean_trivec_hist = Hist1D(zeros(n_tbins), zeros(n_tbins), zeros(n_tbins),
            zeros(n_tbins), zeros(n_tbins), tbins, tbins)

        chain_hist = Hist1D(zeros(n_cbins), zeros(n_cbins), zeros(n_cbins),
            zeros(n_cbins), zeros(n_cbins), cbins, cbins)

    elseif solute_flag == "yes"

        tbins = collect(-1:params.delta_tbins:1)
        n_tbins = length(tbins) - 1

        dbins = collect(0.0:params.delta_dbins:35.0)
        n_dbins = length(dbins) - 1

        cbins = collect(0:1:20)
        n_cbins = length(cbins) - 1

        trivec_hist = Hist2D(zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)),
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), tbins, dbins)
            
        mean_trivec_hist = Hist2D(zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)),
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), tbins, dbins)

        chain_hist = Hist2D(zeros((n_cbins, n_dbins)), zeros((n_cbins, n_dbins)), zeros((n_cbins, n_dbins)),
            zeros((n_cbins, n_dbins)), zeros((n_cbins, n_dbins)), cbins, dbins)

    end

    return trivec_hist, mean_trivec_hist, chain_hist

end

function solute_mass(solute, frame)
    tot_mass = 0.0
    for i in solute
        tot_mass += mass(Atom(frame, i))
    end
    return tot_mass
end

function compute_com_distance(oxy, solute, pos, ucell)
    n_o = length(oxy)
    com = center_of_mass(pos, frame, solute, tot_mass)
    for i in 1:n_o
        com += pos[oxy[i]+1]
    end
    com /= n_o
    return com
end

function center_of_mass(pos, frame, solute, tot_mass)
    com = SVector{3, Float64}(0.0, 0.0, 0.0)
    for i in solute
        a = Atom(frame, i)
        com += mass(a)*pos[i+1]
    end
    com /= tot_mass
    return com
end

function compute_distance(pos_i, pos_j, ucell)
    rij = pos_j - pos_i
    rij -= round.(rij./ucell).*ucell
    return norm(rij)
end

function minimum_distance(pos, ucell, oxy, solute)
    min_dist = Dict{Int, Vector{Float64}}()
    for i in oxy
        for j in solute
            dist = compute_distance(pos[i+1], pos[j+1], ucell)
            if haskey(min_dist, i)
                push!(min_dist[i], dist)
            else
                min_dist[i] = [dist]
            end
        end
        min_dist[i] = [minimum(min_dist[i])]
    end
    return min_dist
end

function find_neighbors!(nhbs, oxy, n_o, pos, params, ucell)
    
    for i in 1:n_o
        o_i = oxy[i]
        
        for j in 1:n_o
            o_j = oxy[j]
            if o_i != o_j
                dist = compute_distance(pos[o_i+1], pos[o_j+1], ucell)
                if dist <= params.dOO
                    if haskey(nhbs, o_i)
                        push!(nhbs[o_i], Neighbors(o_j, dist))
                    else
                        nhbs[o_i] = [Neighbors(o_j, dist)]
                    end
                end
            end
        end

    end

end

function find_neighbors_fast!(nhbs, oxy, n_o, pos, params, ucell)
    
    dOO = params.dOO  # Cache the distance threshold

    @inbounds for i in 1:n_o
        o_i = oxy[i]
        # Loop over unique pairs only: assume symmetric neighbor relationship.
        for j in i+1:n_o
            o_j = oxy[j]
            # Compute the distance between the two oxygen atoms.
            dist = compute_distance(pos[o_i+1], pos[o_j+1], ucell)
            if dist <= dOO
                # Add each as a neighbor of the other.
                push!(nhbs[o_i], Neighbors(o_j, dist))
                push!(nhbs[o_j], Neighbors(o_i, dist))
            end
        end
    end

end

function find_neighbors_faster!(nhbs, oxy, n_o, pos, params, ucell)
    
    dOO = params.dOO  # Cache the distance threshold

    @inbounds for i in 1:n_o
        o_i = oxy[i]
        # Loop over unique pairs only: assume symmetric neighbor relationship.
        for j in i+1:n_o
            o_j = oxy[j]
            # Compute the distance between the two oxygen atoms.
            @fastmath dist = compute_distance(pos[o_i+1], pos[o_j+1], ucell)
            if dist <= dOO
                nhbs[i,j] = dist
                nhbs[j,i] = dist
            end
        end
    end

end

function find_donors!(donors, nhbs, n_o, oxy, frame, params)

    for i in 1:n_o

        o_i = oxy[i]
        h1_i = o_i + 1
        h2_i = o_i + 2

        for j in nhbs[o_i]

            dist = j.dist
            o_j = j.index

            if dist <= params.dOO && o_i != o_j

                # i donor, j acceptor
                angle1 = rad2deg(angle(frame, o_j, o_i, h1_i))
                angle2 = rad2deg(angle(frame, o_j, o_i, h2_i))

                if angle1 < params.angle
                    
                    push!(donors[o_i], HBond(o_i, o_j, h1_i, angle1, dist))
                    
                elseif angle2 < params.angle
                    
                    push!(donors[o_i], HBond(o_i, o_j, h2_i, angle2, dist))
                    
                end
            end
        end

    end

end

function find_donors_chemfiles!(hbonds, n_o, oxy, frame, params)

    for i in 1:n_o

        o_i = oxy[i]
        h1_i = o_i + 1
        h2_i = o_i + 2

        sel_tmp = Selection("name O and (distance(#1, index $o_i) <= 3.5)")
        j = evaluate(sel_tmp, frame)

        for o_j in j

            dist = distance(frame, o_i, o_j)

            if dist <= params.dOO && o_i != o_j

                # i donor, j acceptor
                angle1 = rad2deg(angle(frame, o_j, o_i, h1_i))
                angle2 = rad2deg(angle(frame, o_j, o_i, h2_i))

                if angle1 < params.angle
                    if haskey(hbonds, o_i)
                        push!(hbonds[o_i], HBond(o_i, o_j, h1_i, angle1, dist))
                    else
                        hbonds[o_i] = [HBond(o_i, o_j, h1_i, angle1, dist)]
                    end
                elseif angle2 < params.angle
                    if haskey(hbonds, o_i)
                        push!(hbonds[o_i], HBond(o_i, o_j, h2_i, angle2, dist))
                    else
                        hbonds[o_i] = [HBond(o_i, o_j, h2_i, angle2, dist)]
                    end
                end
            end
        end

    end

end

function adjacency_matrix!(A, donors, n_o, oxy)

    for i in 1:n_o

        o_i = oxy[i]

        if haskey(donors, o_i)

            for hb in donors[o_i]

                o_j = hb.acceptor

                if o_i != o_j

                    j = searchsortedlast(oxy, o_j)
                    A[i, j] = 1
                    A[j, i] = -1

                end

            end

        end

    end

end

function find_hbond_chains!(donors, n_o, oxy, chains)

    for i in 1:n_o
        
        p1 = oxy[i]
            
        for hb1 in donors[p1]
        
            p2 = hb1.acceptor
                        
            for hb2 in donors[p2]
                                
                p3 = hb2.acceptor
                            
                if p3 != p1
                                    
                    for hb3 in donors[p3]
                                                            
                        push!(chains.p1[p1], [hb1, hb2, hb3])                
                        push!(chains.p2[p2], [hb1, hb2, hb3])
                        push!(chains.p3[p3], [hb1, hb2, hb3])
                        p4 = hb3.acceptor                
                        push!(chains.p4[p4], [hb1, hb2, hb3])
                                        
                    end

                end
            end
        
        end
        
    end

end

function trivector(chain, pos, params)
    h1 = chain[1].hydrogen
    h2 = chain[2].hydrogen
    h3 = chain[3].hydrogen

    p1 = chain[1].donor
    p2 = chain[2].donor
    p3 = chain[3].donor

    v1 = (pos[:,h1+1] .- pos[:,p1+1]) / params.dOH
    v2 = (pos[:,h2+1] .- pos[:,p2+1]) / params.dOH
    v3 = (pos[:,h3+1] .- pos[:,p3+1]) / params.dOH

    V1 = v1[1]*e1 + v1[2]*e2 + v1[3]*e3
    V2 = v2[1]*e1 + v2[2]*e2 + v2[3]*e3
    V3 = v3[1]*e1 + v3[2]*e2 + v3[3]*e3

    trivector = V1 ∧ V2 ∧ V3
    bitidx = BitIndex(trivector, 1, 2, 3)
    return trivector[bitidx]
end

function trivector_static(chain, pos, params)
    h1 = chain[1].hydrogen
    h2 = chain[2].hydrogen
    h3 = chain[3].hydrogen

    p1 = chain[1].donor
    p2 = chain[2].donor
    p3 = chain[3].donor

    v1 = (pos[h1+1] - pos[p1+1]) / params.dOH
    v2 = (pos[h2+1] - pos[p2+1]) / params.dOH
    v3 = (pos[h3+1] - pos[p3+1]) / params.dOH

    V1 = v1[1]*e1 + v1[2]*e2 + v1[3]*e3
    V2 = v2[1]*e1 + v2[2]*e2 + v2[3]*e3
    V3 = v3[1]*e1 + v3[2]*e2 + v3[3]*e3

    trivector = V1 ∧ V2 ∧ V3
    bitidx = BitIndex(trivector, 1, 2, 3)
    return trivector[bitidx]
end

function compute_trivectors!(pos, params, chains, trivecs)

    for (p1, p1chains) in chains.p1

        for chain in p1chains

            trivec = trivector_static(chain, pos, params)
            push!(trivecs.p1[p1], trivec)

        end

    end

    for (p2, p2chains) in chains.p2

        for chain in p2chains

            trivec = trivector_static(chain, pos, params)
            push!(trivecs.p2[p2], trivec)

        end

    end

    for (p3, p3chains) in chains.p3

        for chain in p3chains

            trivec = trivector_static(chain, pos, params)
            push!(trivecs.p3[p3], trivec)

        end

    end

    for (p4, p4chains) in chains.p4

        for chain in p4chains

            trivec = trivector_static(chain, pos, params)
            push!(trivecs.p4[p4], trivec)

        end

    end

end

function compute_mean_pi(trivecs)

    mean_trivecs = Dict{Int, Float64}()
    for (key, value) in trivecs
        trivec = mean(value)
        mean_trivecs[key] = trivec
    end

    return mean_trivecs
end

function compute_mean_all(trivecs, oxy)

    mean_trivecs = Dict{Int, Float64}()
    for o_i in oxy
        if haskey(trivecs.p1, o_i) && haskey(trivecs.p2, o_i) && haskey(trivecs.p3, o_i) && haskey(trivecs.p4, o_i)        
            trivec = sum(trivecs.p1[o_i]) + sum(trivecs.p2[o_i]) + sum(trivecs.p3[o_i]) + sum(trivecs.p4[o_i])
            N = length(trivecs.p1[o_i]) + length(trivecs.p2[o_i]) + length(trivecs.p3[o_i]) + length(trivecs.p4[o_i])
            trivec /= N
            mean_trivecs[o_i] = trivec
        end
    end
    
    return mean_trivecs
end

function histogram_trivecs!(counts, data, bins, dbins, min_dist::Dict{Int, Vector{Float64}}=Dict{Int, Vector{Float64}}())

    if typeof(counts) == Vector{Float64}    
    
        if typeof(data) == Dict{Int64, Float64}
            @inbounds for (k, x) in data
                # Only consider x within the histogram range.
                if x >= bins[1] && x <= bins[end]
                    # Find the last index where bins[i] is less than or equal to x.
                    i = searchsortedlast(bins, x)
                    # If x equals the last bin edge, assign it to the final bin.
                    if i == length(bins)
                        i -= 1
                    end
                    counts[i] += 1
                end
            end
        elseif typeof(data) == Dict{Int64, Vector{Float64}}
            @inbounds for (k, x) in data
                for xi in x
                    # Only consider x within the histogram range.
                    if xi >= bins[1] && xi <= bins[end]
                        # Find the last index where bins[i] is less than or equal to x.
                        i = searchsortedlast(bins, xi)
                        # If x equals the last bin edge, assign it to the final bin.
                        if i == length(bins)
                            i -= 1
                        end                        
                        counts[i] += 1
                    end
                end
            end
        end

    elseif typeof(counts) == Matrix{Float64}

        if typeof(data) == Dict{Int64, Float64}
            @inbounds for (k, x) in data
                # Only consider x within the histogram range.
                if x >= bins[1] && x <= bins[end] && min_dist[k][1] >= dbins[1] && min_dist[k][1] <= dbins[end]
                    # Find the last index where bins[i] is less than or equal to x.
                    i = searchsortedlast(bins, x)
                    j = searchsortedlast(dbins, min_dist[k][1])
                    # If x equals the last bin edge, assign it to the final bin.
                    if i == length(bins)
                        i -= 1
                    end
                    if j == length(dbins)
                        j -= 1
                    end
                    counts[i,j] += 1
                end
            end
        elseif typeof(data) == Dict{Int64, Vector{Float64}}
            @inbounds for (k, x) in data
                for xi in x
                    # Only consider x within the histogram range.
                    if xi >= bins[1] && xi <= bins[end] && min_dist[k][1] >= dbins[1] && min_dist[k][1] <= dbins[end]
                        # Find the last index where bins[i] is less than or equal to x.
                        i = searchsortedlast(bins, xi)
                        j = searchsortedlast(dbins, min_dist[k][1])
                        # If x equals the last bin edge, assign it to the final bin.
                        if i == length(bins)
                            i -= 1
                        end
                        if j == length(dbins)
                            j -= 1
                        end
                        counts[i,j] += 1
                    end
                end
            end
        end

    end

end

function histogram_chain_participation!(counts, data, bins, dbins, min_dist::Dict{Int, Vector{Float64}}=Dict{Int, Vector{Float64}}())

    if typeof(counts) == Vector{Float64}    
        @inbounds for (k, x) in data
            n_chains = length(x)
            # Only consider x within the histogram range.
            if n_chains >= bins[1] && n_chains <= bins[end]
                # Find the last index where bins[i] is less than or equal to x.
                i = searchsortedlast(bins, n_chains)
                # If x equals the last bin edge, assign it to the final bin.
                if i == length(bins)
                    i -= 1
                end                        
                counts[i] += 1
            end
        end
    elseif typeof(counts) == Matrix{Float64}        
        @inbounds for (k, x) in data
            n_chains = length(x)
            # Only consider x within the histogram range.
            if n_chains >= bins[1] && n_chains <= bins[end] && min_dist[k][1] >= dbins[1] && min_dist[k][1] <= dbins[end]
                # Find the last index where bins[i] is less than or equal to x.
                i = searchsortedlast(bins, n_chains)
                j = searchsortedlast(dbins, min_dist[k][1])
                # If x equals the last bin edge, assign it to the final bin.
                if i == length(bins)
                    i -= 1
                end
                if j == length(dbins)
                    j -= 1
                end
                counts[i,j] += 1
            end
            
        end
    end

end

function write_chains!(chain_info, f, trivecs, min_dist, chains, params, output)

    for (i, d) in min_dist
        if d[1] <= params.dOO
            if haskey(chains.p1, i)
                for (c, chain) in enumerate(chains.p1[i])
                    oxy1 = chain[1].donor
                    oxy2 = chain[2].donor
                    oxy3 = chain[3].donor
                    oxy4 = chain[3].acceptor
                    trivec = round(trivecs.p1[i][c], digits=4)
                    dist = round(d[1], digits=4)
                    if haskey(chain_info.p1, f)
                        push!(chain_info.p1[f], [trivec, dist, oxy1, oxy2, oxy3, oxy4])
                    else
                        chain_info.p1[f] = [[trivec, dist, oxy1, oxy2, oxy3, oxy4]]
                    end
                    
                end
            end
            if haskey(chains.p2, i)
                for (c, chain) in enumerate(chains.p2[i])
                    oxy1 = chain[1].donor
                    oxy2 = chain[2].donor
                    oxy3 = chain[3].donor
                    oxy4 = chain[3].acceptor
                    trivec = round(trivecs.p2[i][c], digits=4)
                    dist = round(d[1], digits=4)
                    if haskey(chain_info.p2, f)
                        push!(chain_info.p2[f], [trivec, dist, oxy1, oxy2, oxy3, oxy4])
                    else
                        chain_info.p2[f] = [[trivec, dist, oxy1, oxy2, oxy3, oxy4]]
                    end
                end
            end
            if haskey(chains.p3, i)
                for (c, chain) in enumerate(chains.p3[i])
                    oxy1 = chain[1].donor
                    oxy2 = chain[2].donor
                    oxy3 = chain[3].donor
                    oxy4 = chain[3].acceptor
                    trivec = round(trivecs.p3[i][c], digits=4)
                    dist = round(d[1], digits=4)
                    if haskey(chain_info.p3, f)
                        push!(chain_info.p3[f], [trivec, dist, oxy1, oxy2, oxy3, oxy4])
                    else
                        chain_info.p3[f] = [[trivec, dist, oxy1, oxy2, oxy3, oxy4]]
                    end
                end
            end
            if haskey(chains.p4, i)
                for (c, chain) in enumerate(chains.p4[i])
                    oxy1 = chain[1].donor
                    oxy2 = chain[2].donor
                    oxy3 = chain[3].donor
                    oxy4 = chain[3].acceptor
                    trivec = round(trivecs.p4[i][c], digits=4)
                    dist = round(d[1], digits=4)
                    if haskey(chain_info.p4, f)
                        push!(chain_info.p4[f], [trivec, dist, oxy1, oxy2, oxy3, oxy4])
                    else
                        chain_info.p4[f] = [[trivec, dist, oxy1, oxy2, oxy3, oxy4]]
                    end
                end
            end
        end
    end

    open(joinpath(output,"p1_info.bin"), "a") do io
        serialize(io, chain_info.p1)
    end
    open(joinpath(output,"p2_info.bin"), "a") do io
        serialize(io, chain_info.p2)
    end
    open(joinpath(output,"p3_info.bin"), "a") do io
        serialize(io, chain_info.p3)
    end
    open(joinpath(output,"p4_info.bin"), "a") do io
        serialize(io, chain_info.p4)
    end
    
end

function shift_position!(pos, location, i)

    shift = location - pos[:,i]
    pos[:,i] += shift
    pos[:,i+1] += shift
    pos[:,i+2] += shift
    pos[:,i+3] += shift

end

function find_rotation(v, w)
    # Normalize the input vectors.
    a = v / norm(v)
    b = w / norm(w)

    # If the vectors are nearly identical, return the identity.
    if isapprox(a, b; atol=1e-8)
        return Matrix{T}(I, 3, 3)
    end

    # If the vectors are opposite, choose an arbitrary perpendicular axis.
    if isapprox(a, -b; atol=1e-8)
        tmp = abs(a[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        axis = cross(a, tmp)
        axis /= norm(axis)
        # For a rotation by π, the formula simplifies to:
        return -Matrix{T}(I, 3, 3) + 2 * (axis * axis')
    end

    # Compute the rotation axis and angle.
    axis = cross(a, b)
    s = norm(axis)
    axis /= s               # Normalize the axis
    theta = acos(dot(a, b)) # Rotation angle

    # Build the skew-symmetric matrix for the axis.
    K = [  0.0      -axis[3]  axis[2];
         axis[3]      0.0     -axis[1];
        -axis[2]   axis[1]      0.0 ]

    # Rodrigues' rotation formula.
    R = Matrix{Float64}(I, 3, 3) + sin(theta)*K + (1 - cos(theta))*(K*K)
    return R
end


function rotate_to_target!(pos, center, target, i)
    
    pos[:,i+1] -= pos[:,i]
    pos[:,i+2] -= pos[:,i]
    pos[:,i+3] -= pos[:,i]
    pos[:,i] -= pos[:,i]

    source = pos[:,i+1]

    R = find_rotation(source, target)

    pos[:,i+1] = R * pos[:,i+1]
    pos[:,i+2] = R * pos[:,i+2]
    pos[:,i+3] = R * pos[:,i+3]

    pos[:,i] += center
    pos[:,i+1] += center
    pos[:,i+2] += center
    pos[:,i+3] += center

end

function shift_and_rotate!(trajectory, shift, center, target)

    # grab frame, number of atoms
    f = 0
    frame = read_step(trajectory, f)

    # select all water oxygens
    sel_oxy = Selection("name O")
    oxy = Int.(evaluate(sel_oxy, frame))

    # grab current positions
    pos = positions(frame)

    shift_position!(pos, center, oxy[1]+1)
    shift_position!(pos, center.+shift[:,1], oxy[2]+1)
    shift_position!(pos, center.+shift[:,2], oxy[3]+1)
    shift_position!(pos, center.+shift[:,3], oxy[4]+1)

    Trajectory("shifted.pdb", 'w') do traj
        write(traj, frame)
    end

    rotate_to_target!(pos, center, target[:,1], oxy[1]+1)
    rotate_to_target!(pos, center.+shift[:,1], target[:,2], oxy[2]+1)
    rotate_to_target!(pos, center.+shift[:,2], target[:,3], oxy[3]+1)
    rotate_to_target!(pos, center.+shift[:,3], target[:,4], oxy[4]+1)

    Trajectory("rotated.pdb", 'w') do traj
        write(traj, frame)
    end

end