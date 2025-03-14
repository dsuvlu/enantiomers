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

mutable struct Trivectors
    p1::Union{Float64, Vector{Float64}}
    p2::Union{Float64, Vector{Float64}}
    p3::Union{Float64, Vector{Float64}}
    p4::Union{Float64, Vector{Float64}}
    all::Union{Float64, Vector{Float64}}
end

struct Neighbors 
    index::Int
    dist::Float64
end

struct HBSites
    index::Int
    type::String
    donor::Bool
    acceptor::Bool
    hydrogen::Int
end

struct HBond
    donor::Int
    hydrogen::Int
    acceptor::Int
    angle::Float64
    dist::Float64
end

mutable struct Chains
    o1::Int
    o2::Int
    o3::Int
    o4::Int
    dist::Float64
    solute::Int
    trivec::Float64
end

mutable struct Dist 
    dist::Float64
    sol::Int
end

const chain_list = Union{Vector{HBond}, Chains}
mutable struct ChainList
    p1::Dict{Int, Vector{chain_list}}
    p2::Dict{Int, Vector{chain_list}}
    p3::Dict{Int, Vector{chain_list}}
    p4::Dict{Int, Vector{chain_list}}
end

const trivec_list = Union{Float64, Vector{Float64}}
mutable struct TrivecList
    p1::Dict{Int, trivec_list}
    p2::Dict{Int, trivec_list}
    p3::Dict{Int, trivec_list}
    p4::Dict{Int, trivec_list}
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
            arg_type = Bool
            required = true
    end

    return parse_args(s)
end

function initialize_dictionaries(oxy)
    
    nhbs = Dict{Int, Vector{Neighbors}}()
    
    donors = Dict{Int, Vector{HBond}}()
    
    chains = ChainList(Dict{Int, Vector{HBond}}(), Dict{Int, Vector{HBond}}(),
        Dict{Int, Vector{HBond}}(), Dict{Int, Vector{HBond}}())
    
    trivecs = TrivecList(Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}(),
        Dict{Int, Vector{Float64}}(), Dict{Int, Vector{Float64}}())

    min_dist = Dict{Int, Dist}()
    
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
        min_dist[o] = Dist(0.0, 0)
        
    end

    return nhbs, donors, chains, trivecs, min_dist

end

function clear_values!(oxy, nhbs, donors, chains, trivecs, min_dist)

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
        min_dist[o] = Dist(0.0, 0)
        
    end

end

function initialize_histograms(params, solute_flag)

    if !solute_flag

        tbins = collect(-1:params.delta_tbins:1)
        n_tbins = length(tbins) - 1

        cbins = collect(0:1:20)
        n_cbins = length(cbins) - 1

        trivec_hist_total = Hist1D(zeros(n_tbins), 
            zeros(n_tbins), zeros(n_tbins),
            zeros(n_tbins), zeros(n_tbins), tbins, tbins)

        trivec_hist_frame = Hist1D(zeros(n_tbins), 
            zeros(n_tbins), zeros(n_tbins),
            zeros(n_tbins), zeros(n_tbins), tbins, tbins)

        mean_trivec_hist_total = Hist1D(zeros(n_tbins), 
            zeros(n_tbins), zeros(n_tbins),
            zeros(n_tbins), zeros(n_tbins), tbins, tbins)

        mean_trivec_hist_frame = Hist1D(zeros(n_tbins), 
            zeros(n_tbins), zeros(n_tbins),
            zeros(n_tbins), zeros(n_tbins), tbins, tbins)

        chain_hist = Hist1D(zeros(n_cbins), 
            zeros(n_cbins), zeros(n_cbins),
            zeros(n_cbins), zeros(n_cbins), cbins, cbins)

    elseif solute_flag

        tbins = collect(-1:params.delta_tbins:1)
        n_tbins = length(tbins) - 1

        dbins = collect(0.0:params.delta_dbins:35.0)
        n_dbins = length(dbins) - 1

        cbins = collect(0:1:20)
        n_cbins = length(cbins) - 1

        trivec_hist_total = Hist2D(zeros((n_tbins, n_dbins)), 
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)),
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), tbins, dbins)

        trivec_hist_frame = Hist2D(zeros((n_tbins, n_dbins)), 
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)),
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), tbins, dbins)
            
        mean_trivec_hist_total = Hist2D(zeros((n_tbins, n_dbins)), 
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)),
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), tbins, dbins)
        
        mean_trivec_hist_frame = Hist2D(zeros((n_tbins, n_dbins)), 
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)),
            zeros((n_tbins, n_dbins)), zeros((n_tbins, n_dbins)), tbins, dbins)

        chain_hist = Hist2D(zeros((n_cbins, n_dbins)), 
            zeros((n_cbins, n_dbins)), zeros((n_cbins, n_dbins)),
            zeros((n_cbins, n_dbins)), zeros((n_cbins, n_dbins)), cbins, dbins)

    end

    return trivec_hist_total, trivec_hist_frame,
        mean_trivec_hist_total, mean_trivec_hist_frame, chain_hist

end

function initialize_statistics(params, solute_flag)

    if !solute_flag
    
        trivec_sum_frame = Trivectors(0.0, 0.0, 0.0, 0.0, 0.0)
        trivec_sum_total = Trivectors(0.0, 0.0, 0.0, 0.0, 0.0)

        mean_trivec_sum_frame = Trivectors(0.0, 0.0, 0.0, 0.0, 0.0)
        mean_trivec_sum_total = Trivectors(0.0, 0.0, 0.0, 0.0, 0.0)
    
    elseif solute_flag
    
        dbins = collect(0.0:params.delta_dbins:35.0)
        n_dbins = length(dbins) - 1
    
        trivec_sum_frame = Trivectors(zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins))
        trivec_sum_total = Trivectors(zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins))

        mean_trivec_sum_frame = Trivectors(zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins))
        mean_trivec_sum_total = Trivectors(zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins), zeros(n_dbins))

    end

    return trivec_sum_frame, trivec_sum_total, mean_trivec_sum_frame, mean_trivec_sum_total

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

function minimum_distance!(min_dist, pos, ucell, oxy, solute, sol_dist)
    
    for i in oxy
        for j in solute

            sol_dist[j] = compute_distance(pos[i+1], pos[j+1], ucell)
            
        end
        min_pair = findmin(sol_dist)
        min_dist[i] = Dist(min_pair[1], min_pair[2])
    end
    return min_dist
end

function find_solute_hbonds!(sol_hbonds, sol_info, min_dist, params)

    for (k,v) in min_dist
        if v.dist <= params.dOO
            
            sol_i = v.sol
            
            if sol_info[sol_i].acceptor
                o_j = k
                h1_j = k + 1
                h2_j = k + 2

                # sol acceptor, water donor
                angle1 = rad2deg(angle(frame, sol_i, o_j, h1_j))
                angle2 = rad2deg(angle(frame, sol_i, o_j, h2_j))

                if angle1 < params.angle
                    push!(sol_hbonds[sol_i], HBond(o_j, h1_j, sol_i, angle1, v.dist))
                elseif angle2 < params.angle
                    push!(sol_hbonds[sol_i], HBond(o_j, h2_j, sol_i, angle2, v.dist))
                end

            end

            if sol_info[sol_i].donor
                o_j = k
                sol_h = sol_info[sol_i].hydrogen

                # sol donor, water acceptor
                angle0 = rad2deg(angle(frame, o_j, sol_i, sol_h))

                if angle0 < params.angle
                    push!(sol_hbonds[sol_i], HBond(sol_i, sol_h, o_j, angle0, v.dist))
                end

            end

        end
    end
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
                    
                    push!(donors[o_i], HBond(o_i, h1_i, o_j, angle1, dist))
                    
                elseif angle2 < params.angle
                    
                    push!(donors[o_i], HBond(o_i, h2_i, o_j, angle2, dist))
                    
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

function trivector_static(chain, pos, params, e1, e2, e3)
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

function compute_trivectors!(pos, params, chains, trivecs, e1, e2, e3)

    for (p1, p1chains) in chains.p1

        for chain in p1chains

            trivec = trivector_static(chain, pos, params, e1, e2, e3)
            push!(trivecs.p1[p1], trivec)

        end

    end

    for (p2, p2chains) in chains.p2

        for chain in p2chains

            trivec = trivector_static(chain, pos, params, e1, e2, e3)
            push!(trivecs.p2[p2], trivec)

        end

    end

    for (p3, p3chains) in chains.p3

        for chain in p3chains

            trivec = trivector_static(chain, pos, params, e1, e2, e3)
            push!(trivecs.p3[p3], trivec)

        end

    end

    for (p4, p4chains) in chains.p4

        for chain in p4chains

            trivec = trivector_static(chain, pos, params, e1, e2, e3)
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

function reset_statistics!(per_frame, solute_flag)

    if typeof(per_frame.p1) == Float64
        per_frame.p1 = 0.0
        per_frame.p2 = 0.0
        per_frame.p3 = 0.0
        per_frame.p4 = 0.0
        per_frame.all = 0.0
    else
        per_frame.p1 .= 0.0
        per_frame.p2 .= 0.0
        per_frame.p3 .= 0.0
        per_frame.p4 .= 0.0
        per_frame.all .= 0.0
    end

    return per_frame

end

function update_totals!(total, per_frame)

    if typeof(total.p1) == Float64
        total.p1 += per_frame.p1
        total.p2 += per_frame.p2
        total.p3 += per_frame.p3
        total.p4 += per_frame.p4
        total.all += per_frame.all
    else
        total.p1 .+= per_frame.p1
        total.p2 .+= per_frame.p2
        total.p3 .+= per_frame.p3
        total.p4 .+= per_frame.p4
        total.all .+= per_frame.all
    end
    
    return total

end

function histogram_trivecs!(data, sum, counts, bins, dbins, min_dist, 
    solute_flag, mean_flag)

    if !solute_flag
    
        if mean_flag
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
                    sum += x
                end
            end
        elseif !mean_flag
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
                        sum += xi
                    end
                end
            end
        end

    elseif solute_flag

        if mean_flag
            @inbounds for (k, x) in data

                dist = min_dist[k].dist

                # Only consider x within the histogram range.
                if x >= bins[1] && x <= bins[end] && dist >= dbins[1] && dist <= dbins[end]
                    # Find the last index where bins[i] is less than or equal to x.
                    i = searchsortedlast(bins, x)
                    j = searchsortedlast(dbins, dist)
                    # If x equals the last bin edge, assign it to the final bin.
                    if i == length(bins)
                        i -= 1
                    end
                    if j == length(dbins)
                        j -= 1
                    end
                    counts[i,j] += 1
                    sum[j] += x
                end
            end
        elseif !mean_flag
            @inbounds for (k, x) in data

                dist = min_dist[k].dist

                for xi in x
                    # Only consider x within the histogram range.
                    if xi >= bins[1] && xi <= bins[end] && dist >= dbins[1] && dist <= dbins[end]
                        # Find the last index where bins[i] is less than or equal to x.
                        i = searchsortedlast(bins, xi)
                        j = searchsortedlast(dbins, dist)
                        # If x equals the last bin edge, assign it to the final bin.
                        if i == length(bins)
                            i -= 1
                        end
                        if j == length(dbins)
                            j -= 1
                        end
                        counts[i,j] += 1
                        sum[j] += xi
                    end
                end
            end
        end

    end

end

function histogram_chain_participation!(data, counts, bins, dbins, min_dist, solute_flag)

    if !solute_flag
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
    elseif solute_flag
        @inbounds for (k, x) in data
            n_chains = length(x)
            dist = min_dist[k].dist
            # Only consider x within the histogram range.
            if n_chains >= bins[1] && n_chains <= bins[end] && dist >= dbins[1] && dist <= dbins[end]
                # Find the last index where bins[i] is less than or equal to x.
                i = searchsortedlast(bins, n_chains)
                j = searchsortedlast(dbins, dist)
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

    for (i, v) in min_dist
        if v.dist <= params.dOO
            
            for (c, chain) in enumerate(chains.p1[i])
                oxy1 = chain[1].donor
                oxy2 = chain[2].donor
                oxy3 = chain[3].donor
                oxy4 = chain[3].acceptor
                trivec = trivecs.p1[i][c]
                dist = v.dist
                sol = v.sol                
                if haskey(chain_info.p1, f)
                    push!(chain_info.p1[f], Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec))
                else
                    chain_info.p1[f] = [Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec)]
                end
                
            end
            
            
            for (c, chain) in enumerate(chains.p2[i])
                oxy1 = chain[1].donor
                oxy2 = chain[2].donor
                oxy3 = chain[3].donor
                oxy4 = chain[3].acceptor
                trivec = trivecs.p2[i][c]
                dist = v.dist
                sol = v.sol 
                if haskey(chain_info.p2, f)
                    push!(chain_info.p2[f], Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec))
                else
                    chain_info.p2[f] = [Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec)]
                end
            
            end
            
            
            for (c, chain) in enumerate(chains.p3[i])
                oxy1 = chain[1].donor
                oxy2 = chain[2].donor
                oxy3 = chain[3].donor
                oxy4 = chain[3].acceptor
                trivec = trivecs.p3[i][c]
                dist = v.dist
                sol = v.sol 
                if haskey(chain_info.p3, f)
                    push!(chain_info.p3[f], Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec))
                else
                    chain_info.p3[f] = [Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec)]
                end
                
            end
            
            
            for (c, chain) in enumerate(chains.p4[i])
                oxy1 = chain[1].donor
                oxy2 = chain[2].donor
                oxy3 = chain[3].donor
                oxy4 = chain[3].acceptor
                trivec = trivecs.p4[i][c]
                dist = v.dist
                sol = v.sol 
                if haskey(chain_info.p4, f)
                    push!(chain_info.p4[f], Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec))
                else
                    chain_info.p4[f] = [Chains(oxy1, oxy2, oxy3, oxy4, dist, sol, trivec)]
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

function write_solute_hbonds!(sol_hbonds, f, output)
    open(joinpath(output,"sol_hbonds.bin"), "a") do io
        serialize(io, [f sol_hbonds])
    end
end

function write_trivec_frame!(trivec_sum_frame, trivec_hist_frame, 
    mean_trivec_sum_frame, mean_trivec_hist_frame, output)

    open(joinpath(output,"p1_trivec.bin"), "a") do io
        serialize(io, trivec_sum_frame.p1./vec(sum(trivec_hist_frame.p1, dims=1)))
    end
    open(joinpath(output,"p2_trivec.bin"), "a") do io
        serialize(io, trivec_sum_frame.p2./vec(sum(trivec_hist_frame.p2, dims=1)))
    end
    open(joinpath(output,"p3_trivec.bin"), "a") do io
        serialize(io, trivec_sum_frame.p3./vec(sum(trivec_hist_frame.p3, dims=1)))
    end
    open(joinpath(output,"p4_trivec.bin"), "a") do io
        serialize(io, trivec_sum_frame.p4./vec(sum(trivec_hist_frame.p4, dims=1)))
    end

    open(joinpath(output,"p1_mean_trivec.bin"), "a") do io
        serialize(io, mean_trivec_sum_frame.p1./vec(sum(mean_trivec_hist_frame.p1, dims=1)))
    end
    open(joinpath(output,"p2_mean_trivec.bin"), "a") do io
        serialize(io, mean_trivec_sum_frame.p2./vec(sum(mean_trivec_hist_frame.p2, dims=1)))
    end
    open(joinpath(output,"p3_mean_trivec.bin"), "a") do io
        serialize(io, mean_trivec_sum_frame.p3./vec(sum(mean_trivec_hist_frame.p3, dims=1)))
    end
    open(joinpath(output,"p4_mean_trivec.bin"), "a") do io
        serialize(io, mean_trivec_sum_frame.p4./vec(sum(mean_trivec_hist_frame.p4, dims=1)))
    end

end

function solute_donors_acceptors(frame, solute)

    sol_info = Dict{Int, HBSites}()

    top = Topology(frame)
    bond = Int.(bonds(top))
    sel_all_solute = Selection("not resname HOH")
    all_solute = Int.(Chemfiles.evaluate(sel_all_solute, frame))
    
    for col in eachcol(bond)

        a = col[1]
        b = col[2]
        a_name = name(Atom(frame, a))
        a_type = type(Atom(frame, a))
        b_name = name(Atom(frame, b))
        b_type = type(Atom(frame, b))
        if a in solute && b_type == "H"
            println("Bond between $a_name, $a_type, $a and $b_name, $b_type, $b")
            if a_type == "N"
                sol_info[a] = HBSites(a, a_type, true, false, b)
            elseif a_type == "O"
                sol_info[a] = HBSites(a, a_type, true, true, b)
            end
        end
    end

    for k in solute
        if !haskey(sol_info, k)
            k_type = type(Atom(frame, k))
            sol_info[k] = HBSites(k, k_type, false, true, 0)
        end
    end

    return sol_info

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