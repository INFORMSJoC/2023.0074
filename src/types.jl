using Parameters

@with_kw mutable struct Config
    data::Dict{Symbol, Any} = Dict{Symbol, Any}()
end

# methods for dot notation access (e.g. config.m = 10)
Base.getproperty(c::Config, s::Symbol) = s === :data ? getfield(c, :data) : getfield(c, :data)[s]
Base.setproperty!(c::Config, s::Symbol, v) = s === :data ? setfield!(c, :data, v) : (getfield(c, :data)[s] = v)
Base.propertynames(c::Config) = (:data, keys(getfield(c, :data))...)

# methods for dictionary-style access (e.g. config[:m] = 10)
Base.getindex(c::Config, key) = getfield(c, :data)[Symbol(key)]
Base.setindex!(c::Config, v, key) = (getfield(c, :data)[Symbol(key)] = v)

mutable struct cMCMDProblem
    n_nodes::Int
    n_edges::Int
    n_commodities::Int
    n_days::Int
    n_clusters::Int
    k::Int               # Sparsity param
    A::Array{Float64}   # n_nodes x n_edges matrix
    b::Array{Float64}   # n_commodities x n_nodes vector
    c::Array{Float64}   # n_edges cost vector
    d::Vector{Float64}  # n_edges cost vector
    u::Vector{Float64}  # n_edges capacity vector    

    lb::Vector{Float64} # n_edges lower bound / feasible start for z
    ub::Vector{Float64} # n_edges upper bound for z

    # Some parameters
    γ::Float64                      # regularization
    sampling_rate::Float64          # divide days by s_rate each iteration
    sampling_rate_kelley::Float64   # divide days by s_rate each iteration
    
    # EDGES 
    edge_map::Dict{Int, Tuple{Int, Int, Int}} # TRACK edges before and after sparsifying - mostly needed for interpreting results
    outgoing_edges::Dict{Int, Vector{Int}}    # for e in outgoing_edges[i], e is an outgoing edge from node i
    incoming_edges::Dict{Int, Vector{Int}}    # for e in incoming_edges[i], e is an incoming edge to node i
    old_to_new_map::Dict{Int, Int}            # old_to_new_map[i] = j, means that edge i is now edge j
end


# NOTE: THE OFFSET INCLUDES THE dot(∇obj, -z0) TERM
mutable struct Cut
    obj::Float64
    ∇obj::Array{Float64}
    status::Symbol
    # status::MathOptInterface.TerminationStatusCode
end


mutable struct Dual
    α::Array{Float64}
    λ::Float64
    βl::Array{Float64}
    βu::Array{Float64}
    ρ::Array{Float64}
    w::Array{Float64}
    ofv::Float64
    status::Symbol
end

  
mutable struct kelleyCutInfo
    addRootnodeCuts::Bool
    kelleyPrimalEpochs::Int
    ε::Float64
    method::String
end



mutable struct ClusteringInfo
    cluster_partition::Vector{Vector{Int}}
    clustering::KmeansResult
end


mutable struct PrimalSolution
    support::Vector{Int}
    z::Vector{Float64}
    # t::Float64
    value::Float64
    # offset::Float64
    offset::Array{Float64}
    slope::Array{Float64}
    isbinary::Bool
    method::String
    R_sample::Vector{Int}
    C_samples::Vector{Vector{Int}}
    P::Dict{String, Vector{Int}}
end

mutable struct InfeasibleCut
    z::Vector{Float64}
    R_sample::Vector{Int}
    p::Array{Float64}
    b::Array{Float64}
    infeas_r::Vector{Int}
end

mutable struct DataLogger
    data::Dict{String, Any}
end