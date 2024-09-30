function create_config(parsed_args)
    config = Config()

    # each KEY (command-line argument name) corresponds to a TUPLE of the form (variable_name, processing_function)
    # entries only exist for arguments that need to be PROCESSED OR RENAMED before being passed to the Config constructor
    process_args = Dict(
        "nodes" => ("m", x -> trunc(Int, 10 * x)),
        "days" => ("nR", x -> trunc(Int, x)),
        "commodities" => ("nC", x -> trunc(Int, x)),
        "sample" => ("R_div", x -> trunc(Int, x)),
        "sample_kelley" => ("R_div_kelley", x -> trunc(Int, x)),
        "n_clusters" => ("n_clusters", x -> trunc(Int, x)),
        "time" => ("TIME_LIMIT", identity),
        "seed" => ("randomSeed", identity),
        "gamma" => ("γ", identity),
        "benchmark" => ("isBenchmark", identity),
        "useMosek" => ("useMosek", x -> x == 1 ? true : false),
    )
    
    for (arg, val) in parsed_args
        if haskey(process_args, arg)
            config_name, process_func = process_args[arg]
            config[config_name] = process_func(val)
        else
            config[arg] = val
        end
    end
    
    # POST-PROCESSING
    config.γ = config.γ
    config.n = config.m * config.m
    config.k = config.m <= 20 ? config.n - config.m :
               config.m <= 40 ? round(Int, config.m^2 / 4) :
               round(Int, config.m^2 / 10)
               
    config.rootCuts = config.R_div_kelley == 0 ? 0 : config.rootCuts

    if config.n_clusters_ratio < 0.99
        config.n_clusters = min(10, max(1, ceil(config.n_clusters_ratio * config.nR)))
    end

    if config.R_div == 0 
        config.method = "naive"
        config.rootCuts = 0
    elseif config.R_div == 1 && !config.isBenchmark
        config.R_div_kelley = 1
    end

    config.kelleyCuts = kelleyCutInfo((config.rootCuts > 0), config.rootCuts, 1e-6, config.method_kelley)

    # VER
    config.VER = construct_ver(config)

    # Limit nC up to the number of nodes m. Since each node is responbiel for up to one commodity
    config.nC = min(config.nC, config.m)
    config.corr_pcnt = config.corr_pcnt == 0 ? 0 : config.corr_pcnt # by default, its 0.0 (float), but the files are like: "0-1000", "0.2-1000", ...

    if config.method == "scp_adapt"
        config.rootCuts = 0
        config.kelleyCuts = kelleyCutInfo(false, 0, 1e-6, config.method_kelley)
        config.R_div = 1
        config.is_magnanti_wong_cut = false
    end
    return config
end


function construct_ver(config)
    VER = if config.method in ["scp_slim", "scp_kcut", "scp_fat", "scp_adapt", "scp_hybrid"]
        "$(uppercase(config.ver_prefix))_" * join(uppercase.([split(config.method, "_")[2], split(config.method_kelley, "_")[2]]), "-")
    else
        "$(uppercase(config.ver_prefix))_DETERMINISTIC"
    end

    # Add suffixes based on configuration
    suffixes = String[]
    config.isBenchmark && pushfirst!(suffixes, "R$(config.Rx).$(config.Ry)_r$(config.nR)")
    config.useMosek && push!(suffixes, "MOSEK")
    config.useMosek == false && push!(suffixes, "GUROBI")
    config.isBenchmark && push!(suffixes, "CORR$(config.corr_pcnt)")
    config.use_si_vi && push!(suffixes, "SI_VI")
    config.use_vi_1 && push!(suffixes, "VI_1")
    config.use_vi_2 && push!(suffixes, "VI_2")
    config.use_avg_scenario && push!(suffixes, "AVG_SCEN")
    config.use_partial_cuts && push!(suffixes, "PARTIAL")
    config.is_magnanti_wong_cut && push!(suffixes, "MWC")
    config.use_z0 && push!(suffixes, "Z0")
    config.use_file_demands && push!(suffixes, "FD")

    return join([VER, suffixes...], "_")
end


function setWorkspace(args;step = 10)
    parsed_args = parse_commandline(args)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    println("@@@@@@@@@@@@@@@ Parsed args: @@@@@@@@@@@@@@@")
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    config = create_config(parsed_args)

    println("**********************************************")
    println("*************** Configuration: ***************")
    println("**********************************************")
    for field in propertynames(config)
        if field != :data
            println("  $field: $(config[field])")
        end
    end


    println("useMosek is $(config.useMosek)")
    println("corr_pcnt is $(config.corr_pcnt)")
    println("METHOD = $(config.method)")
    println("VER = $(config.VER)")
    
    return config
end

function parse_commandline(args)
    s = ArgParseSettings()

    @add_arg_table s begin
        "--jobid", "--jid" 
            help = "job id corresponding to this run"
            arg_type = String
            default = "00000"
        "--nodes", "-m"
            help = "option for number of nodes (m=1 means 10 nodes)"
            arg_type = Int
            default = 1
        "--seed", "--randomSeed"
            help = "option for random seed"
            arg_type = Int
            default = 1
        "--commodities"
            help = "option for number of commodities"
            arg_type = Int
            default = 5
        "--Rx"
            help = "option for Rx"
            arg_type = Int
            default = 6
        "--Ry"
            help = "option for Ry"
            arg_type = Int
            default = 5
        "--demand_lb"
            help = "option for demand lower multiplier"
            arg_type = Float64
            default = 0.99
        "--demand_ub"
            help = "option for demand upper multiplier"
            arg_type = Float64
            default = 1.01
        "--gamma"
            help = "option for γ"
            arg_type = Float64
            default = 1.01
        "--sample", "-s"
            help = "option for sampling rate"
            arg_type = Int
            default = 4
        "--sample_kelley", "--s_kelley"
            help = "option for sampling rate"
            arg_type = Int
            default = 4
        "--time", "-t"
            help = "option for timeout threshold"
            arg_type = Int
            default = 3600
        "--days", "-R", "--nR"
            help = "option for R days"
            arg_type = Int
            default = 10
        "--n_clusters"
            help = "option for n_clusters days"
            arg_type = Int
            default = 1
        "--n_clusters_ratio"
            help = "option for n_clusters_ratio -> n_clusters = max(1, ceil(config.nR * ratio))"
            arg_type = Float64
            default = 1.0
        "--rootCuts"
            help = "option for Root Node Cuts"
            arg_type = Int
            default = 0
        "--ver_prefix", "--prefix", "--ver"
            help = "VER prefix"
            arg_type = String
            default = ""
        "--useMosek", "--Mosek"
            help = "flag to use mosek (default) or gurobi for inner problem"
            arg_type = Int
            default = 1
        "--corr_pcnt"
            help = "option for corr_pcnt"
            arg_type = Float64
            default = 0.0
        "--benchmark"
            help = "flag to indicate its a R benchmark"
            action = :store_true
        "--use_z0"
            help = "flag to indicate whether to use a feasible spanning tree as a starting point"
            action = :store_true
        "--use_file_demands"
            help = "flag to indicate whether to use pre-generated file demands or generate new ones on R instances"
            action = :store_true
        "--use_si_vi"
            help = "flag to indicate whether to use strong inequality VIs in root node kelley cuts"
            action = :store_true
        "--use_partial_cuts"
            help = "flag to indicate whether to use partial cuts in infeasibility cuts"
            action = :store_true
        "--use_vi_1"
            help = "flag to indicate whether to use VIs 1 in root node kelley cuts"
            action = :store_true
        "--use_vi_2"
            help = "flag to indicate whether to use VIs 2 in root node kelley cuts"
            action = :store_true
        "--use_avg_scenario"
            help = "flag to indicate whether to use avg scenarion in 2nd stage"
            action = :store_true
        "--is_magnanti_wong_cut"
            help = "flag to indicate whether to use magnanti and wong cuts"
            action = :store_true
        # Methods arguments
        "--method", "--md"
            help = "scp method (e.g. scp_fat or scp_slim)"
            arg_type = String
            default = "scp_slim"
        "--method_kelley", "--md_k" 
            help = "kelley cuts method"
            arg_type = String
            default = "scp_slim"
        "--slim_repeats", "--sr"
            help = "number of repeats for outter"
            arg_type = Int
            default = 4
        "--ntries"
            help = "number of tries to generate feasible synthetic instance"
            arg_type = Int
            default = 50
        "--k_nn"
            help = "k nearest neighbors for clustering"
            arg_type = Int
            default = 6
        "--verbose_logging"
            help = "if 0, run quietly, otherwise print logs"
            arg_type = Int
            default = 0
        "--round_z0"
            help = "round z0 to 0 or 1 during lazy callbacks on integer nodes"
            action = :store_false
    end
    return parse_args(args, s)
end




function gapcalculate(obj, bnd)
    gap = abs(obj - bnd) / obj
    return gap
end

#######################################
# FUNCTIONALITY FOR MODULAR LOGGING 
#######################################
function log2dict(config, time, time_setup, time_opt, objval, objtrue, objbnd, gap, gaptrue, confidence_adjusted_bound_gap, nCuts, unregObj, objtrue_iters)
    log_dict = Dict("VER" => config.VER, "randomSeed" => config.randomSeed, "m" => config.m, "nC" => config.nC, 
                "nR" => config.nR, "R_div_kelley" => config.R_div_kelley, "kelleyPrimalEpochs" => config.kelleyCuts.kelleyPrimalEpochs,
                "time" => time, "time_setup" => time_setup, "time_opt" => time_opt, "objtrue" => objtrue, 
                "gap" => gap, "time" => time, "objtrue" => objtrue, "gap" => gap, 
                "R_div" => config.R_div, "objval" => objval, "gaptrue" => gaptrue, "objbnd"=>objbnd, "nCuts"=>nCuts, "n_clusters"=>config.n_clusters,
                "unregObj"=>unregObj, "gamma"=>config.γ, "jobid"=>config.jobid, "confidence_adjusted_bound_gap"=>confidence_adjusted_bound_gap, "objtrue_iters"=>objtrue_iters)
    return log_dict
end


function create_objfile(;mydict=Dict(), type = "naive", mode = "create")
    if mode == "create"
        if !isfile("data/benchmarks/OBJ_$(uppercase(type)).csv")
            println("CREATING FILE data/benchmarks/OBJ_$(uppercase(type)).csv")
            open("data/benchmarks/OBJ_$(uppercase(type)).csv","w") do fp
                try
                    println(fp,"VER, seed, m, nC, nR, gamma, R_div, R_div_kelley, kelleyPrimalEpochs, n_clusters, nCuts, time, time_setup, time_opt, objval, objtrue, objunreg, objbnd, gap, gaptrue, confidence_adjusted_bound_gap, objtrue_iters")
                catch
                    println(fp,"ERROR")
                end
            end
        else
            println("FILE ALREADY EXISTS (data/benchmarks/OBJ_$(uppercase(type)).csv)")
        end
    elseif mode == "append"
        open("data/benchmarks/OBJ_$(uppercase(type)).csv","a") do fp
            try
                println(fp,"$(mydict["VER"]), $(mydict["randomSeed"]), $(mydict["m"]), $(mydict["nC"]), $(mydict["nR"]), $(mydict["gamma"]), $(mydict["R_div"]),"* 
                    "$(mydict["R_div_kelley"]), $(mydict["kelleyPrimalEpochs"]), $(mydict["n_clusters"]), $(mydict["nCuts"]), $(mydict["time"]), $(mydict["time_setup"]), $(mydict["time_opt"]),"* 
                    "$(mydict["objval"]), $(mydict["objtrue"]), $(mydict["unregObj"]), $(mydict["objbnd"]), $(mydict["gap"]), $(mydict["gaptrue"]), $(mydict["confidence_adjusted_bound_gap"]), $(mydict["objtrue_iters"])")
            catch
                println(fp,"ERROR")
            end
        end
    end
end




function verify_solutions(Zopt_naive, Zopt_solver; tol = 0.05)
    differences = findall(abs.(Zopt_naive .- Zopt_solver) .> tol)
    if length(differences) == 0
        println("The two Solutions are the same")
        second_sum = sum(Zopt_solver)
        println("Total Edges Built: $(second_sum)")
    else
        println("($(length(differences))) Differences spotted : $(differences)")
        first_sum = sum(Zopt_naive)
        second_sum = sum(Zopt_solver)
        println("Edges Built (First) : $(first_sum) \nEdges Built (Second) : $(second_sum)")
    end 
end



function read_scenarios(file_path)
    scenarios = []
    open(file_path) do file
        readline(file) # Skip the first line
        for line in eachline(file)
            parts = split(line, '\t') # Split the line by tab
            probability = parse(Float64, parts[1]) # Parse the probability
            demands = parse.(Float64, parts[2:end]) # Parse the remaining parts as demands
            push!(scenarios, demands) # Add the scenario to the list
        end
    end
    return scenarios
end




function check_vi(zopt, cmcmd_prob)
    edges = findall(cmcmd_prob.ub .> 0)
    m = cmcmd_prob.n_nodes
    nC = cmcmd_prob.n_commodities
    nR = cmcmd_prob.n_days
    b = cmcmd_prob.b
    
    # VI 1: ORIGIN/DESTINATION INEQUALITIES (OD)
    origins = [findall(x -> x < 0, b[:,i,1]) for i in 1:nC]
    destinations = [findall(x -> x > 0, b[:,i,1]) for i in 1:nC]
    is_vi1_feasible = true
    for k in 1:nC
        for origin_node in origins[k]
            lhs = sum(zopt[(origin_node-1)*m + j] for j in 1:m if j != origin_node)
            if lhs < 1
                println("OD Constraint for origin node $(origin_node) in commodity $(k) is violated.")
                is_vi1_feasible = false
            else
                # println("OD Constraint for origin node $(origin_node) in commodity $(k) is satisfied.")
            end
        end

        for destination_node in destinations[k]
            lhs = sum(zopt[(i-1)*m + destination_node] for i in 1:m if (i != destination_node) & ((i-1)*m+destination_node ∈ edges))
            if lhs < 1
                println("OD Constraint for destination node $(destination_node) in commodity $(k) is violated.")
                is_vi1_feasible = false
            else
                # println("OD Constraint for destination node $(destination_node) in commodity $(k) is satisfied.")
            end
        end
    end

    
    # VI 2: NETWORK CONNECTIVITY CUTS
    intermediary_nodes = [i for i in 1:m if all(b[i,k,r] == 0 for k in 1:nC, r in 1:nR)]
    is_vi2_feasible = true
    for i in intermediary_nodes
        incoming_edges = [(j-1)*m + i for j in 1:m if (j != i) & ((j-1)*m + i ∈ edges)]
        outgoing_edges = [(i-1)*m + j for j in 1:m if (j != i) & ((i-1)*m + j ∈ edges)]

        # Check constraint for incoming edges
        for e_minus in incoming_edges
            if zopt[e_minus] > sum(zopt[e_plus] for e_plus in outgoing_edges)
                println("Incoming edge constraint violated for intermediary node i = $(i), edge = $(e_minus)")
                is_vi2_feasible = false
            end
        end
        
        # Check constraint for outgoing edges
        for e_plus in outgoing_edges
            if zopt[e_plus] > sum(zopt[e_minus] for e_minus in incoming_edges)
                println("Outgoing edge constraint violated for intermediary node i = $(i), edge = $(e_plus)")
                is_vi2_feasible = false
            end
        end
    end
    
    return all([is_vi1_feasible, is_vi2_feasible])
end


function sample_dict_elements(P, P_div)
    sampled_keys = sample(collect(keys(P)), Int(ceil(length(P)/P_div)), replace=false)
    return Dict(s => P[s] for s in sampled_keys)
end


function generate_core_point(cmcmd_prob;)
    k, n = cmcmd_prob.k, cmcmd_prob.n_edges
    ϵ = 0.0001
    # z_core .= k/n - ϵ
    z_core = ones(n) .* max(min(k/n - ϵ, 1), ϵ)
    return z_core
end

function verify_core_point_simple(cmcmd_prob, z_core)
    k, lb, ub, n = cmcmd_prob.k, cmcmd_prob.lb, cmcmd_prob.ub, cmcmd_prob.n_edges

    # edges = findall(ub .> 0)
    z0 = lb
    constraints_satisfied = []

    # CONS 1
    ϵ = 0.001
    constr_satisfied = (sum(z_core[e] for e in 1:n) <= k - ϵ)
    push!(constraints_satisfied, constr_satisfied)

    # CONS 2
    open_edges = [e for e in 1:n if lb[e] == 0 && ub[e] == 1]
    z_core_is_frac = [z_core[e] > 0 && z_core[e] < 1 for e in open_edges]
    push!(constraints_satisfied, all(z_core_is_frac))

    # generate 1 point in MP Hull
    if all(constraints_satisfied)
        println("z_core is a Core Point")
        return true
    else
        println("z_core is NOT a Core Point")
        return false
    end
end


function compare_structs(s1, s2)
    comparisons = []
    for field in fieldnames(typeof(s1))
        comparison = (getfield(s1, field) == getfield(s2, field))
        push!(comparisons, comparison)
        if !comparison
            println("different field = $field")
            println("s1 = $(getfield(s1, field))")
            println("s2 = $(getfield(s2, field))")
        end
    end
    
    if all(comparisons)
        println("Structs are the same")
    else
        println("Structs are different")
        println("Number of equal fields = $(sum(comparisons))")
        println("Number of different fields = $(length(comparisons) - sum(comparisons))")
    end
    return comparisons
end


function safe_round(z0; digits = 5)
    z0_round = round.(z0; digits = digits)
    if any(z0_round .!= 0 .&& z0_round .!= 1)
        divergent_vals = findall(z0_round .!= 0 .&& z0_round .!= 1)
        @warn "z0_round is not binary:\nDivergent values at indices $(findall(z0_round .!= 0 .&& z0_round .!= 1)), with values $(z0_round[findall(z0_round .!= 0 .&& z0_round .!= 1)])"
    end
    return z0_round
end



function create_edge_mappings(m, ub)
    available_edges = findall(ub .> 0)
    N_available = length(available_edges)
    println("Available edges: ", available_edges)

    edge_map = Dict{Int, Tuple{Int, Int, Int}}()
    outgoing_edges = Dict{Int, Vector{Int}}()
    incoming_edges = Dict{Int, Vector{Int}}()
    old_to_new_map = Dict{Int, Int}()

    for (new_index, old_index) in enumerate(available_edges)
        start_node = div(old_index - 1, m) + 1
        end_node = mod(old_index - 1, m) + 1
        edge_map[new_index] = (start_node, end_node, old_index)
        old_to_new_map[old_index] = new_index
        println("New index $new_index corresponds to old index $old_index ($(start_node) -> $(end_node))")
        start_node == end_node && @warn "Self-loop detected at edge $old_index ($start_node -> $end_node))"
        if !haskey(outgoing_edges, start_node)
            outgoing_edges[start_node] = []
        end
        push!(outgoing_edges[start_node], new_index)
        
        if !haskey(incoming_edges, end_node)
            incoming_edges[end_node] = []
        end
        push!(incoming_edges[end_node], new_index)
    end

    return available_edges, edge_map, outgoing_edges, incoming_edges, old_to_new_map
end


function calc_mem(nodes)
    memory = Int(round(nodes * 2/3))
    cores = Int(min(round(memory /4 + 1), 35))
    println("supercloud RAM requested: $(cores * 4)")
    return cores * 4
end

function ij2n_index(i,j,m)
    return ((i-1)*m+1:i*m)[j]    
end

function n2ij_index(n,m)
    return div(n,m, RoundUp), mod(n-1,m)+1
end
  
function r_line2nums(line)
    return [parse(Int, m.match) for m in eachmatch(r"\d+", line)]
end


