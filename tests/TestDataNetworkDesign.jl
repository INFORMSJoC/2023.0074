include("../src/nd_core.jl")

function compact_summary(results)
    println("\nCompact Summary (values rounded to 3 decimals):")
    println("Method                                        | Obj Value    | True Obj     | Obj Bnd      | Gap     | True Gap | Time   | Edges  |conf_bnd_gap | outter_iters")
    println("----------------------------------------------|--------------|--------------|--------------|---------|----------|--------|--------|-------------|-------------")
    
    function safe_round(x; digits=3)
        if x isa Number && !isnan(x)
            return round(x, digits=digits)
        else
            return x
        end
    end
    
    method_order = sort([key for key in keys(results) if key != "naive"])
    method_order = ["naive", method_order...]
    for key in method_order
        if haskey(results, key)
            value = results[key]
            method = rpad(key, 45)
            objval = lpad(string(safe_round(value["objval"])), 12)
            objtrue = lpad(string(safe_round(value["objtrue"])), 12)
            objbnd = lpad(string(safe_round(value["objbnd"])), 12)
            gap = lpad(string(safe_round(value["gap"])), 7)
            gaptrue = lpad(string(safe_round(value["gaptrue"])), 8)
            time = lpad(string(safe_round(value["time"])), 6)
            edges = lpad(string(safe_round(sum(value["zopt"]))), 6)
            conf_bnd_gap = lpad(string(safe_round(value["confidence_adjusted_bound_gap"])), 11)
            outter_iters = lpad(string(safe_round(value["objtrue_iters"])), 11)
            
            println("$method | $objval | $objtrue | $objbnd | $gap | $gaptrue | $time | $edges | $conf_bnd_gap | $outter_iters")
        end
    end

end


function solution_checks(results)
    println("\nSolution Checks:")
    
    naive_zopt = results["naive"]["zopt"]
    
    # 1. Check if deterministic of each method produces same zopt as naive
    println("\n1. Deterministic solutions matching naive:")
    for (key, value) in results
        if endswith(key, "_det") && key != "naive"
            matches = all(isapprox.(value["zopt"], naive_zopt, atol=1e-3))
            println("  $key: $(matches ? "Matches" : "Differs")")
        end
    end
    
    # 2. Check if stochastic of each method produces same results as deterministic
    println("\n2. Stochastic solutions matching corresponding deterministic:")
    for (key, value) in results
        if endswith(key, "_stoch")
            det_key = replace(key, "_stoch" => "_det")
            if haskey(results, det_key)
                matches = all(isapprox.(value["zopt"], results[det_key]["zopt"], atol=1e-3))
                println("  $key: $(matches ? "Matches" : "Differs")")
            end
        end
    end
end

function compare_solutions(results; do_verify=true)
    compact_summary(results)
    if do_verify
        solution_checks(results)
    end
end

function update_cmcmd_prob(cmcmd_prob, config)
    cmcmd_prob_dict = Dict(
        :sampling_rate => :R_div,
        :sampling_rate_kelley => :R_div_kelley,
        :n_nodes => :m,
        :n_edges => :n,
        :n_days => :nR, 
        :n_clusters => :n_clusters,
        :k => :k,
    )
    for (cmcmd_key, config_key) in cmcmd_prob_dict
        if haskey(config.data, config_key)
            # check if config and cmcmd_prob have differing values for the same key
            cmcmd_value = getproperty(cmcmd_prob, cmcmd_key)
            config_value = config[config_key]
            if cmcmd_value != config_value
                println("[cMCMD Update: config.$(config_key) --> cmcmd_prob.$(cmcmd_key)] Updating $(cmcmd_key) from $(cmcmd_value) to $(config_value)")
                setproperty!(cmcmd_prob, cmcmd_key, config_value)
            end
        end
    end
    return cmcmd_prob
end

function generate_custom_configs(base_config::Dict, param_arrays::Dict=Dict(); differentiating_configs::Vector{<:Dict}=Dict[])
    custom_configs = Dict()

    if !isempty(param_arrays)
        # Generate configs based on parameter arrays
        param_names = collect(keys(param_arrays))
        param_values = collect(values(param_arrays))
        
        for combination in Iterators.product(param_values...)
            config = deepcopy(base_config)
            config_name = "config"
            for (name, value) in zip(param_names, combination)
                config[name] = value
                config_name *= "_$(name)$(value)"
            end
            custom_configs[config_name] = config
        end
    elseif !isempty(differentiating_configs)
        # Generate configs based on differentiating configs
        for (i, diff_config) in enumerate(differentiating_configs)
            config = deepcopy(base_config)
            for (key, value) in diff_config
                config[key] = value
            end
            config_name = get(diff_config, "name", "config_$i")
            custom_configs[config_name] = config
        end
    else
        # If no variations provided, just use the base config
        custom_configs["base_config"] = base_config
    end

    return custom_configs
end


function main(custom_configs; solve_func_used="test_infeas", catch_error=true)
    all_results = Dict()
    
    for (config_name, config) in custom_configs
        println("\n\n\n\n\n\n\n\n\n\n---------------------------------------------------------------------------------------------------------------------")
        println("\nRunning experiment for configuration: $config_name \n")
        println("---------------------------------------------------------------------------------------------------------------------\n")
        results, cmcmd_prob = run_experiment_with_config(config, solve_func_used=solve_func_used, catch_error=catch_error)
        all_results[config_name] = results
        
        println("\nResults for configuration: $config_name")
        compare_solutions(results)
    end
    
    return all_results
end

function run_experiment_with_config(custom_config; solve_func_used="test_infeas", catch_error=true)
    base_config = setWorkspace(ARGS)    
    # Override base config with custom config
    for (key, value) in custom_config
        setproperty!(base_config, Symbol(key), value)
    end
    base_config.n = base_config.m * base_config.m    
    methods = get(custom_config, "methods", [])
    kelley_methods = get(custom_config, "kelley_methods", ["scp_slim"])
    cuts_vals = get(custom_config, "cuts_vals", [0])
    is_deterministic_vals = get(custom_config, "is_deterministic_vals", [true])
    println("[TEST EXP] CONFIG IS: \n\t$(base_config.data)")
    println("[TEST EXP] (Custom) CONFIG WAS: \n\t$custom_config")
    
    # Generate problem
    Random.seed!(base_config.randomSeed)
    
    cmcmd_prob = cMCMDGetProblem(base_config)
    results = Dict()
    # naive solution
    naive_config = deepcopy(base_config)
    naive_config.method = "naive"
    naive_result = solve_problem(cmcmd_prob, naive_config, solve_func_used=solve_func_used, catch_error=catch_error)
    results["naive"] = naive_result
    
    # PRINT
    println("[TEST EXP] RAN NAIVE:")
    println("\t\tγ = $(cmcmd_prob.γ), n_days = $(cmcmd_prob.n_days), is_cmcmd_feasible = $(cMCMD_is_feasible(cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.lb, cmcmd_prob.u))")
    println("\t\tTotal Edges for cmcmd_prob.ub = $(sum(cmcmd_prob.ub)), Total Edges for cmcmd_prob.lb = $(sum(cmcmd_prob.lb))")
    println("\t\tTotal Cost for Naive = $(naive_result["objtrue"]), Total Edges for Naive = $(sum(naive_result["zopt"]))")

    for method in methods
        for kelley_method in kelley_methods
            for cuts in cuts_vals
                for is_deterministic in is_deterministic_vals
                    config = deepcopy(base_config)
                    config.method = method
                    config.method_kelley = (method == "scp_hybrid" ? "scp_fat" : method)
                    config.kelleyCuts = kelleyCutInfo((cuts > 0), cuts, 1e-6, config.method_kelley)
                    config.R_div = is_deterministic ? 1 : 2
                    config.R_div_kelley = is_deterministic ? 1 : 2
                    cmcmd_prob = update_cmcmd_prob(cmcmd_prob, config)
                    println("----------------------------------------------")
                    println("[TEST EXP] Running $(config.method) with $(config.method_kelley) cuts: $(config.kelleyCuts.kelleyPrimalEpochs) and deterministic: $(is_deterministic) [R_div = $(config.R_div), R_div_kelley = $(config.R_div_kelley)]")
                    println("----------------------------------------------")

                    result = solve_problem(cmcmd_prob, config, solve_func_used=solve_func_used, catch_error = catch_error)
                    key = "$(method)_kelley-$(kelley_method)_cuts-$(cuts)_$(is_deterministic ? "det" : "stoch")"
                    results[key] = result
                end
            end
        end
    end
    
    return results, cmcmd_prob
end

function solve_problem(cmcmd_prob, config; solve_func_used="test_infeas", catch_error = true)
    println("\t details: config.method = $(config.method), config.method_kelley = $(config.method_kelley), config.kelleyCuts.kelleyPrimalEpochs = $(config.kelleyCuts.kelleyPrimalEpochs), config.R_div = $(config.R_div), config.R_div_kelley = $(config.R_div_kelley)")
    solve_func = solve_func_used=="test_infeas" ? SCPNetworkOpt : SCPNetworkOpt
    Random.seed!(config.randomSeed)
    objval, zopt, model, gap, confidence_adjusted_bound_gap, objtrue_iters, objbnd, objtrue, gaptrue, time_solve = NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
    if catch_error==false
        time_solve = @elapsed objval, zopt, model, gap, time_setup, time_opt, confidence_adjusted_bound_gap, objtrue_iters, objbnd = solve_func(cmcmd_prob, config)
        objtrue, _, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt)
        gaptrue = gapcalculate(objtrue, objbnd)
    else
        try
            time_solve = @elapsed objval, zopt, model, gap, time_setup, time_opt, confidence_adjusted_bound_gap, objtrue_iters, objbnd = solve_func(cmcmd_prob, config)
            objtrue, _, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt)
            gaptrue = gapcalculate(objtrue, objbnd)
        catch e
            println("***********************************")
            println("ERROR OCCURRED WHILE SOLVING PROBLEM:")
            println("\t details: config.method = $(config.method), config.method_kelley = $(config.method_kelley), config.kelleyCuts.kelleyPrimalEpochs = $(config.kelleyCuts.kelleyPrimalEpochs), config.R_div = $(config.R_div), config.R_div_kelley = $(config.R_div_kelley)")
            println(" Error message: ", e)
            # throw the error if catch_error=false
            println(stacktrace(catch_backtrace()))
            sleep(5)
            objval, zopt, gap, time_solve, confidence_adjusted_bound_gap, objtrue_iters, objbnd = NaN, NaN, NaN, NaN, NaN, NaN, NaN
        end
    end
    return Dict(
        "VERSION" => "$(config.method)_kelley-$(config.kelleyCuts.method)_cuts-$(config.kelleyCuts.kelleyPrimalEpochs)_R$(config.R_div)",
        "objval" => objval,
        "objbnd" => objbnd,
        "objtrue" => objtrue,
        "gap" => gap,
        "gaptrue" => gaptrue,
        "time" => time_solve,
        "zopt" => zopt, 
        "confidence_adjusted_bound_gap" => confidence_adjusted_bound_gap,
        "objtrue_iters" => objtrue_iters
    )

end

#######################################
# SYNTHETIC INSTANCES WITH Z0
#######################################
base_config_synth = Dict(
    "m" => 10,
    "nR" => 5,
    "nC" => 3,
    "TIME_LIMIT" => 60,
    "use_z0" => true,
    "randomSeed" => 42,
    "methods" => ["scp_slim", "scp_fat"],
    "isBenchmark" => false,
    "use_file_demands" => false,
    "kelley_methods" => ["scp_slim", "scp_fat"],
    "cuts_vals" => [0, 20],
    "is_deterministic_vals" => [true, false],
    "useMosek" => true,
    "verbose_logging" => 0,
)


param_arrays = Dict(
    "nC" => [3, 5, 10],
    "verbose_logging" => 0,
    "round_z0" => [true],
)

custom_configs_synth = generate_custom_configs(base_config_synth, param_arrays)

#######################################
# BENCHMARK NO_Z0 INSTANCES WITH GAMMA 
#######################################

base_config_bench = Dict(
    "Rx" => 6,
    "Ry" => 1,
    "nR" => 5,
    "isBenchmark" => true,
    "use_file_demands" => true,
    "TIME_LIMIT" => 60,
    "use_z0" => false,
    "randomSeed" => 43,
    "methods" => ["scp_slim", "scp_hybrid"],
    "cuts_vals" => [0, 20],
    "is_deterministic_vals" => [false],
)

param_arrays = Dict(
    "γ" => [0.1, 0.5, 10, 1000000],
    "verbose_logging" => [1],
    "round_z0" => [true],
)


custom_configs_bench = generate_custom_configs(base_config_bench, param_arrays)

  #####################################
 ############ RUN EXPERIMENTS ########
#####################################
solve_for_bench = false
# xx = Ref{Any}()
if solve_for_bench
    solve_func_used = "test_infeas"
    custom_configs = custom_configs_bench
else
    solve_func_used = "main"
    custom_configs = custom_configs_synth
end

catch_error = false
all_results = main(custom_configs, solve_func_used=solve_func_used, catch_error=catch_error)

all_results = deepcopy(all_results)
for (config_name, results) in all_results
    println("\nOverall results for configuration: $config_name")
    compare_solutions(results, do_verify=false)
end

