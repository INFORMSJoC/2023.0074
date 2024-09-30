println("Load Packages ...")
include("src/nd_core.jl")

function main(args)
    #  SET CONFIGURATION BASED ON CMD LINE ARGS
    config = setWorkspace(args)
    create_objfile(type = config.method*"_main_"*config.kelleyCuts.method*"_root", mode = "create")

    # PROBLEM GENERATION
    println("[Main] TIME LIMIT = $(config.TIME_LIMIT)")    
    time_feas = time()
    Random.seed!(config.randomSeed)
    cmcmd_prob = cMCMDGetProblem(config)

    # SCP METHOD -> SCP METHOD ROOT
    println("[Main] METHOD = $(config.method)")
    println("[Main] KELLEY METHOD = $(config.kelleyCuts.method)")
    time_scp     = @elapsed objval_scp, zopt_scp, model_scp, gap_scp, time_setup_scp, time_opt_scp, confidence_adjusted_bound_gap, objtrue_iters, objbnd_scp = SCPNetworkOpt(cmcmd_prob, config)
    objtrue_scp, _, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt_scp)
    gaptrue_scp = gapcalculate(objtrue_scp, objbnd_scp)
    
    # EXIT AND LOG OVERALL STUFF
    println("[DataNetworkDesign Exit] \n\t\tobjval_scp = $(objval_scp), objtrue_scp = $(objtrue_scp), gap_scp = $(gap_scp), gaptrue_scp = $(gaptrue_scp),\n\t\tconfidence_adjusted_bound_gap = $(confidence_adjusted_bound_gap), nCuts = $(nCuts), γ = $(config.γ), jobid = $(config.jobid), objtrue_iters = $(objtrue_iters)")
    log_scp = log2dict(config, time_scp, time_setup_scp, time_opt_scp, objval_scp, objtrue_scp, objbnd_scp, gap_scp, gaptrue_scp, confidence_adjusted_bound_gap, nCuts, 0, objtrue_iters)
    create_objfile(;mydict = log_scp, type = config.method*"_main_"*config.kelleyCuts.method*"_root", mode = "append")
end


println("START")
main(ARGS)
println("DONE")