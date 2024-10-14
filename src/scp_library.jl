# SCP SUITE
function SCPNetworkOpt(cmcmd_prob, config)
    start_time, time_setup, time_remaining_setup, time_opt = time(), 0, 0, 0
    if config.method == "naive" # Naive Solver
        objval, zopt, model, _, gap = cMCMDNaiveSolver(cmcmd_prob, config)
        return objval, zopt, model, gap, 0, (time() - start_time), 0, 0, objective_bound(model)
    end
    
    method, kelleyCuts, use_vi_2, use_vi_1, use_si_vi, TIME_LIMIT,  is_magnanti_wong_cut, use_avg_scenario =
                                config.method, config.kelleyCuts, config.use_vi_2, config.use_vi_1, config.use_si_vi, 
                                config.TIME_LIMIT, config.is_magnanti_wong_cut, config.use_avg_scenario
    A, b, c, d, u, k, γ, lb, ub, m, n, nC, nR, R_div, n_clusters, R_div_kelley =
                cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.c, cmcmd_prob.d, cmcmd_prob.u, cmcmd_prob.k, cmcmd_prob.γ, cmcmd_prob.lb,
                cmcmd_prob.ub, cmcmd_prob.n_nodes, cmcmd_prob.n_edges, cmcmd_prob.n_commodities, cmcmd_prob.n_days,
                cmcmd_prob.sampling_rate, cmcmd_prob.n_clusters, cmcmd_prob.sampling_rate_kelley

    # SETUP TRACKING VARIABLES
    confidence_adjusted_bound_gap = Inf
    z_last = zeros(n)
    gap_last = 1
    objval_last = 1
    bnd_last = -1
    iter_count = 0
    objtrue_iters = 0 
    clusterInfo = cluster_days(cmcmd_prob, config=config)
    cluster_partition = clusterInfo.cluster_partition
    global nCuts = 0 # for logging purposes only
    global z_core = generate_core_point(cmcmd_prob)
    global P, b_P
    println("[SCPNetworkOpt] Configuration:\n\t\t\tMethod: $(method), \n\t\t\tRoot Node Cuts: $(kelleyCuts), 
            \t\tSampling: $(Int(ceil(nR/R_div)))/$(nR) Days at every iteration (R_div = $(R_div), nR = $(nR)), $(Int(ceil(nR/R_div_kelley)))/$(nR) Days at every Root Node Cut iteration (R_div_kelley = $(R_div_kelley), nR = $(nR)
            \t\tUse VI 1: $(use_vi_1)\n\t\t\tUse VI 2: $(use_vi_2)\n\t\t\tUse SI VI: $(use_si_vi)\n\t\t\tUse Avg Scenario: $(use_avg_scenario)\n\t\t\tTime Limit: $(TIME_LIMIT)\n\t\t\tMagnanti-Wong Cut: $(is_magnanti_wong_cut)\n\t\t\tScenario Partition: $(cluster_partition)")
    
    
    objs = zeros(nR)
    ∇objs = zeros(nR, n)
    obj = 0
    ∇obj = zeros(n)
    if method == "scp_kcut"
        objs = zeros(n_clusters)
        ∇objs = zeros(n_clusters, n)
    end


    # INITIALIZE MODEL
    m1 = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(m1, "OutputFlag", 1) 
    set_optimizer_attribute(m1, "TimeLimit", TIME_LIMIT)
    set_optimizer_attribute(m1, "LazyConstraints", 1) 
    set_optimizer_attribute(m1, "MIPGap", 0.01)

    @variable(m1, z[e in 1:n], Bin)

    # Add constraints
    z0 = lb
    if config.use_z0
        println("--- [SCP] using z0 as lower bound ---")
        @constraint(m1,[e in 1:n], z[e] .>= lb[e])
        for e in 1:n 
            JuMP.set_start_value(z[e], lb[e])
        end
    else
        println("--- [SCP] using 0 as lower bound ---")
        @constraint(m1,[e in 1:n], z[e] .>= 0)
        for e in 1:n 
            JuMP.set_start_value(z[e], ub[e])
        end    
        z0 = ub
    end
    @constraint(m1,[e in 1:n], z[e] .<= ub[e])
    @constraint(m1, sum(z[e] for e in 1:n) <= k)

    # use_avg_scenario = true
    if use_avg_scenario
        println("\n\n\n [AVERAGE SCENARIO] \n\n\n")
        # x ∈ R^(|E| × K)
        @variable(m1, x_avg[e in 1:n,1:nC]>=0)
        # y ∈ R^(|E|)
        @variable(m1, y_avg[e in 1:n])
        # average across R dimension, and drop the R dimension
        b_avg = dropdims(mean(b, dims=3), dims=3)
        # y_avg = ∑x_avg
        @constraint(m1, [e in 1:n], y_avg[e] == sum(x_avg[e,c] for c = 1:nC))
        # Ax_avg = b
        @constraint(m1, [i=1:m, c=1:nC], dot(A[i,1:n], x_avg[:,c]) == b_avg[i,c])
        # y_avg ≤ u
        @constraint(m1, [r=1:nR], y_avg[:] .<= u[1:n])
        # y_avg ≤ M ⋅ z
        @constraint(m1, [e in 1:n], y_avg[e] <= u[e] * z[e])    
    end
    
    # ADD VALID INEQUALITIES IF APPLICABLE
    add_vi_constraints(m1, cmcmd_prob, z, use_vi_1=use_vi_1, use_vi_2=use_vi_2, use_nci=true)

    # PERFORM 1 DETERMINISTIC CUT IN THE BEGINNING
    obj0 = NaN
    println("[Deterministic Cut @ Beginning] method = $(method)")
    add_infeasibility_cut_det = false
    if method == "scp_slim"
        @variable(m1,t >= 0)
        @objective(m1, Min, (dot(c[1:n],z)) + t)
        obj0, ∇obj0 , feas_status = cMCMDCutting_plane(cmcmd_prob, config, z0)
        add_infeasibility_cut_det = any(isnan, obj0)
        if !add_infeasibility_cut_det # optimality cut
            @constraint(m1, t >= obj0 + dot(∇obj0, z - z0[1:n]))
        end
    elseif method == "scp_fat" || method == "scp_hybrid"
        @variable(m1,t[1:nR] >= 0)
        @objective(m1, Min, (dot(c[1:n],z)) + sum(t))
        obj0, ∇obj0  = cMCMDCutting_plane(cmcmd_prob, config, z0, returnObjs=true)
        add_infeasibility_cut_det = any(isnan, obj0)
        if !add_infeasibility_cut_det # optimality cut
            @constraint(m1, [r=1:nR], t[r] >= obj0[r] + dot(∇obj0[r,:], z - z0[1:n]))
        end
    elseif method == "scp_kcut"
        obj0 = zeros(n_clusters)
        ∇obj0 = zeros(n_clusters, n)
        @variable(m1,t[1:n_clusters] >= 0)
        @objective(m1, Min, (dot(c[1:n],z)) + sum(t))

        obj0_fat, ∇obj0_fat = cMCMDCutting_plane(cmcmd_prob, config, ub, returnObjs=true)
        add_infeasibility_cut_det = any(isnan, obj0_fat)
        if !add_infeasibility_cut_det # optimality cut
            for cluster = 1:n_clusters
                cluster_indices = cluster_partition[cluster]
                obj0[cluster] = sum(obj0_fat[cluster_indices])
                ∇obj0[cluster, :] = sum(∇obj0_fat[cluster_indices, :], dims=1)
                @constraint(m1, t[cluster] >= obj0[cluster] + dot(∇obj0[cluster, :], z - ub[1:n]))
            end
        end
    elseif method == "scp_adapt"
        @variable(m1,t[1:nR] >= 0)
        @objective(m1, Min, (dot(c[1:n],z)) + sum(t))
        obj0, ∇obj0  = cMCMDCutting_plane(cmcmd_prob, config, z0, returnObjs=true)

        @constraint(m1, [r=1:nR], t[r] >= obj0[r] + dot(∇obj0[r,:], z - z0[1:n]))
        P = Dict("0" => collect(1:nR)) #Partition
        b_P = Dict(s => mean(b[:,:,P[s]],dims=3) for s in keys(P)) #Average b across all scenarions within each subset of P
    end

    if add_infeasibility_cut_det
        println("\t\t[Deterministic Cut @ Beginning] Adding infeasibility cut")
        p_certificate, β_certificate, infeas_r = get_infeas_certificate(cmcmd_prob, z0)
        if !(any(isnan, p_certificate) & any(isnan, β_certificate))
            for r in infeas_r
                r_tilde = r
                if (any(β_certificate[:,r] .> 0) || any(p_certificate[:,:,r] .!= 0))
                    println("\t\t\t Adding certificate constraint for r = $(r) (r_tilde = $(r_tilde))")
                    @constraint(m1, sum(sum(p_certificate[:,k,r] .* b[:,k,r_tilde]) for k in 1:nC) - sum(z[e] * u[e] * β_certificate[e, r]  for e in 1:n) <= 0)
                end
            end
        end
        @constraint(m1, sum((1- z0[e]) * z[e] for e in 1:n) >= 1)
        @constraint(m1, sum(z0[e] * (1 - z[e]) for e in 1:n) + sum((1- z0[e]) * z[e] for e in 1:n) >= 1)      
    end
    
    ######################################### 
    ############### KELLEY ################## 
    ######################################### 
    if kelleyCuts.addRootnodeCuts
        global z_kelley_relax,  cutPoolKelley, cutPoolInfeasKelley, kelleyPath, lb_kelley, ub_kelley, kelley_time = kelleyPrimal(cmcmd_prob, config, method = kelleyCuts.method, stabilizationPoint = copy(lb), kelleyPrimalEpochs = kelleyCuts.kelleyPrimalEpochs, ε = kelleyCuts.ε, clusterInfo = clusterInfo);
    end

    if kelleyCuts.addRootnodeCuts
        if length(cutPoolKelley) > 0
            # Enable possibility of different method for root node cuts and main body cuts
            println("[Kelley Cuts] Adding $(length(cutPoolKelley)) Root Node Cuts\n\t\t Kelley Method : $(kelleyCuts.method)\n\t\t Main Method : $(method)")
            if (kelleyCuts.method == "scp_fat") && (method == "scp_slim")
                @variable(m1,tau[1:nR] >= 0)
                @constraint(m1, t >= sum(tau))
            elseif (kelleyCuts.method == "scp_fat") && (method == "scp_kcut")
                @variable(m1,tau[1:nR] >= 0)
                for c in 1:n_clusters
                    @constraint(m1, t[c] >= sum(tau[i] for i in cluster_partition[c]))
                end
            elseif (kelleyCuts.method == "scp_kcut") && (method == "scp_slim")
                @variable(m1,tau[1:n_clusters] >= 0)
                @constraint(m1, t >= sum(tau))
            end

            for c in cutPoolKelley
                # ROOT NODE FAT
                if c.method == "scp_fat" && method == "scp_fat"
                    for r in c.R_sample
                        @constraint(m1, t[r] >= c.offset[r] + dot(c.slope[r,:], z))
                    end
                elseif c.method == "scp_fat" && method == "scp_slim"
                    for r in c.R_sample
                        @constraint(m1, tau[r] >= c.offset[r] + dot(c.slope[r,:], z))
                    end
                elseif c.method == "scp_fat" && method == "scp_kcut"
                    for r in c.R_sample
                        @constraint(m1, tau[r] >= c.offset[r] + dot(c.slope[r,1:n], z))
                    end
                # ROOT NODE HYBRID
                elseif c.method == "scp_fat" && method == "scp_hybrid" # There is no hybrid method for root node, since it is equivalent to fat
                    for r in c.R_sample
                        @constraint(m1, t[r] >= c.offset[r] + dot(c.slope[r,:], z))
                    end
                    excluded_indices = setdiff(1:size(c.slope, 1), c.R_sample)
                    total_sum = sum([c.offset[r] + dot(c.slope[r,:], z) for r in excluded_indices], init=0)
                    @constraint(m1, sum(t[excluded_indices]) >= total_sum)
                # ROOT NODE SLIM
                elseif c.method == "scp_slim" && method == "scp_slim"
                    @constraint(m1, t >= c.offset[1] + dot(c.slope[:], z))
                elseif c.method == "scp_slim" && method == "scp_fat"
                    @constraint(m1, sum(t) >= c.offset[1] + dot(c.slope[:], z))
                elseif c.method == "scp_slim" && method == "scp_kcut" # Technically, same as above, but more clear to have it separately CONSOLIDATE
                    @constraint(m1, sum(t) >= c.offset[1] + dot(c.slope[:], z))
                elseif c.method == "scp_slim" && method == "scp_hybrid"
                    println(" [KELLEY] --- ADDING HYBRID SLIM CUTS")
                    @constraint(m1, sum(t) >= c.offset[1] + dot(c.slope[1:n], z))
                # ROOT NODE KCUT
                elseif c.method == "scp_kcut" && method == "scp_kcut"
                    C_active = findall(x -> length(x) > 0, c.C_samples) # Clusters that actually appear in the sample
                    for cluster in C_active
                        @constraint(m1, t[cluster] >= c.offset[cluster] + dot(c.slope[cluster,:], z))
                    end
                elseif c.method == "scp_kcut" && method == "scp_fat"
                    C_active = findall(x -> length(x) > 0, c.C_samples) # Clusters that actually appear in the sample
                    for cluster in C_active
                        @constraint(m1, sum(t[c] for c in cluster_partition[cluster]) >= c.offset[cluster] + dot(c.slope[cluster,:], z))
                    end
                elseif c.method == "scp_kcut" && method == "scp_slim"
                    C_active = findall(x -> length(x) > 0, c.C_samples) # Clusters that actually appear in the sample
                    for cluster in C_active #1
                        @constraint(m1, tau[cluster] >= c.offset[cluster] + dot(c.slope[cluster,:], z))
                    end
                end
            end
        end
        if length(cutPoolInfeasKelley) > 0
            println("\t[Kelley Cuts] Adding $(length(cutPoolInfeasKelley)) Root Node Infeasibility Cuts")
            for c in cutPoolInfeasKelley
                for r in 1:size(c.b)[2]
                    if (any(c.b[:,r] .> 0) || any(c.p[:,:,r] .!= 0))
                        r_tilde = c.R_sample[r]
                        println(" \t[Kelley]  Adding Infeasibility Cut for r = $(r) (r_tilde = $(r_tilde))")
                        @constraint(m1, sum(c.p[:,:,r] .* b[:,:,r_tilde]) <= sum(z[e] * u[e] * c.b[e, r] for e in 1:n))
                    end
                end
            end
        end
    end


    ##################################################################################################################
    # Outer approximation method for Convex Integer Optimization (CIO)
    function Newcut(cb, cb_where)
        if cb_where == GRB_CB_MIPSOL
            global cutPoolLazy_added
            if !cutPoolLazy_added
                if length(cutPoolLazy) > 0
                    println("[Newcut] Adding $(length(cutPoolLazy)) optimality lazy cuts from previous outer iterations")
                    for lazyCut in cutPoolLazy
                        if lazyCut.method == "scp_slim" && method == "scp_slim"
                            con = @build_constraint(t >= lazyCut.offset[1] + dot(lazyCut.slope, z))
                            MOI.submit(m1, MOI.LazyConstraint(cb), con)
                        elseif lazyCut.method == "scp_fat" && method == "scp_fat"
                            for r in lazyCut.R_sample
                                con = @build_constraint(t[r] >= lazyCut.offset[r] + dot(lazyCut.slope[r,:], z))
                                MOI.submit(m1, MOI.LazyConstraint(cb), con)
                            end
                        elseif lazyCut.method == "scp_hybrid" && method == "scp_hybrid"
                            for r in lazyCut.R_sample
                                con = @build_constraint(t[r] >= lazyCut.offset[r] + dot(lazyCut.slope[r,:], z))
                                MOI.submit(m1, MOI.LazyConstraint(cb), con)
                            end
                            excluded_indices = setdiff(1:size(lazyCut.slope, 1), lazyCut.R_sample)
                            total_sum = sum([lazyCut.offset[r] + dot(lazyCut.slope[r,:], z) for r in excluded_indices], init=0)
                            con = @build_constraint(sum(t[excluded_indices]) >= total_sum)
                            MOI.submit(m1, MOI.LazyConstraint(cb), con)
                        elseif lazyCut.method == "scp_kcut" && method == "scp_kcut"
                            C_active = findall(x -> length(x) > 0, lazyCut.C_samples) # Clusters that actually appear in the sample
                            for cluster in C_active
                                con = @build_constraint(t[cluster] >= lazyCut.offset[cluster] + dot(lazyCut.slope[cluster,:], z))
                                MOI.submit(m1, MOI.LazyConstraint(cb), con)
                            end
                        elseif lazyCut.method == "scp_adapt" && method == "scp_adapt"
                            P_sample = lazyCut.P
                            for (idx, s) in enumerate(keys(P_sample))
                                indices_s = P_sample[s]
                                con = @build_constraint(sum(t[indices_s])/length(indices_s) >= lazyCut.offset[idx] + dot(lazyCut.slope[idx,:], z))
                                MOI.submit(m1, MOI.LazyConstraint(cb), con)
                            end
                        end
                    end
                end
                if length(cutPoolInfeas) > 0
                    println("[Newcut] Adding $(length(cutPoolInfeas)) infeasibility lazy cuts from previous outer iterations")
                    for c in cutPoolInfeas
                        for r in 1:size(c.b)[2]
                            if (any(c.b[:,r] .> 0) || any(c.p[:,:,r] .!= 0))
                                r_tilde = c.R_sample[r]
                                con = @build_constraint(sum(c.p[:,:,r] .* b[:,:,r_tilde]) <= sum(z[e] * u[e] * c.b[e, r] for e in 1:n))
                                MOI.submit(m1, MOI.LazyConstraint(cb), con)
                            end
                        end
                    end
                end
                cutPoolLazy_added = true
            end
            # PUSH STUFF
            Gurobi.load_callback_variable_primal(cb, cb_where)
            z_cur = zeros(length(z0))
            z_cur[1:n] = [callback_value(cb, z[e]) for e in 1:n]
            config.round_z0 && (z_cur = safe_round.(z_cur, digits=3))
            primal_bound = Ref{Cdouble}()
            dual_bound = Ref{Cdouble}()
            GRBcbget(cb, cb_where, GRB_CB_MIPSOL_OBJBST, primal_bound)
            GRBcbget(cb, cb_where, GRB_CB_MIPSOL_OBJBND, dual_bound)
            objbst = copy(primal_bound[])   # Incumbent objective 
            objbnd = copy(dual_bound[])     # Objective Bound
            mip_gap = min(gapcalculate(objbst, objbnd), 1) 

            status = MOI.get(JuMP.backend(m1), MOI.CallbackNodeStatus(cb))
            if status == MOI.CALLBACK_NODE_STATUS_INTEGER
                is_zcur_feasible = cMCMD_is_feasible(A, b, z_cur, u)
                println(" [Newcut] Status : $(status) --- Saving z_last\n\t\tis_zcur_feasible = $(is_zcur_feasible) --- ($(sum(z_cur))/$(length(z_cur)))")
                if is_zcur_feasible
                    z_last, gap_last, objval_last, bnd_last = copy(z_cur), copy(mip_gap), copy(objbst), max(objbnd, bnd_last)
                end
            end

            # SAMPLING STUFF
            global R_sample = sample(1:nR, Int(ceil(nR/R_div)), replace=false)
            R_sample = sort(R_sample)
            println("\t\t\tR_sample : $(sort(R_sample))")
        
            # GENERATE CUT
            add_infeasibility_cut = false
            if method == "scp_slim"
                nCuts = nCuts + 1
                if is_magnanti_wong_cut
                    z_core = 0.5 .* (z_cur .+ z_core)
                    if !cMCMD_is_feasible(A, b, z_core, u)
                        z_core = z_cur
                    end
                    println(" [M&W from Newcut] is_magnanti_wong_cut = $(is_magnanti_wong_cut)")
                end
                println(" [Newcut] Calculating Cut")
                obj, ∇obj, feas_status = cMCMDCutting_plane(cmcmd_prob, config, z_cur, R_sample = R_sample, is_magnanti_wong_cut = is_magnanti_wong_cut, z0_magnanti = z_core)
                add_infeasibility_cut = any(isnan, obj)
                # if add_infeasibility_cut and config.use_partial_cuts 
                if add_infeasibility_cut && config.use_partial_cuts 
                    # Adding optimality cuts for the other scenarios
                    R_sample_partial = setdiff(R_sample, R_sample[infeas_r]) # feasible scenarios
                    if length(R_sample_partial) >= 1
                        println("\t\t[use_partial_cuts] Adding $(length(R_sample_partial)) optimality cuts for infeasible z")
                        obj_partial, ∇obj_partial, feas_status = cMCMDCutting_plane(cmcmd_prob, config, z_cur, R_sample = R_sample_partial, is_magnanti_wong_cut = is_magnanti_wong_cut, z0_magnanti = z_core)
                        if isnan(obj_partial)
                            throw(ArgumentError("\t\t\t[use_partial_cuts] NaN obj_partial"))
                        else
                            con = @build_constraint(t >= obj_partial + dot(∇obj_partial, z - z_cur[1:n]))
                            MOI.submit(m1, MOI.LazyConstraint(cb), con)
                            # Store Lazy Single-Cut in Cut Pool
                            s0 = PrimalSolution( findall(z_cur .> 0), z_cur, obj_partial + dot(c,z_cur), [obj_partial - dot(∇obj_partial, z_cur[1:n])], ∇obj_partial, false, method, R_sample, [], Dict())
                            global countLazy, cutPoolPartial
                            push!(cutPoolPartial, s0)
                        end
                    end
                end
                
                if !add_infeasibility_cut
                    con = @build_constraint(t >= obj + dot(∇obj, z - z_cur[1:n]))
                    MOI.submit(m1, MOI.LazyConstraint(cb), con)
                    println("\t[Newcut] Adding Slim Cut")
                    # Store Lazy Single-Cut in Cut Pool
                    s0 = PrimalSolution( findall(z_cur .> 0), z_cur, obj + dot(c,z_cur), [obj - dot(∇obj, z_cur[1:n])], ∇obj, false, method, R_sample, [], Dict())
                    global countLazy, cutPoolLazy
                    push!(cutPoolLazy, s0)
                    countLazy = countLazy + 1
                    countLazy % 10 == 0 ? println("countlazy = $(countLazy)") : nothing
                end
            elseif method == "scp_fat"
                nCuts = nCuts + length(R_sample)
                objs, ∇objs  = cMCMDCutting_plane(cmcmd_prob, config, z_cur, R_sample = R_sample, returnObjs = true)
                add_infeasibility_cut = any(isnan, objs)
                if !add_infeasibility_cut
                    for r in R_sample
                        con = @build_constraint(t[r] >= objs[r] + dot(∇objs[r,:], z - z_cur[1:n]))
                        MOI.submit(m1, MOI.LazyConstraint(cb), con)
                    end
                    # Store Lazy Multi-Cut in Cut Pool
                    global countLazy, cutPoolLazy
                    countLazy = countLazy + length(R_sample)
                    s0 = PrimalSolution( findall(z_cur .> 0), z_cur, sum(objs[r] for r in 1:nR) + dot(c,z_cur), [objs[r] - dot(∇objs[r,:], z_cur[1:n]) for r in 1:nR], ∇objs, false, method, R_sample, [], Dict())
                    push!(cutPoolLazy, s0)
                end
            elseif method == "scp_hybrid"
                nCuts = nCuts + 1
                objs, ∇objs = cMCMDCutting_plane(cmcmd_prob, config, z_cur, R_sample = R_sample, returnObjs=true)
                add_infeasibility_cut = any(isnan, objs)
                if add_infeasibility_cut && config.use_partial_cuts
                    # Add partial optimality cuts for the other scenarios
                    R_sample_partial = setdiff(R_sample, R_sample[infeas_r])
                    if length( R_sample_partial ) >= 1
                        println("Adding some optimality cuts for infeasible z's: ", length( R_sample_partial) , " cuts")
                        objs_partial, ∇objs_partial = cMCMDCutting_plane(cmcmd_prob, config, z_cur, R_sample = R_sample_partial, returnObjs=true)
                        for r in R_sample_partial
                            con = @build_constraint(t[r] >= objs_partial[r] + dot(∇objs_partial[r,:], z - z_cur[1:n]))
                            MOI.submit(m1, MOI.LazyConstraint(cb), con)
                        end
                        excluded_indices = setdiff(1:size(∇objs_partial, 1), R_sample_partial)
                        total_sum = sum([objs_partial[r] + dot(∇objs_partial[r,:], z - z_cur[1:n]) for r in excluded_indices], init=0)
                        con = @build_constraint(sum(t[excluded_indices]) >= total_sum)
                        MOI.submit(m1, MOI.LazyConstraint(cb), con)
                        # Store Lazy Multi-Cut in Cut Pool
                        global countLazy, cutPoolPartial
                        countLazy = countLazy + length(R_sample_partial)
                        s0 = PrimalSolution( findall(z_cur .> 0), z_cur, sum(objs_partial[r] for r in 1:nR) + dot(c,z_cur), [objs_partial[r] - dot(∇objs_partial[r,:], z_cur[1:n]) for r in 1:nR], ∇objs_partial, false, method, R_sample_partial, [], Dict())
                        push!(cutPoolPartial, s0)
                    end
                end
                if !add_infeasibility_cut
                    println(" [Newcut] --- ADDING HYBRID CUT")
                    nCuts = nCuts + length(R_sample) + 1
                    for r in R_sample
                        con = @build_constraint(t[r] >= objs[r] + dot(∇objs[r,:], z - z_cur[1:n]))
                        MOI.submit(m1, MOI.LazyConstraint(cb), con)
                    end
                    excluded_indices = setdiff(1:size(∇objs, 1), R_sample)
                    total_sum = sum([objs[r] + dot(∇objs[r,:], z - z_cur[1:n]) for r in excluded_indices], init=0)
                    con = @build_constraint(sum(t[excluded_indices]) >= total_sum)
                    MOI.submit(m1, MOI.LazyConstraint(cb), con)
                    # Store Lazy Multi-Cut in Cut Pool
                    global countLazy, cutPoolLazy
                    countLazy = countLazy + length(R_sample)
                    s0 = PrimalSolution( findall(z_cur .> 0), z_cur, sum(objs[r] for r in 1:nR) + dot(c,z_cur), [objs[r] - dot(∇objs[r,:], z_cur[1:n]) for r in 1:nR], ∇objs, false, method, R_sample, [], Dict())
                    push!(cutPoolLazy, s0)
                end
            elseif method == "scp_kcut"
                global C_samples = [findall(y -> y in R_sample, x) for x in cluster_partition]
                global C_active = findall(x -> length(x) > 0, C_samples)
                nCuts = nCuts + length(C_active)
                for cluster in C_active
                    objs[cluster], ∇objs[cluster,:] , feas_status = cMCMDCutting_plane(cmcmd_prob, config, z_cur, b_override = b[:,:,cluster_partition[cluster]], R_sample = C_samples[cluster])
                    if isnan(objs[cluster])
                        throw(DomainError(objs, "Infeasibility for KCut (not being handled yet)"))
                        con = @build_constraint(sum(z_cur[e] * (1 - z[e]) for e in 1:n) + sum((1- z_cur[e]) * z[e] for e in 1:n) >= 1)
                        MOI.submit(m1, MOI.LazyConstraint(cb), con)
                    else
                        con = @build_constraint(t[cluster] >= objs[cluster] + dot(∇objs[cluster,:], z - z_cur[1:n]))
                        MOI.submit(m1, MOI.LazyConstraint(cb), con)
                    end
                end
                # Store Lazy Multi-Cut in Cut Pool
                global countLazy, cutPoolLazy
                countLazy = countLazy + length(C_active)
                s0 = PrimalSolution( findall(z_cur .> 0), z_cur, sum(objs[cluster] for cluster in 1:n_clusters) + dot(c,z_cur), [objs[cluster] - dot(∇objs[cluster,:], z_cur[1:n]) for cluster in 1:n_clusters], ∇objs, false, method, R_sample, C_samples, Dict())
                push!(cutPoolLazy, s0)
            elseif method == "scp_adapt"
                nCuts = nCuts + 1
                # P_active: the keys of P (partition) that appear in R_sample
                # P_samples : the indices of the scenarios in each partition that appear in R_sample
                P_sample = sample_dict_elements(P, R_div)
                objs_p = Dict{String, Float64}
                ∇objs_p = Dict{String, Vector{Float64}}
                objs_p = Dict(s => 0.0 for s in keys(P_sample))
                ∇objs_p = Dict(s => zeros(n) for s in keys(P_sample))
                for s in keys(P_sample)
                    indices_s = P[s] # could be P_sample[s] but its the same
                    obj, ∇obj = cMCMDCutting_plane(cmcmd_prob, config, z_cur, b_override=b_P[s])
                    add_infeasibility_cut = isnan(obj)
                    if !add_infeasibility_cut
                        objs_p[s] = obj
                        ∇objs_p[s] = ∇obj[:]
                        con = @build_constraint(sum(t[indices_s]) / length(indices_s) >= obj + dot(∇obj, z .- z_cur[1:n]))
                        MOI.submit(m1, MOI.LazyConstraint(cb), con)
                        global countLazy, cutPoolLazy
                        countLazy = countLazy + length(P)
                        s0 = PrimalSolution( findall(z_cur .> 0), z_cur, sum(objs_p[s] for s in keys(P_sample)) + dot(c,z_cur), [objs_p[s] - dot(∇objs_p[s], z_cur[1:n]) for s in keys(P_sample)], Matrix(hcat([∇objs_p[s] for s in keys(P_sample)]...)'), false, method, R_sample, [], P_sample)
                        push!(cutPoolLazy, s0)
                    end
                end
            end

            if add_infeasibility_cut
                p_certificate, β_certificate, infeas_r = get_infeas_certificate(cmcmd_prob, z_cur, R_sample = R_sample)
                println("[Newcut] Adding Infeasibility Cut\n\t\tR_sample = $(R_sample)\n\t\tinfeas_r = $(R_sample[infeas_r])")
                if !(any(isnan, p_certificate) & any(isnan, β_certificate)) # check get_infeas_certificate terminated ok
                    for r in infeas_r
                        if (any(β_certificate[:,r] .> 0) || any(p_certificate[:,:,r] .!= 0))
                            r_tilde = R_sample[r]
                            LHS_check = sum(p_certificate[:,:,r] .* b[:,:,r_tilde]) - sum(z_cur[e] * u[e] * β_certificate[e, r] for e in 1:n)
                            println("\t\t\tAdding cut for r = $(r) (r_tilde = $(r_tilde)) [LHS_check at z_cur = $(LHS_check)]")
                            con = @build_constraint(sum(p_certificate[:,:,r] .* b[:,:,r_tilde]) - sum(z[e] * u[e] * β_certificate[e, r] for e in 1:n) <= 0)
                            MOI.submit(m1, MOI.LazyConstraint(cb), con)
                        end
                    end
                    infeas_cut = InfeasibleCut(z_cur, R_sample, p_certificate, β_certificate, infeas_r)
                    global cutPoolInfeas
                    push!(cutPoolInfeas, infeas_cut)
                end
            end
            
            if method == "scp_slim"
                global t_cur = callback_value(cb, t);                        
                global objVal_cur = dot(c,z_cur) + t_cur;
                global objVal_inner = dot(c,z_cur) + obj;
                println("[Newcut End Logging]\n  --objbst = $(objbst)\n  --objbnd = $(objbnd)\n  --mip_gap = $(mip_gap)\n  --objVal_cur = $(objVal_cur)\n  --objVal_inner = $(objVal_inner)")
            end
        end
    end

    MOI.set(m1, Gurobi.CallbackFunction(), Newcut)

    time_setup = time() - start_time
    time_remaining_setup = TIME_LIMIT - time_setup
    println("Model Setup Complete in t = $(time_setup)s")
    println("TIME LIMIT: $(TIME_LIMIT)")
    println("TIME SPENT IN SETUP: $(time_setup)")
    println("TIME REMAINING AFTER SETUP: $(time_remaining_setup)")
    if time_remaining_setup <= 0.05*TIME_LIMIT
        println("took too long to setup, setting to 1/4 the TIME_LIMIT")
        set_optimizer_attribute(m1, "TimeLimit", TIME_LIMIT/4)
    else
        set_optimizer_attribute(m1, "TimeLimit", time_remaining_setup)
    end
    


    # Solve the model and get the optimal solutions
    objtrue_iters = 0 
    objtrue_maxiters = 20
    global countLazy = 0
    global cutPoolLazy = PrimalSolution[]
    global cutPoolPartial = PrimalSolution[]
    global cutPoolInfeas = InfeasibleCut[]

    objtrue = -1
    objgaptrue = 1
    objgap = 1
    objval = -Inf
    zopt = zeros(n)
    objbnd = -1
    confidence_adjusted_bound_gap = +1
    # time_setup, time_opt

    for iter in 1:objtrue_maxiters
        global cutPoolLazy_added = false # Toggle whether the cutpool has been added during the Single-tree BnB
        optimize!(m1)
        objtrue = -1
        objgaptrue = 1
        objgap = 1
        objval = -Inf

        terminated_ok = ((termination_status(m1) == MOI.OPTIMAL) || 
                        (termination_status(m1) == MOI.SLOW_PROGRESS) || 
                        termination_status(m1) == MOI.TIME_LIMIT) && 
                        has_values(m1)
        println("terminated with values: $(terminated_ok) --- termination_status(m1) = $(termination_status(m1))")
        if termination_status(m1) == MOI.TIME_LIMIT
            println("TIME LIMIT REACHED")
        elseif termination_status(m1) == MOI.SLOW_PROGRESS
            println("SLOW PROGRESS")
        end
        time_spent = time() - start_time
        time_remaining = TIME_LIMIT - time_spent
        println("TIME LIMIT: $(TIME_LIMIT)")
        println("TIME SPENT: $(time_spent)")
        println("TIME REMAINING: $(time_remaining)")
        

        if terminated_ok
            println("IN TERMINATED_OK PART")
            zopt = value.(z)
            is_zopt_feasible = cMCMD_is_feasible(A, b, zopt, u);
            println("[Outer Loop] Zopt Feasibility: $(is_zopt_feasible)")
            objval = objective_value(m1)
            objbnd = objective_bound(m1)
            objgap = abs(objval - objbnd) / abs(objval) # Gap calculation according to Gurobi                

            if !is_zopt_feasible
                println("[Outer Loop] Adding zopt infeas cuts")
                @constraint(m1, sum(zopt[e] * (1 - z[e]) for e in 1:n) + sum((1- zopt[e]) * z[e] for e in 1:n) >= 1)
                @constraint(m1, sum((1- zopt[e]) * z[e] for e in 1:n) >= 1)
                is_z_last_feasible = cMCMD_is_feasible(A, b, z_last, u);
                if is_z_last_feasible
                    zopt = z_last;
                    objval = objval_last
                    objgap = gapcalculate(objval, objbnd)
                    objtrue, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt)
                    objgaptrue = gapcalculate(objtrue, objbnd)
                    println("[Outer Loop] Model timed out, but last z is feasible:")
                    println(" \t\tobjgap = $(objgap), objbnd = $(objbnd), objtrue = $(objtrue), objgaptrue = $(objgaptrue)")
                else
                    println("[Outer Loop] Adding z_last infeas cuts")
                    @constraint(m1, sum(z_last[e] * (1 - z[e]) for e in 1:n) + sum((1- z_last[e]) * z[e] for e in 1:n) >= 1)
                    @constraint(m1, sum((1- z_last[e]) * z[e] for e in 1:n) >= 1)
                    if time_remaining > 0
                        set_optimizer_attribute(m1, "TimeLimit", time_remaining)
                        continue
                    end
                end
            end
            objtrue, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt)
            objgaptrue = gapcalculate(objtrue, objbnd)
            println("[Outer Loop] Finished Iteration $(iter):\n\t\tObjval = $(objval) - Objbnd = $(objbnd) - Objgap = $(objgap) - Objtrue = $(objtrue) - Objgaptrue = $(objgaptrue)")
        else
            println("IN THE OTHER SIDE")
            is_scp_feasible = cMCMD_is_feasible(A, b, z_last, u);
            println("IS_SCP_FEASIBLE: $(is_scp_feasible) [For z_last]")
            println("Z_Last has dim: $(size(z_last))")
            if is_scp_feasible
                zopt = z_last;
                objval = objval_last
                objbnd = bnd_last
                objgap = gap_last
                objtrue, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt)
                objgaptrue = gapcalculate(objtrue, objbnd)
                println("[Outer Loop] Model timed out, but last z is feasible")
                println(" --- objgap = $(objgap), objbnd = $(objbnd), objtrue = $(objtrue), objgaptrue = $(objgaptrue)")
            else
                zopt = -1;
                objval = -1;
                objgap = -1;
                println("[Outer Loop] Infeasible")
                @constraint(m1, sum(z_last[e] * (1 - z[e]) for e in 1:n) + sum((1- z_last[e]) * z[e] for e in 1:n) >= 1)      
                @constraint(m1, sum((1- z_last[e]) * z[e] for e in 1:n) >= 1)
                if time_remaining > 0
                    set_optimizer_attribute(m1, "TimeLimit", time_remaining)
                    continue
                else
                    zopt = cmcmd_prob.ub;
                    objtrue, _, _, _ = get_opt_x_cmcmd(cmcmd_prob, zopt)
                    objval = objtrue
                    objbnd = objective_bound(m1)
                    objgap = gapcalculate(objval, objbnd)
                    objgaptrue = gapcalculate(objtrue, objbnd)
                    println("Timed out with outer loop infeasible - returning Ones")
                    println(" --- objgap = $(objgap), objbnd = $(objbnd), objtrue = $(objtrue), objgaptrue = $(objgaptrue)")
                    time_opt = time() - start_time - time_setup
                    return objval, zopt, m1, objgap, time_setup, time_opt, -1, objtrue_iters, objbnd
                end
            end
        end
        println("Keeping on")

        # ################# 
        # CONFIDENCE EXIT
        # ################# 
        confidence_level = 0.90
        Z = quantile(Normal(0,1), (1 + confidence_level) / 2) # ~1.64

        # Calculate the upper confidence bound at the given confidence level
        global R_sample
        W = length(R_sample)
        println("W length is $(W)")
        
        objs, _ , feas_status = cMCMDCutting_plane(cmcmd_prob, config, zopt, R_sample = R_sample, returnObjs = true)
        if any(isnan.(objs))
            println("[Keeping On] - Zopt Infeasible - Continuing")
            @constraint(m1, sum(zopt[e] * (1 - z[e]) for e in 1:n) + sum((1- zopt[e]) * z[e] for e in 1:n) >= 1)
            @constraint(m1, sum((1- zopt[e]) * z[e] for e in 1:n) >= 1)
            if time_remaining > 0
                set_optimizer_attribute(m1, "TimeLimit", time_remaining)
                continue
            end
        end

        U_tilde = dot(c,zopt) + nR/W * sum(objs[R_sample])
        s_U = sqrt(1/(W-1) * sum((objs[R_sample] .- mean(objs[R_sample])).^2))
        upper_confidence_bound = U_tilde + Z * s_U / sqrt(W)
        lower_confidence_bound = U_tilde - Z * s_U / sqrt(W) # Only for logging purposes - not needed anywhere


        # Calculate the confidence-adjusted bound gap
        confidence_adjusted_bound_gap = (upper_confidence_bound - objval)/upper_confidence_bound
        e_parameter = -0.0001

        println("[Outter Loop] Checking Confidence Bound for Termination... \n\tObjgap: $(objgap), ObjGapTrue: $(objgaptrue) ---\n\tObj: $(objval), ObjTrue: $(objtrue)")
        # Check if the confidence-adjusted bound gap is negative
        if ((confidence_adjusted_bound_gap < e_parameter) && (method != "scp_adapt")) || (R_div == 1)
            println("[Outter Loop] Terminate the optimization process after $(objtrue_iters+1) outer iterations.\n\tObjgap: $(objgap), ObjGapTrue: $(objgaptrue)")
            println("\t Termination Conditions: \n\t\tconfidence_adjusted_bound_gap = $(confidence_adjusted_bound_gap), method = $(method), R_div = $(R_div)")
            println("\tupper_confidence_bound: $(upper_confidence_bound), objval: $(objval), confidence_adjusted_bound_gap = $(confidence_adjusted_bound_gap) ---")
            println("\tU_tilde: $(U_tilde), s_U: $(s_U), mean(objs[R_sample]) = $(mean(objs[R_sample])), mean(objs) = $(mean(objs))")
            break
        end
        println("---> NOT BREAKING ---")
        println("     upper_confidence_bound: $(upper_confidence_bound), objval: $(objval), confidence_adjusted_bound_gap = $(confidence_adjusted_bound_gap) ---")
        println("     U_tilde: $(U_tilde), s_U: $(s_U), mean(objs[R_sample]) = $(mean(objs[R_sample])), mean(objs) = $(mean(objs)) ---")
        println("     Iter Count: $(iter_count) -- Objgap: $(objgap), ObjGapTrue: $(objgaptrue) ---")
        println("     obj true ITERS: $(objtrue_iters) / $(objtrue_maxiters)")
        println("<--- ")
        
        if time_remaining <= 0.05*TIME_LIMIT
            println("BREAKING FOR TIME: $(time_remaining)S REMAINING")
            break
        end

        # ADD ONE DETERMINISTIC CUT AT LAST SOLUTION
        if method == "scp_slim"
            obj0, ∇obj0 , feas_status = cMCMDCutting_plane(cmcmd_prob, config, zopt)
            @constraint(m1, t >= obj0 + dot(∇obj0, z - zopt[1:n]))
        elseif method == "scp_fat"
            obj0 = zeros(nR)
            ∇obj0 = zeros(nR, n)
            obj0, ∇obj0  = cMCMDCutting_plane(cmcmd_prob, config, zopt, returnObjs=true)
            @constraint(m1, [r in 1:nR], t[r] >= obj0[r] + dot(∇obj0[r,:], z - zopt[1:n]))
        elseif method == "scp_hybrid"
            obj0 = zeros(nR)
            ∇obj0 = zeros(nR, n)
            obj0, ∇obj0  = cMCMDCutting_plane(cmcmd_prob, config, zopt, returnObjs=true)
            @constraint(m1, [r in 1:nR], t[r] >= obj0[r] + dot(∇obj0[r,:], z - zopt[1:n]))
        elseif method == "scp_kcut"
            obj0 = zeros(n_clusters)
            ∇obj0 = zeros(n_clusters, n)
            # obj0 = zeros(nR)
            # ∇obj0 = zeros(nR, n)
            for cluster = 1:n_clusters
                # obj0[cluster], ∇obj0[cluster,:] , feas_status = cMCMDCutting_plane(A,b[:,:,cluster_partition[cluster]],d,u,zopt,k,γ, ub)
                obj0[cluster], ∇obj0[cluster,:] , feas_status = cMCMDCutting_plane(cmcmd_prob, config, zopt, b_override=b[:,:,cluster_partition[cluster]])
                @constraint(m1, t[cluster] >= obj0[cluster] + dot(∇obj0[cluster,:], z - zopt[1:n]))
            end
        elseif method == "scp_adapt"
            LB = objective_value(m1)
            UB_list = Dict(s => Inf for s in keys(P))
            z_cur = zeros(n)
            z_cur[1:n] = 1 .* (value.(z) .> 0) #Binarize
            s_to_split = "-1"; s_left = []; s_right = []
            for s in keys(sort(P, by=x->length(x), rev=true)) #Go through all the subregions in P, by decreasing size
                indices_s = P[s]
                # objs, ∇objs, duals = cMCMDCutting_plane_newJuMP(A,b[:,:,indices_s],d,u,z_cur,k,γ,ub, returnObjs=true, adjustOffset=false)
                # objs, ∇objs, _, duals = cMCMDCutting_plane(A,b[:,:,indices_s],d,u,z_cur,k,γ,ub, returnObjs=true, adjustOffset=false)
                objs, ∇objs, _, duals = cMCMDCutting_plane(cmcmd_prob, config, z_cur, b_override=b[:,:,indices_s], returnObjs=true, adjustOffset=false)
                UB_list[s] = sum(objs)
                for (idxr, r) in enumerate(indices_s)
                    @constraint(m1, t[r] >= objs[idxr] + dot(∇objs[idxr,:], z - z_cur[1:n]))
                end

                offset_of_the_average = dot( mean(duals[:,:,r] for r in 1:length(indices_s)), mean(b[:,:,r] for r in indices_s))
                average_of_offsets = mean( dot(duals[:,:,r], b[:,:,indices_s[r]]) for r in 1:length(indices_s))

                if s_to_split <= "0" && (offset_of_the_average != average_of_offsets)
                    #Reformating of the dual variables
                    duals_matrix = zeros(size(duals,1)*size(duals,2), size(duals,3))
                    [duals_matrix[:,r] = vec(duals[:,:,r]) for r in 1:length(indices_s)]
                    #Splitting into 2 clusters
                    R = kmeans(duals_matrix, 2; maxiter=100, display=:none)
                    cluster_assign = assignments(R) # get the assignments of points to clusters
                    s_to_split = s
                    s_left = indices_s[findall(cluster_assign .== 1)]
                    s_right = indices_s[findall(cluster_assign .== 2)]
                end
            end
            UB = dot(c, z_cur) + sum(UB_list[s] for s in keys(UB_list))
            if s_to_split >= "0"
                delete!(P, s_to_split)
                P[s_to_split*"_0"] = s_left
                P[s_to_split*"_1"] = s_right
                b_P = Dict(s => mean(b[:,:,P[s]],dims=3) for s in keys(P)) #Average b across all scenarions within each subset of P
            end
            if s_to_split < "0" || (UB-LB)/UB < 0.01 || time_remaining <= 0.05*TIME_LIMIT
                break
            end
        end

        # ADD LAZY CUTS AS CONSTRAINTS - Note: this adds last lazy cut as a normal constraint in order to reset lazy constraints / allow re-optimization, then within newcut we add all stored lazy cuts as lazy constraints again
        if length(cutPoolLazy) > 0
            # ADD THE LATEST LAZY CONSTRAINT AS NORMAL CONSTRAINT TO RESET THE LAZY CONSTRAINTS
            cLazy = cutPoolLazy[length(cutPoolLazy)]
            if cLazy.method == "scp_slim" && method == "scp_slim"
                @constraint(m1, t >= cLazy.offset[1] + dot(cLazy.slope, z))
            elseif cLazy.method == "scp_fat" && method == "scp_fat"
                # for r = 1:nR
                for r in cLazy.R_sample
                    @constraint(m1, t[r] >= cLazy.offset[r] + dot(cLazy.slope[r,:], z))
                end
            elseif cLazy.method == "scp_hybrid" && method == "scp_hybrid"
                # for r = 1:nR
                for r in cLazy.R_sample
                    @constraint(m1, t[r] >= cLazy.offset[r] + dot(cLazy.slope[r,:], z))
                end
                excluded_indices = [i for i in 1:nR if i ∉ cLazy.R_sample]
                total_sum = sum([cLazy.offset[r] + dot(cLazy.slope[r,:], z) for r in excluded_indices], init=0)
                @constraint(m1, sum(t[excluded_indices]) >= total_sum)
            elseif cLazy.method == "scp_kcut" && method == "scp_kcut"
                C_active = findall(x -> length(x) > 0, cLazy.C_samples) # Clusters that actually appear in the sample
                for cluster in C_active
                    @constraint(m1, t[cluster] >= cLazy.offset[cluster] + dot(cLazy.slope[cluster,:], z))
                end
            elseif cLazy.method == "scp_adapt" && method == "scp_adapt"
                P_sample = cLazy.P
                for (idx, s) in enumerate(keys(P_sample))
                    indices_s = P_sample[s]
                    @constraint(m1, sum(t[indices_s]) / length(indices_s) >= cLazy.offset[idx] + dot(cLazy.slope[idx, :], z))
                end
            end
        end

        R_div = max(Int(ceil(R_div/1.5)), 1)
        println("New R_div = $(R_div)")
        objtrue_iters += 1
        set_optimizer_attribute(m1, "TimeLimit", time_remaining)
    end

    time_opt = time() - start_time - time_setup

    terminated_ok = ((termination_status(m1) == MOI.OPTIMAL) || 
                        (termination_status(m1) == MOI.SLOW_PROGRESS) || 
                        termination_status(m1) == MOI.TIME_LIMIT) && 
                        has_values(m1)
    println("METHOD = $(method) - Terminated with values: $(terminated_ok)")       
    objgap = gapcalculate(objval, objbnd)
    return objval, zopt, m1, objgap, time_setup, time_opt, confidence_adjusted_bound_gap, objtrue_iters, objbnd
end


function get_opt_x_cmcmd(cmcmd_prob,zopt) # RAW NAIVE SOLVER FOR MC MULTIDAY
    A, b, c, d, u, k, γ, lb, ub, m, n, nC, nR = cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.c, 
                cmcmd_prob.d, cmcmd_prob.u, cmcmd_prob.k, cmcmd_prob.γ, cmcmd_prob.lb, cmcmd_prob.ub, 
                cmcmd_prob.n_nodes, cmcmd_prob.n_edges, cmcmd_prob.n_commodities, cmcmd_prob.n_days
    
    m0 = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(m0, "OutputFlag", 1) 
    set_optimizer_attribute(m0, "TimeLimit", 900)
    set_optimizer_attribute(m0, "FeasibilityTol", 1e-5)

    # x ∈ R^(|E| × K)
    @variable(m0, x[e in 1:n,1:nC,1:nR]>=0)
    # y ∈ R^(|E|)
    @variable(m0, y[e in 1:n, 1:nR])
    # y = ∑x
    @constraint(m0, [e in 1:n, r=1:nR], y[e,r] == sum(x[e,c,r] for c = 1:nC))
    # Ax = b
    @constraint(m0, [i=1:m, c=1:nC, r=1:nR], dot([A[i,e] for e in 1:n], [x[e,c,r] for e in 1:n]) == b[i,c,r])
    # y ≤ u
    @constraint(m0, [e in 1:n, r=1:nR], y[e,r] .<= u[e] * zopt[e])
    
    @objective(m0, Min, dot([c[e] for e in 1:n],zopt[1:n]) + sum(dot([d[e] for e in 1:n],y[:,r]) for r in 1:nR) +  1/(2*γ) * sum(sum(y[e, r]^2 for e in 1:n) for r = 1:nR))
    
    optimize!(m0)
    status = termination_status(m0)
    try 
        xopt = value.(x)
        objval = objective_value(m0)
        println("[get_opt_x_cmcmd] Evaluation Solved\n\t\tTermination Status: $(termination_status(m0))\n\t\tObjective Value: $(objval)")
        return objval, zopt, m0, xopt, 0
    catch
        zopt, xopt, objval = -1, -1, -1
        println("[get_opt_x_cmcmd] Infeasible")
        return objval, zopt, m0, xopt, 0
    end
end


function MCMD_is_feasible(A, b, lb)

    feasibility = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(feasibility, "OutputFlag", 0) 
    set_optimizer_attribute(feasibility, "TimeLimit", 600) 

    m , n = size(A)
    nC = size(b)[2]
    nR = size(b)[3]
    bigM = sum([sum(abs.(b[:,:,r]))/2 for r in 1:nR])

    @variable(feasibility, x[e in 1:n, 1:nC, 1:nR] >= 0)
    @constraint(feasibility, [i=1:m, c=1:nC, r=1:nR], dot([A[i,e] for e in 1:n], [x[e,c, r] for e in 1:n]) == b[i,c, r])
    @constraint(feasibility, [e in 1:n], x[e, :, :] .<= lb[e] * bigM)
    @objective(feasibility, Max, 0.)
    optimize!(feasibility)
    return (termination_status(feasibility) == MOI.OPTIMAL)
end


function cMCMD_is_feasible(A, b, lb, u)
    m , n = size(A)
    nC = size(b)[2]
    nR = length(size(b)) == 2 ? 1 : size(b)[3]
    
    feasibility = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(feasibility, "OutputFlag", 0) 
    set_optimizer_attribute(feasibility, "TimeLimit", 600) 
    
    # x ∈ R^(|E| × K)
    @variable(feasibility, x[1:n, 1:nC, 1:nR] >= 0)
    # y ∈ R^(|E|)
    @variable(feasibility, y[1:n,1:nR])
    # y = ∑x
    @constraint(feasibility, [e in 1:n, r=1:nR], y[e,r] == sum(x[e,c,r] for c = 1:nC))
    # Ax = b
    @constraint(feasibility, [i=1:m, c=1:nC, r=1:nR], dot([A[i,e] for e in 1:n], [x[e,c, r] for e in 1:n]) == b[i,c,r])
    # y ≤ u
    @constraint(feasibility, [e in 1:n, r=1:nR], y[e,r] .<= u[e] * lb[e])

    @objective(feasibility, Min, 0.)
    optimize!(feasibility)
        
    return (termination_status(feasibility) == MOI.OPTIMAL)
end



# NAIVE SOLVER
function cMCMDNaiveSolver(cmcmd_prob, config; relax = false)
    A, b, c, d, u, k, γ, lb, ub, m, n, nC, nR = cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.c, 
            cmcmd_prob.d, cmcmd_prob.u, cmcmd_prob.k, cmcmd_prob.γ, cmcmd_prob.lb, cmcmd_prob.ub, 
            cmcmd_prob.n_nodes, cmcmd_prob.n_edges, cmcmd_prob.n_commodities, cmcmd_prob.n_days

    m0 = Model(optimizer_with_attributes(
            Gurobi.Optimizer, "TimeLimit" => config.TIME_LIMIT))
    
    # x ∈ R^(|E| × K)
    @variable(m0, x[e in 1:n,1:nC,1:nR]>=0)
    # y ∈ R^(|E|)
    @variable(m0, y[e in 1:n, 1:nR])
    # z ∈ Z^(|E|)
    if relax
        @variable(m0, z[e in 1:n] >=0)
    else
        @variable(m0, z[e in 1:n], Bin)
    end
    JuMP.set_start_value.(z, lb[1:n]) 
    if config.use_z0
        println("using z0")
        @constraint(m0,[e in 1:n], z[e] .>= lb[e])
    end
    @constraint(m0,[e in 1:n], z[e] .<= ub[e])
    # ∑z ≤ k
    @constraint(m0, sum(z[e] for e in 1:n) <= k)
    # y = ∑x
    @constraint(m0, [e in 1:n, r=1:nR], y[e,r] == sum(x[e,c,r] for c = 1:nC))
    # Ax = b
    @constraint(m0, [i=1:m, c=1:nC, r=1:nR], dot(A[i,1:n], x[:,c,r]) == b[i,c,r])
    # y ≤ u
    @constraint(m0, [e in 1:n, r=1:nR], y[e, r] <= u[e] * z[e])    
    
    @objective(m0, Min, dot(c[1:n],z) + sum(dot(d[1:n],y[:,r]) for r in 1:nR) +  1/(2*γ) * sum(sum(y[e, r]^2 for e in 1:n) for r = 1:nR))
    
    optimize!(m0)
    
    zopt = value.(z)
    xopt = value.(x)
    objval = objective_value(m0)
    
    if relax
        objgap = 0;
    else
        objbnd = objective_bound(m0)
        objgap = abs(objval - objbnd) / abs(objval) # Gap calculation according to Gurobi
    end
   
    println("Model Solved")
    return objval, zopt, m0, xopt, objgap, objbnd
end


function cMCMDGetNetwork(config::Config; ntries_inner = 5, sparse_ratio = nothing, k_nn = nothing)
    m, nC, k, nR, ntries, n = config.m, config.nC, config.k, config.nR, config.ntries, config.n
    Random.seed!(config.randomSeed)
    u = zeros(n)
    for i = 1:ntries
        A, b, c, d, lb, ub, edge_map, outgoing_edges, incoming_edges, old_to_new_map = MCMDGetNetwork(m,nC,k, nR,sparse_ratio = nothing, k_nn = k_nn)
        n = size(A)[2]
        b_tmp = maximum(b, dims=3)
        mult_factor = 4
        B = sum([ sum(b_tmp[b_tmp[:, j] .> 0, j]) for j = 1:nC])*mult_factor*(1+m/10)
        n_edges = sum(lb)
        for j = 1:ntries_inner
            u = vec(round.(Int, (2*rand(n) .+ 3.)*B/n_edges/5))
            is_cmcmd_feasible = cMCMD_is_feasible(A, b, lb, u)
            if is_cmcmd_feasible
                println("Found feasibility after $(i) outer and $(j) inner tries")
                return A, b, c, d, lb, ub, u, edge_map, outgoing_edges, incoming_edges, old_to_new_map   # SHOULD RETURN HERE
            elseif mod(i,Int(floor(ntries/5))) == 0
                println("Did Not Find Feasibility after $(i) tries")         
            elseif i == ntries
                println("Did Not Find Feasibility after $(i) tries - ever")         
            end
        end
    end
    return A, b, c, d, lb, ub, u, edge_map, outgoing_edges, incoming_edges, old_to_new_map
end

function cMCMDGetProblem(config::Config)
    if !config.isBenchmark
        A, b, c, d, lb, ub, u, edge_map, outgoing_edges, incoming_edges, old_to_new_map = cMCMDGetNetwork(config, k_nn = config.k_nn)
        config.n = length(edge_map)
        config.k = min(config.k, config.n)
        config.nC = size(b)[2]
        cmcmd_prob = cMCMDProblem(config.m, config.n, config.nC, config.nR, config.n_clusters, config.k, 
                                  A, b, c, d ./ config.nR, u, lb, ub, config.γ, config.R_div, config.R_div_kelley, 
                                  edge_map, outgoing_edges, incoming_edges, old_to_new_map)
    else
        cmcmd_prob = getBenchmarkProblem(config)
    end
    return cmcmd_prob
end

function getBenchmarkProblem(config::Config)
    try_n = 0
    cmcmd_prob=nothing
    while try_n < config.ntries
        try_n += 1
        A, b, c, d, lb, ub, u, m, n, nC, edge_map, outgoing_edges, incoming_edges, old_to_new_map = get_benchmark(config)
        config.n = length(edge_map)
        config.k = min(config.k, config.n)
        config.nC = nC
        cmcmd_prob = cMCMDProblem(m, n, nC, config.nR, config.n_clusters, config.k, 
                                  A, b, c, d ./ config.nR, u, lb, ub, config.γ * config.nR, config.R_div, config.R_div_kelley, 
                                  edge_map, outgoing_edges, incoming_edges, old_to_new_map)
        
        if cMCMD_is_feasible(cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.lb, cmcmd_prob.u)
            println("(Outer - cMCMDGetProblem) Feasibility found after $(try_n) outer tries")
            break
        elseif !config.use_z0
            println("(Outer - cMCMDGetProblem) Not using z0, no feasibility check (try_n = $(try_n))")
            break
        end
    end
    return cmcmd_prob
end



function MCMDGetNetwork(m = 20, nC = 2, k = 200, nR = 2;sparse_ratio = nothing, k_nn = nothing)
    # Define no of connections
    n = m * m
    
    # DEFINE b - nR
    b = round.(Int, 5 .+ 20 .* rand(m,nC,nR));
    perm = shuffle(1:m)
    # One node at random provides for all the others (column wise) -> sum(b[:,i]) = 0
    for j in 1:nR
        for i in 1:nC
            b[perm[i], i, j] = - (sum(b[:,i, j]) - b[perm[i], i, j])
        end
    end

    # DEFINE A
    A = zeros(m, n)     #Flow conservation matrix
    for i in 1:m
        for j in 1:m
            if j != i
                k = (i-1)*m + j
                A[i,k] = 1
                k = (j-1)*m + i
                A[i,k] = -1
            end
        end
    end

    # DEFINE lb - Initial connections
    lb = zeros(n)
    list = [(i,j) for i in 1:m for j in (i+1):m]
    cover = []
    while true  #Find a spanning tree
        subset = shuffle(list)[1:2*m]
        cover = Set{Int}()
        for e in subset
            push!(cover, e[1])
            push!(cover, e[2])
        end
        if length(cover) >= m
            for e in subset
                lb[(m*(e[2]-1) + e[1])] = 1.
                lb[(m*(e[1]-1) + e[2])] = 1.
            end
            break
        end
    end

    repetitions = 1
    if m >= 100
        repetitions = Int(floor(m/10))
    end

    println("doing $(repetitions) repetitions")
    while !MCMD_is_feasible(A, b, lb)    #Add edges until first_stage_feasible
        println("enter is feasible")
        for repetition in 1:repetitions
            i = rand(findall(lb .< .5))
            lb[i] = 1
        end
    end

    @assert MCMD_is_feasible(A, b, lb)
    
    ## Costs
    Nodes_locations = rand(m,2)
    c = 1 .+ 3 .* rand(n) #construction cost U(1,4)
    c[findall(lb .> .5)] .= 0. #No cost for existing edges
    c *= 20

    d = zeros(n) #transportation cost
    for i in 1:m
        for j in (i+1):m
            z = Nodes_locations[i,:] - Nodes_locations[j,:]
            d[(m*(j-1) + i)] = LinearAlgebra.norm(z)
            d[(m*(i-1) + j)] = d[(m*(j-1) + i)]
        end
    end
    d *= 10

    # DEFINE ub - All connections
    ub = ones(n)
    ub[[(m*(i-1) + i) for i in 1:m]] .= 0. #Do not allow for loops
    lb[[(m*(i-1) + i) for i in 1:m]] .= 0. #Do not allow for loops

    # Randomly disable sparse_ratio % of nodes other than initial feasible ones
    isNaiveSparse = !isnothing(sparse_ratio)
    isKnnSparse = !isnothing(k_nn)
    if isNaiveSparse
        println("enter naive")
        ub_keep = sample((1:n)[ub .> lb], Int(ceil(length((1:n)[ub .> lb])*sparse_ratio)), replace=false)
        ub_kill = setdiff((1:n)[ub .> lb], ub_keep)
        ub[ub_kill] .= 0
    elseif isKnnSparse  # Disable All but KNN of each node, except for initial feasible lb
        println("enter knn")
        for i in 1:m
            node_indices = ((i-1)*m+1:i*m)
            node_ub = (node_indices)[ub[node_indices] .> lb[node_indices]]
            node_dists = d[node_ub]
            node_knn = partialsortperm(node_dists, 1:min(k_nn, length(node_dists)))
            node_knn_edges = node_ub[node_knn]
            # d[node_knn_edges]
            ub_kill = setdiff(node_indices[ub[node_indices] .> lb[node_indices]], node_knn_edges)
            ub[ub_kill] .= 0
        end
    end


    # SELECT EDGES AND REDEFINE EVERYTHING
    available_edges, edge_map, outgoing_edges, incoming_edges, old_to_new_map = create_edge_mappings(m, ub)

    A_filtered = A[:, available_edges]
    b_filtered = b
    lb_filtered = lb[available_edges]
    ub_filtered = ub[available_edges]
    c_filtered = c[available_edges]
    d_filtered = d[available_edges]

    return A_filtered, b_filtered, c_filtered, d_filtered, lb_filtered, ub_filtered, edge_map, outgoing_edges, incoming_edges, old_to_new_map
end


function get_benchmark(config::Config)
    Rx, Ry, nR, demand_lb, demand_ub, ntries, use_file_demands, use_z0, corr_pcnt = config.Rx, config.Ry, config.nR, 
                    config.demand_lb, config.demand_ub, config.ntries, config.use_file_demands, config.use_z0, config.corr_pcnt
    
    Rx = lpad(string(Rx), 2, '0')
    words = readlines("./benchmarks/R/r$(Rx).$(Ry).dow", keep=true)
  
    # GET NODES, ARCS, COMMODITIES 
    m, n_edges, nC = r_line2nums(words[2])
    println("Benchmark Problem with : (m, n_edges, nC) = $m, $n_edges, $nC")
    n = m * m 
    A = zeros(m, n) 
    ub = zeros(m*m)
    c = zeros(n)
    d = zeros(n)
    u = zeros(n)
    i_list = []
    j_list = []
    
    
    # DEFINE A
    A = zeros(m, n)     #Flow conservation matrix
    for i in 1:m
        for j in 1:m
            if j != i
                k = (i-1)*m + j
                A[i,k] = 1
                k = (j-1)*m + i
                A[i,k] = -1
            end
        end
    end
    # CONSTRUCT c, d, u, ub
    for n_edge = 1:n_edges
        i, j, d_ij, u_ij, c_ij, = r_line2nums(words[2+n_edge])
        push!(i_list, i)
        push!(j_list, j)
        n_index = ij2n_index(i,j,m)
        # edges, costs, capacity become original direction
        ub[n_index] = 1
        c[n_index] = c_ij
        d[n_index] = d_ij
        u[n_index] = u_ij
    end
    
    # CONSTRUCT LB
    lb = zeros(n)
    list = [(i_list[arc],j_list[arc]) for arc in 1:n_edges]
    cover = []
    if config.use_z0
        println("Using z0")
        while true  #Find a spanning tree
            subset = shuffle(list)[1:2*m]
            cover = Set{Int}()
            for e in subset
                push!(cover, e[1])
                push!(cover, e[2])
            end
            if length(cover) >= m
                for e in subset
                    lb[(m*(e[2]-1) + e[1])] = 1.
                    lb[(m*(e[1]-1) + e[2])] = 1.
                end
                break
            end
        end
    else
        println(" --- Not using z0 ---")
    end    
    # CONSTRUCT b (DEMANDS)
    if use_file_demands
        println("Using file demands")
        scenarios_file_path = "./benchmarks/instances-R/r$(Rx)-$(corr_pcnt)-1000"
        println("READING: $(scenarios_file_path)")
        all_scenarios = read_scenarios(scenarios_file_path)
        selected_scenarios = all_scenarios[1:nR]

         # CONSTRUCT b scenarios (DEMANDS)
         b = zeros(m, nC, nR)
         for r in 1:nR
             demands = selected_scenarios[r]
             for n_commodity in 1:nC
                origin_i, dest_j, = r_line2nums(words[2+n_edges+n_commodity])
                b[origin_i, n_commodity, r] = demands[n_commodity]
                b[dest_j, n_commodity, r] = -demands[n_commodity]
             end
         end
    else
        println("Using uniform demands")
        # CONSTRUCT b NOMINAL (DEMANDS)
        b_nominal = zeros(m, nC)
        for n_commodity in 1:nC
            origin_i, dest_j, b_ij = r_line2nums(words[2+n_edges+n_commodity])
            b_nominal[origin_i, n_commodity] = - b_ij
            b_nominal[dest_j, n_commodity] = b_ij
        end
    
        # CONSTRUCT b scenarios (DEMANDS)
        if (config.demand_lb isa Nothing)
            config.demand_lb = 0.99
        end
        if (demand_ub isa Nothing)
            demand_ub = 1.01        
        end
    
        α = zeros(nC, nR)
        b = zeros(m, nC, nR)
    
        for r in 1:nR
            n_try = 0
            correlated_scenarios = sample(1:nR, Int(ceil(nR * corr_pcnt)), replace = false)
    
            while n_try <= ntries
                n_try = n_try + 1
                for k in 1:nC
                    α[k, r] = rand(Uniform(config.demand_lb, demand_ub))
                    b[:, k, r] = round.(α[k, r] .* b_nominal[:,k])
                end

                β = rand(Uniform(config.demand_lb, demand_ub), nC)
                for r in correlated_scenarios
                    for k in 1:nC
                        α[k, r] = β[k] + rand(Uniform(-(1-config.demand_lb)/4, (demand_ub-1)/4))
                        b[:, k, r] = round.(α[k, r] .* b_nominal[:,k])
                    end
                end
    
                if cMCMD_is_feasible(A, b[:, :, r], ub, u)
                    break
                end
            end    
            println(cMCMD_is_feasible(A, b[:, :, correlated_scenarios], ub, u))
            if cMCMD_is_feasible(A, b[:, :, correlated_scenarios], ub, u)
                println("found feasible after $(n_try) tries")
            end
        end
    end

    # CONSTRUCT LB (Initial spanning tree)
    repetitions = 1
    println("doing $(repetitions) repetitions")
    while !cMCMD_is_feasible(A, b, lb, u) && config.use_z0   #Add edges until first_stage_feasible
        println("enter is feasible for benchmark")
        println("doing $(repetitions) repetitions:")
        for repetition in 1:repetitions
            println(" try $(repetition)")
            i = rand(findall((lb .< .5) .& (ub .> .5)))
            lb[i] = 1
        end
    end

    #SELECT EDGES AND REDEFINE EVERYTHING
    available_edges, edge_map, outgoing_edges, incoming_edges, old_to_new_map = create_edge_mappings(m, ub)
    A_filtered = A[:, available_edges]
    b_filtered = b
    lb_filtered = lb[available_edges]
    ub_filtered = ub[available_edges]
    c_filtered = c[available_edges]
    d_filtered = d[available_edges]
    n = length(available_edges)
    u_filtered = u[available_edges]
  
    return A_filtered, b_filtered, c_filtered, d_filtered, lb_filtered, ub_filtered, u_filtered, m, n, nC, edge_map, outgoing_edges, incoming_edges, old_to_new_map
end



function cluster_days(cmcmd_prob; cluster_target = nothing, config=nothing)
    n_clusters = cmcmd_prob.n_clusters
    n_nodes = cmcmd_prob.n_nodes
    n_commodities = cmcmd_prob.n_commodities
    n_days = cmcmd_prob.n_days
    b = cmcmd_prob.b
    if cluster_target isa Nothing
        println("Cluster target not specified, using demands instead")
        demand_vectors = zeros(n_nodes * n_commodities, n_days)
        # make a dataset with nR points with (m * nR) dimensions
        r = 1
        [demand_vectors[:, r] = vec(b[:,:,r]) for r in 1:n_days]
        cluster_target = demand_vectors
    elseif size(cluster_target)[2] != n_days
        println("Cluster target dim 2 not matching number of days\n (size(cluster_target)[2] = $(size(cluster_target)[2]), n_days = $(n_days))\n using demands instead")
        demand_vectors = zeros(n_nodes * n_commodities, n_days)
        [demand_vectors[:, r] = vec(b[:,:,r]) for r in 1:n_days]
        cluster_target = demand_vectors
    end
    
    # cluster into 4 clusters using K-means
    R = kmeans(cluster_target, n_clusters; maxiter=200, display=:final)
    
    cluster_assign = assignments(R) # get the assignments of points to clusters
    cluster_c = counts(R) # get the cluster sizes
    cluster_M = R.centers # get the cluster centers
    cluster_partition = [findall(cluster_assign .== i) for i in 1:n_clusters]
    return ClusteringInfo(cluster_partition, R)
end



function cMCMDCutting_plane(cmcmd_prob, config, z0; b_override = nothing, R_sample = nothing, adjustOffset = false, returnObjs = false, 
                            is_magnanti_wong_cut = false, z0_magnanti = nothing, use_si_vi = false)
    # THROW WARNING IF NOT 0 .<= Z0 .<= 1 
    !all(0 .<= z0 .<= 1) && @warn "z0 not in [0,1]"
    # CLIP z0 BETWEEN 0 AND 1 
    z0 = max.(0.0, min.(1.0, z0))
    A,b,d,u,k,γ,ub = cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.d, cmcmd_prob.u, cmcmd_prob.k, cmcmd_prob.γ, cmcmd_prob.ub
    if !(b_override isa Nothing) # override cmcmd_prob.b with b_override
        b = b_override
    end
    is_stochastic = !(R_sample isa Nothing)
    cur_ver = is_stochastic ? "Stochastic" : "Deterministic"
    m = size(A)[1]
    n = size(A)[2]
    nC = size(b)[2]
    nR = size(b)[3]
    println()
    println("[enter cMCMDCP]")
    nR_tilde = nR
    if is_stochastic
        nR_tilde = length(R_sample)
    end
    free_ones = findall(z0 .> 1e-10) #Variables >0  => appear in the objective

    if config.useMosek
        println("\n\tUsing Mosek")
        m2 = Model(optimizer_with_attributes(Mosek.Optimizer))
        if config.verbose_logging == 0
            set_optimizer_attribute(m2, "QUIET", true)
        end
        set_optimizer_attribute(m2, "MIO_MAX_TIME", 7200.0) 
    else
        println("\n\tUsing GRB")
        m2 = Model(optimizer_with_attributes(Gurobi.Optimizer))
        set_optimizer_attribute(m2, "OutputFlag", config.verbose_logging) 
        set_optimizer_attribute(m2, "TimeLimit", 7200) 
    end
    
    # Add variables
    @variable(m2, p[1:m, 1:nC, 1:nR_tilde])  # Nodes × K × R
    @variable(m2, α[e in 1:n, 1:nR_tilde])        # Edges × R
    @variable(m2, β[e in 1:n, 1:nR_tilde] >= 0)   # Edges × R
    
    # VI 3: SI
    SI_violated = false; b_si = zeros(n, nC, nR_tilde);
    println("USE SI VIOLATION: $(use_si_vi)")
    if use_si_vi
        # CHECK IF SI IS VIOLATED AT ALL
        destinations = [findall(x -> x > 0, b[:,i,1]) for i in 1:nC]
        b_si = [min.(u[e], sum(b[i, k, r] for i in destinations[k])) for e in 1:n, k in 1:nC, r in (is_stochastic ? R_sample : 1:nR_tilde)]
        _, _, _, xopt, _ = get_opt_x_cmcmd(cmcmd_prob, z0)
        if xopt == -1 # get_opt_x_cmcmd infeasible -> cMCMDCutting_plane infeasible
            return NaN, NaN, false
        end
        SI_violated = !all([(xopt[e,k,r]<=b_si[e,k,r] * z0[e]) for e in 1:n, k in 1:nC, r in 1:nR_tilde])
        if SI_violated
            println("ADDING S.I. V.I.")
            @variable(m2, lambda[e in 1:n, 1:nC, 1:nR_tilde] >= 0)   # Edges × R    
        else
            println("NO S.I. VIOLATION")
        end
    end

    # Add constraints
    for r in 1:nR_tilde
        for c=1:nC
            for e in free_ones
                if use_si_vi & SI_violated
                    @constraint(m2, 0 >=  dot(A[:,e], p[:,c,r]) - d[e] + α[e, r] - lambda[e, c, r])
                else
                    @constraint(m2, 0 >=  dot(A[:,e], p[:,c,r]) - d[e] + α[e, r])
                end
            end
        end
    end
    
    function calc_objective_expr_line(z0, u, β, p, b, α, nR_tilde)
        n = length(u)
        free_ones = findall(z0 .> 1e-10)
        return -sum(dot([u[e] for e in 1:n], β[:, r] .* z0[1:n]) for r in 1:nR_tilde, init=0) + 
               sum(p .* b, init=0) - 
               (γ / 2) * sum([z0[e] * α[e, r]^2 for r in 1:nR_tilde for e in free_ones], init=0)
    end
        
    # OBJECTIVE 
    @variable(m2, t[e in 1:n, 1:nR_tilde])
    @constraint(m2, [e in free_ones, r in 1:nR_tilde], t[e,r] >= α[e, r] + β[e, r])
    @constraint(m2, [e in free_ones, r in 1:nR_tilde], t[e,r] >= -(α[e, r] + β[e, r]))
    objective_expr = calc_objective_expr_line(z0, u, β, p, (is_stochastic ? b[:,:,R_sample] : b), t, nR_tilde)
    @objective(m2, Max, objective_expr)

    if use_si_vi & SI_violated
        println("[cMCMDCutting_Plane] Adjusting Objective Value for SI VI")
        previous_objective = objective_function(m2)
        new_objective = previous_objective - sum(lambda[e, k, r] * b_si[e, k, r] * z0[e] for e in 1:n, k in 1:nC, r in 1:nR_tilde)
        set_objective_function(m2, new_objective)
    end
    # OPTIMIZE
    optimize!(m2)
    terminated_ok = (termination_status(m2) == MOI.OPTIMAL)
    terminated_slowprog =  (termination_status(m2) == MOI.SLOW_PROGRESS)
    terminated_timelimit = (termination_status(m2) == MOI.TIME_LIMIT && has_values(m2))
    if terminated_slowprog
        terminated_ok = cMCMD_is_feasible(cmcmd_prob.A, cmcmd_prob.b[:,:,(is_stochastic ? R_sample : 1:nR)], z0, cmcmd_prob.u)
        println("[cMCMDCutting_plane SLOW_PROGRESS] feasibility = $(terminated_ok)")
    end

    # CHECK FOR SI VI 
    if !terminated_ok && (use_si_vi & SI_violated)
        println(" SI VI Termination Status: $(termination_status(m2))")
        set_objective_function(m2, previous_objective)
        @constraint(m2, [e in 1:n, c in 1:nC, r in 1:nR_tilde], lambda[e, c, r] == 0)
        optimize!(m2)
        SI_violated = false # Assume SI was never introduced
        terminated_ok = ((termination_status(m2) == MOI.OPTIMAL) || 
                        (termination_status(m2) == MOI.SLOW_PROGRESS) || 
                        termination_status(m2) == MOI.TIME_LIMIT) && 
                        has_values(m2)
        println(" SI VI --> NEW TERMINATION STATUS = $(termination_status(m2))")
    end

    # CHECK FOR MAGNANTI & WONG CUT
    added_mw_cut = false
    if terminated_ok && (is_magnanti_wong_cut ? !verify_core_point_simple(cmcmd_prob, z0) : false) # IF z0 is a core point - NO NEED TO GENERATE CUT
        added_mw_cut = true
        println(" [ADDING M & W CUT] $(cur_ver) Version")
        obj_bar = objective_value(m2)
        # SET NEW M&W OBJECTIVE
        objective_expr_magnanti = calc_objective_expr(z0_magnanti, u, β, p, (is_stochastic ? b[:,:,R_sample] : b), α, nR_tilde, 1:n)
        set_objective_function(m2, objective_expr_magnanti)
        
        # SET NEW M&W CONSTRAINT
        objective_expr_bar = calc_objective_expr(z0, u, β, p, (is_stochastic ? b[:,:,R_sample] : b), α, nR_tilde, 1:n)
        @constraint(m2, objective_expr_bar >= obj_bar)

        optimize!(m2)
        terminated_ok = ((termination_status(m2) == MOI.OPTIMAL) || 
                        (termination_status(m2) == MOI.SLOW_PROGRESS) || 
                        termination_status(m2) == MOI.TIME_LIMIT) && 
                        has_values(m2)        
        println(" [M & W CUT STATUS] $(termination_status(m2))")
    end


    feas_status = true
    if terminated_ok
        println("[cMCMDCutting_plane Termination] optimal")
    elseif terminated_timelimit
        println("[cMCMDCutting_plane Termination] suboptimal but with feasible solution")
    else
        println("[cMCMDCutting_plane Termination] infeasible")
        println("  \tTermination Status: $(termination_status(m2))\n\tR_sample = $(R_sample)")
        feas_status = false
    end
  
    if !(terminated_ok || terminated_timelimit)
        return NaN, NaN, feas_status
    else
        if is_stochastic
            println("----------------------------------------")
            println(" [cMCMDCutting_plane Exit - Stochastic] Terminated Ok")
            println("\tTermination status: $(termination_status(m2))\n\tOBJ VALUE -> $(objective_value(m2))\n\tR_sample = $(R_sample)")
            p_val_tilde = value.(p)
            β_val_tilde = value.(β)
            α_val_tilde = value.(α)
            λ_val_tilde = (use_si_vi & SI_violated) ? value.(lambda) : zeros(n, nC, nR_tilde)

            α_val_tilde[setdiff(1:n, free_ones),:] = -[max(maximum([dot(A[:,e], p_val_tilde[:,c,r]) + λ_val_tilde[e, c, r] for c in 1:nC]) - d[e] , 0) for e in setdiff(1:n, free_ones), r in 1:nR_tilde]
            #= Just for clarity, this is equivalent to:
                for r in 1:nR_tilde
                    for e in setdiff(1:n, free_ones)
                        # maximum ( max_k A^T p^[k] - d[e], 0)
                        max_k = max(maximum([dot(A[:,e], p_val_tilde[:,c,r]) for c in 1:nC]) - d[e] , 0)
                        α_val_tilde[e, r] = max_k
                    end
                end    (END COMMENT)=#
            
            α_val,β_val,p_val,λ_val = average_cut(cmcmd_prob, config, R_sample, α_val_tilde, β_val_tilde, p_val_tilde; λ_val_tilde = λ_val_tilde);
            
            if (use_si_vi & SI_violated)
                b_si_val = zeros(n, nC, nR)
                b_si_val[:,:,R_sample] = b_si
            end
            objs  = zeros(nR)
            ∇objs = zeros(nR, n)            
            for r = 1:nR
                objs[r] = - dot([u[e] for e in 1:n],β_val[:,r] .* z0[1:n]) + (sum(p_val[:,:,r] .* b[:,:,r])) - γ/2 * sum(z0[1:n] .* ((α_val[:, r].+β_val[:, r]).^2) )
                ∇objs[r,:] = - γ/2 * ((α_val[:, r].+β_val[:, r]).^2) .- (u[1:n] .* β_val[:,r])
                if (use_si_vi & SI_violated & (r in R_sample))
                    objs[r] -= sum(λ_val[:,:,r] .* b_si_val[:,:,r] .* z0)
                end
            end

            # calc obj 
            obj = sum(objs)
            ∇obj = sum(∇objs[r,:] for r in 1:nR)
            println("OBJ = $obj, objective_value = $(objective_value(m2))")

            # calc offset if needed
            if adjustOffset
                for r in 1:nR 
                    objs[r] -= dot(∇objs[r,:], z0[1:n])
                end
                obj = obj - dot(∇obj, z0[1:n])
            end

            println("[cMCMDCutting_plane Exit - Stochastic] OBJ = $obj, objective_value = $(objective_value(m2))\n********************************")
            if !returnObjs
                return obj, ∇obj, feas_status, p_val
            else
                return objs, ∇objs, feas_status, p_val
            end
        else
            println("----------------------------------------")
            println("[cMCMDCutting_plane Exit - Determinsitic]")
            p_val = value.(p)
            β_val = value.(β)
            α_val = value.(α)
            λ_val = (use_si_vi & SI_violated) ? value.(lambda) : zeros(n, nC, nR_tilde)
  
            α_val[setdiff(1:n, free_ones),:] = -[max(maximum([dot(A[:,e], p_val[:,c,r]) + λ_val[e, c, r] for c in 1:nC]) - d[e] , 0) for e in setdiff(1:n, free_ones), r in 1:nR]            
  
            obj = objective_value(m2)
            ∇obj = - γ/2 * sum((α_val .+ β_val).^2, dims=2)

            # calc objs vector
            objs = zeros(nR)
            ∇objs = zeros(nR, n)            
            for r = 1:nR
                objs[r] = - dot([u[e] for e in 1:n],β_val[:,r] .* z0[1:n]) + (sum(p_val[:,:,r] .* b[:,:,r])) - γ/2 * sum(z0[1:n] .* ((α_val[:, r].+β_val[:, r]).^2) )
                ∇objs[r,:] = - γ/2 * ((α_val[:, r].+β_val[:, r]).^2) .- (u[1:n] .* β_val[:,r])
                objs[r] += (use_si_vi & SI_violated) ? sum(λ_val[:,:,r] .* b_si_val[:,:,r] .* z0) : 0
            end

            # calc offset if needed
            if adjustOffset
                obj = obj - dot(∇obj, z0[1:n])
                for r in 1:nR 
                    objs[r] -= dot(∇objs[r,:], z0[1:n])
                end
            end

            if !returnObjs
                return obj, ∇obj, feas_status, p_val
            else
                return objs, ∇objs, feas_status, p_val
            end
        end
    end
end


function get_infeas_certificate(cmcmd_prob, z0; R_sample = nothing)
    A, b, d, k, ub, u, γ = cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.d, cmcmd_prob.k, cmcmd_prob.ub, cmcmd_prob.u, cmcmd_prob.γ
    m, n, nC, nR = size(A)[1], size(A)[2], size(b)[2], size(b)[3]
    nR_tilde = nR
    if !(R_sample isa Nothing)
        nR_tilde = length(R_sample)
        b = cmcmd_prob.b[:,:,R_sample];
    end
    free_ones = findall(z0 .> 1e-10) #Variables >0  => appear in the objective

    m_cert = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(m_cert, "QUIET", true)
    set_optimizer_attribute(m_cert, "MIO_MAX_TIME", TIME_LIMIT/10)

    @variable(m_cert, p[1:m, 1:nC, 1:nR_tilde])
    @variable(m_cert, β[e in 1:n, 1:nR_tilde] >= 0)
    ϵ = 1

    # Add constraints
    #    ∑_k <p[:,k,r], b[:,k,r]>   -   ∑_e (z0[e] * u[e] * β[e, r]) > 0
    @constraint(m_cert, [r in 1:nR_tilde], sum(sum(p[:,k,r] .* b[:,k,r]) for k in 1:nC) - sum(z0[e] * u[e] * β[e, r]  for e in 1:n) <= ϵ)
    @constraint(m_cert, [r in 1:nR_tilde, c in 1:nC, e in 1:n], 0 >= dot(A[:,e], p[:,c,r]) - β[e, r])

    @objective(m_cert, Max, sum(sum(sum(p[:,k,r] .* b[:,k,r]) for k in 1:nC) - sum(z0[e] * u[e] * β[e, r]  for e in 1:n) for r in 1:nR_tilde))

    
    optimize!(m_cert)
    terminated_ok = (termination_status(m_cert) == MOI.OPTIMAL) || (termination_status(m_cert) == MOI.SLOW_PROGRESS)

    obj_partial = []
    for r in 1:nR_tilde
        obj_r = value.(sum(sum(p[:,k,r] .* b[:,k,r]) for k in 1:nC) - sum(z0[e] * u[e] * β[e, r]  for e in 1:n))
        push!(obj_partial, obj_r)
    end
    infeas_r = findall(obj_partial .> ϵ/1000)
    println("[infeas_certificate] \n\tObj_partial = $obj_partial\n\tInfeas_r = $infeas_r\n\tTermination status = $(termination_status(m_cert))")
    terminated_ok = ((termination_status(m_cert) == MOI.OPTIMAL) || 
                        (termination_status(m_cert) == MOI.SLOW_PROGRESS) || 
                        termination_status(m_cert) == MOI.TIME_LIMIT) && 
                        has_values(m_cert)
    
    p_val = terminated_ok ? value.(p) : NaN
    β_val = terminated_ok ? value.(β) : NaN
    infeas_r = terminated_ok ? infeas_r : NaN
    !terminated_ok && @warn "[infeas_certificate]\n\t\tcertification problem is infeasible"
    return p_val, β_val, infeas_r
end

function average_cut(cmcmd_prob, config, R_sample, α_val_tilde, β_val_tilde, p_val_tilde; λ_val_tilde = nothing)
    # SET DUAL_TILDE VALUES
    m, n, nC, nR = config.m, config.n, config.nC, config.nR
    p_val = zeros(m, nC, nR)
    p_val[:,:, R_sample] = p_val_tilde

    β_val = zeros(n, nR)
    β_val[:, R_sample] = β_val_tilde

    α_val = zeros(n, nR)
    α_val[:, R_sample] = α_val_tilde      

    ############
    # AVG method_average_cut
    ############
    p_val[:,:,Not(R_sample)] .= mean(p_val[:,:, R_sample], dims = 3)
    β_val[:,Not(R_sample)] .= mean(β_val[:, R_sample], dims = 2)
    α_val[:,Not(R_sample)] .= mean(α_val[:, R_sample], dims = 2)
    if !(λ_val_tilde isa Nothing)
        λ_val = zeros(n, nC, nR)
        λ_val[:,:, R_sample] = λ_val_tilde[:,:,:]
        λ_val[:,:,Not(R_sample)] .= mean(λ_val[:,:, R_sample], dims = 3)
    end
    return α_val,β_val,p_val,λ_val
end


function kelleyPrimal(cmcmd_prob, config; ε=1e-6, kelleyPrimalEpochs=20, method = "scp_slim", rootCutsLim=kelleyPrimalEpochs, 
                        stabilizationPoint=copy(cmcmd_prob.lb), clusterInfo = clusterInfo)
    A, b, c, d, u, k, γ, lb, ub, n, nR, n_clusters, R_div_kelley = 
            cmcmd_prob.A, cmcmd_prob.b, cmcmd_prob.c, cmcmd_prob.d, cmcmd_prob.u, cmcmd_prob.k, cmcmd_prob.γ, cmcmd_prob.lb, 
            cmcmd_prob.ub, cmcmd_prob.n_edges, cmcmd_prob.n_days, cmcmd_prob.n_clusters, cmcmd_prob.sampling_rate_kelley
    cluster_partition = clusterInfo.cluster_partition
    global cutPool = PrimalSolution[]
    global cutPoolInfeas = InfeasibleCut[]
    zPath = []
    z_memory = zeros(n)
  
    offset = 0
    slope = zeros(n)
    offsets = zeros(nR)
    slopes = zeros(nR, n)
    if method == "scp_kcut"
        slopes = zeros(n_clusters, n)
    end
  
    #Root node
    TIME_LIMIT_KELLEY = TIME_LIMIT * (1/2)
    rootmodel = Model(optimizer_with_attributes(Gurobi.Optimizer, 
                "TimeLimit" => TIME_LIMIT_KELLEY,
                "OutputFlag" => 0)
            )
      
    start_time_kelley = time()
  
    @variable(rootmodel, ξ[1:n])
    # Add constraints
    println("--- [Kelley] using $(config.use_z0 ? "z0" : "0") as lower bound ---")
    lower_bound = config.use_z0 ? lb : zeros(n)
    stabilizationPoint = config.use_z0 ? stabilizationPoint : copy(cmcmd_prob.ub)
    @constraint(rootmodel, [e=1:n], lower_bound[e] <= ξ[e])
    @constraint(rootmodel, [e=1:n], ξ[e] <= ub[e])
    @constraint(rootmodel, sum(ξ) <= k)
    JuMP.set_start_value.(ξ, stabilizationPoint)


    UB = Inf; LB = -Inf
    if method == "scp_slim"
        @variable(rootmodel,θ >= 0)
        @objective(rootmodel, Min, dot(c,ξ) + θ)
        offset, slope  = cMCMDCutting_plane(cmcmd_prob, config, stabilizationPoint, adjustOffset = true)
        @constraint(rootmodel, θ >= offset + dot(slope, ξ))
        UB = dot(c,stabilizationPoint) + offset + dot(slope,stabilizationPoint[1:n])
    elseif method == "scp_fat" || method == "scp_hybrid" # hybrid is scp_fat at root node
        @variable(rootmodel, θ[1:nR] >= 0)
        @objective(rootmodel, Min, (dot(c[1:n],ξ)) + sum(θ))
        offset, slope, _, _ = cMCMDCutting_plane(cmcmd_prob, config, stabilizationPoint, adjustOffset=true, returnObjs=true)
        for r = 1:nR
            @constraint(rootmodel, θ[r] >= offset[r] + dot(slope[r,:], ξ))
        end
        UB = dot(c,stabilizationPoint) + sum(offset[r] + dot(slope[r,:],stabilizationPoint[1:n]) for r in 1:nR)
    elseif method == "scp_kcut"
        @variable(rootmodel,θ[1:n_clusters] >= 0)
        @objective(rootmodel, Min, (dot(c[1:n],ξ)) + sum(θ))
        offset = zeros(n_clusters)
        slope = zeros(n_clusters, n)
        for cluster = 1:n_clusters
            offset[cluster], slope[cluster,:]  = cMCMDCutting_plane(cmcmd_prob, config, stabilizationPoint, b_override=b[:,:,cluster_partition[cluster]], adjustOffset = true)
            @constraint(rootmodel, θ[cluster] >= offset[cluster] + dot(slope[cluster,:], ξ))
        end
        UB = dot(c,stabilizationPoint) + sum(offset[cluster] + dot(slope[cluster,:],stabilizationPoint[1:n]) for cluster in 1:n_clusters)
    end
  
    λ = .2
    δ = 2*ε
    consecutiveNonImprov_1 = 0
    consecutiveNonImprov_2 = 0
    start_time = time()
    rootCutsSense = 1
    println("doing $(kelleyPrimalEpochs) epochs")
    for epoch in 1:kelleyPrimalEpochs
        optimize!(rootmodel)
        status = termination_status(rootmodel)
        if status != MOI.OPTIMAL
            println("[Kelley] Epoch NOT MOI.OPTIMAL\n\tStatus = $(status)\n\tRoot node gap: $(abs(UB-LB))")
            kelleyPrimalEpochs = epoch
            break
        end
                      
        println("[Kelley] Termination Status: $(status)")
        zstar = zeros(length(stabilizationPoint))
        zstar[1:n] .= value.(ξ)
        θstar = value.(θ)
        z_memory[:] .= zstar[:]
        stabilizationPoint .+= zstar; stabilizationPoint ./= 2
  
        push!(zPath, zstar) #0 > LB - obj > - delta
        if (LB - objective_value(rootmodel))/abs(LB) >= -1e-5
            if consecutiveNonImprov_1 == 5
                consecutiveNonImprov_2 += 1
            else
                consecutiveNonImprov_1 += 1
            end
        else
            if consecutiveNonImprov_1 < 5
                consecutiveNonImprov_1 = 0
            end
            consecutiveNonImprov_2 = 0
        end
        LB = max(LB, objective_value(rootmodel))
  
        if consecutiveNonImprov_1 == 5
            λ = 1
        elseif consecutiveNonImprov_2 == 5
            δ = 0.
        end
  
        z0 = (epoch % kelleyPrimalEpochs == 1) ? zstar[:] : λ*zstar .+ (1-λ)*stabilizationPoint .+ δ
        z0 .= min.(1., max.(z0, 0.)) #ensure between 0 and 1
    
        # SAMPLING
        global R_sample, C_samples
        R_sample = sample(1:nR, Int(ceil(nR/R_div_kelley)), replace=false)
  
        add_infeasibility_cut = false
        if method == "scp_slim"
            z0[z0 .< 1e-4] .= 0
            offset, slope  = cMCMDCutting_plane(cmcmd_prob, config, z0, R_sample = R_sample, adjustOffset = true, use_si_vi = config.use_si_vi)
            if config.use_si_vi && isnan(offset)
                println(" *** NaN offset with SI VI *** ")
                offset, slope  = cMCMDCutting_plane(cmcmd_prob, config, z0, R_sample = R_sample, adjustOffset = true)
            end

            add_infeasibility_cut = any(isnan, offset)
            if !add_infeasibility_cut # optimality cut
                @constraint(rootmodel, θ >= offset + dot(slope,ξ)) # adjusted offset, no need for z0
                if epoch % kelleyPrimalEpochs >= 1
                    UB = min(UB, dot(c, z0) + offset + dot(slope, z0[1:n]))
                end
            end
        elseif method == "scp_fat" || method == "scp_hybrid"
            z0[z0 .< 1e-4] .= 0
            offsets, slopes = cMCMDCutting_plane(cmcmd_prob, config, z0, R_sample = R_sample, adjustOffset = true, returnObjs = true)
            add_infeasibility_cut = any(isnan, offsets)
            if !add_infeasibility_cut
                for r in R_sample
                    @constraint(rootmodel, θ[r] >= offsets[r] + dot(slopes[r,:], ξ))
                end        
                if epoch % kelleyPrimalEpochs >= 1
                    UB = min(UB, dot(c,z0) + sum(θstar[Not(R_sample)]) + sum(offsets[r] + dot(slopes[r,:],z0[1:n]) for r in R_sample))
                end
            end
        elseif method == "scp_kcut"
            z0[z0 .< 2e-4] .= 0
            C_samples = [findall(y -> y in R_sample, x) for x in cluster_partition]
            C_active = findall(x -> length(x) > 0, C_samples)
            for cluster in C_active
                offsets[cluster], slopes[cluster,:]  = cMCMDCutting_plane(cmcmd_prob, config, z0, b_override=b[:,:,cluster_partition[cluster]], R_sample = C_samples[cluster], adjustOffset = true)
                @constraint(rootmodel, θ[cluster] >= offsets[cluster] + dot(slopes[cluster,:], ξ))
            end
            if epoch % kelleyPrimalEpochs >= 1
                UB = min(UB, dot(c,z0) + sum(θstar[Not(C_active)]) + sum(offsets[cluster] + dot(slopes[cluster,:],z0[1:n]) for cluster in C_active))
            end
        end

        if add_infeasibility_cut # infeasibility cut
            println("[kelley] Adding infeasibility cut")
            p_certificate, β_certificate, infeas_r = get_infeas_certificate(cmcmd_prob, z0, R_sample = R_sample)
            if !(any(isnan, p_certificate) || any(isnan, β_certificate)) # check get_infeas_certificate terminated ok                
                # if no non-zero elements in p_certificate or β_certificate, add type 1 infeasibility cut
                for r in 1:size(β_certificate)[2]
                    if (any(β_certificate[:,r] .> 0) || any(p_certificate[:,:,r] .!= 0)) # add type 2 certificate constraint
                        r_tilde = R_sample[r]
                        @constraint(rootmodel, sum(p_certificate[:,:,r] .* b[:,:,r_tilde]) <= sum(ξ[e] * u[e] * β_certificate[e, r] for e in 1:n))
                    end
                end
                # SAVE INFEAS CUTS KELLEY
                infeas_cut = InfeasibleCut(z0, R_sample, p_certificate, β_certificate, infeas_r)
                push!(cutPoolInfeas, infeas_cut)
            else # add type 1 infeasibility cut
                throw(ArgumentError("[Kelley] ERROR"))
            end
        end

        println("\t[kelley] rootCutsSense = $(rootCutsSense), epoch = $(epoch), rootCutsLim = $(rootCutsLim), kelleyPrimalEpochs = $(kelleyPrimalEpochs), offset = $(offset)")
        if (( (rootCutsSense > 0)&(epoch <= rootCutsLim) )||( (rootCutsSense < 0)&(epoch >= kelleyPrimalEpochs-rootCutsLim+1) )) && (!add_infeasibility_cut)
            if method == "scp_slim"
                s0 = PrimalSolution( findall(z0 .> 0), z0, offset + dot(slope,z0[1:n]) + dot(c,z0), [offset], slope, false, method, R_sample, [], Dict())
            elseif method == "scp_fat" || method == "scp_hybrid"
                s0 = PrimalSolution( findall(z0 .> 0), z0, sum(offsets[r] + dot(slopes[r,:],z0[1:n]) for r in 1:nR) + dot(c,z0), offsets, slopes, false, method, R_sample, [], Dict())
            elseif method == "scp_kcut"
                s0 = PrimalSolution( findall(z0 .> 0), z0, sum(offsets[cluster] + dot(slopes[cluster,:],z0[1:n]) for cluster in 1:n_clusters) + dot(c,z0), offsets, slopes, false, method, R_sample, C_samples, Dict())
            end
            push!(cutPool, s0)
        end
  
        if abs(UB-LB)/abs(UB) <= 1e-6 || consecutiveNonImprov_2 >= 10
            println("\t[kelley] Breaking with Root node gap: ", abs(UB-LB))
            kelleyPrimalEpochs = epoch
            break
        end
  
        ##############################################
        time_spent_kelley = time() - start_time_kelley
        time_remaining_kelley = TIME_LIMIT_KELLEY - time_spent_kelley
        println("(Kelley) TIME LIMIT: $(TIME_LIMIT_KELLEY) - TIME SPENT: $(time_spent_kelley) - TIME REMAINING: $(time_remaining_kelley)")
        if time_remaining_kelley <= 0.05*TIME_LIMIT_KELLEY
            println("(Kelley) BREAKING FOR TIME: $(time_remaining_kelley)S REMAINING")
            break
        else
            set_optimizer_attribute(rootmodel, "TimeLimit", time_remaining_kelley)
        end
        ##############################################
  
  
    end # --- end for epoch in epochs
  
    solve_time = time() - start_time  
    println("\t[kelley] Kelley Primal: Generated $kelleyPrimalEpochs cuts\n\t\tBest upper bound: $UB\n\t\tBest lower bound: $LB")
    return min.(max.(z_memory, lb), ub), cutPool, cutPoolInfeas, zPath, LB, UB, solve_time
end


function add_vi_constraints(m1, cmcmd_prob, z; use_vi_1=true, use_vi_2=true, use_nci=true, print_level=0)
    b = cmcmd_prob.b
    nC = cmcmd_prob.n_commodities
    nR = cmcmd_prob.n_days
    m = cmcmd_prob.n_nodes
    edge_map = cmcmd_prob.edge_map
    outgoing_edges = cmcmd_prob.outgoing_edges
    incoming_edges = cmcmd_prob.incoming_edges
    u = cmcmd_prob.u

    # VI 1: ORIGIN/DESTINATION INEQUALITIES (OD)
    vi_1_outgoing_added = 0
    vi_1_incoming_added = 0
    if use_vi_1
        origins = [findall(x -> x > 0, b[:,i,1]) for i in 1:nC]
        destinations = [findall(x -> x < 0, b[:,i,1]) for i in 1:nC]
        for k in 1:nC
            println("ADDING (OD) VIs FOR COMMODITY k = $(k): \n  Origin Nodes: $(origins[k]) \n  Destination Nodes: $(destinations[k])")
            for origin_node in origins[k]
                if haskey(outgoing_edges, origin_node)
                    @constraint(m1, sum(z[e] for e in outgoing_edges[origin_node]) >= 1)
                    vi_1_outgoing_added += 1
                end
            end
            for destination_node in destinations[k]
                if haskey(incoming_edges, destination_node)
                    @constraint(m1, sum(z[e] for e in incoming_edges[destination_node]) >= 1)
                    vi_1_incoming_added += 1
                end
            end
        end
    end

    # VI 2: NETWORK CONNECTIVITY CUTS
    vi_2_outgoing_added = 0
    vi_2_incoming_added = 0
    if use_vi_2
        intermediary_nodes = [i for i in 1:m if all(b[i,k,r] == 0 for k in 1:nC, r in 1:nR)]
        println("Intermediary Nodes: $(intermediary_nodes)")
        for i in intermediary_nodes
            if haskey(incoming_edges, i) && haskey(outgoing_edges, i)
                # Constraint for incoming edges
                println("For intermediary node i = $(i), adding Network Connectivity constraints:\n Incoming Edges: ")
                for e_minus in incoming_edges[i]
                    println("  $(z[e_minus]) <= $(sum(z[e_plus] for e_plus in outgoing_edges[i]))")
                    @constraint(m1, z[e_minus] <= sum(z[e_plus] for e_plus in outgoing_edges[i]))
                    vi_2_incoming_added += 1
                end
                
                # Constraint for outgoing edges
                println(" Outgoing Edges: ")
                for e_plus in outgoing_edges[i]
                    println("  $(z[e_plus]) <= $(sum(z[e_minus] for e_minus in incoming_edges[i]))")
                    @constraint(m1, z[e_plus] <= sum(z[e_minus] for e_minus in incoming_edges[i]))
                    vi_2_outgoing_added += 1
                end
            end
        end
    end

    # Cardinality network connectivity inequalities (NCIs) -- redundant with VI 1 if h = 0
    nci_outgoing_added = 0
    nci_incoming_added = 0
    if use_nci
        for i in 1:m 
            dmax_out = maximum([sum([b[i,k,r] for k in 1:nC if b[i,k,r] > 0]) for r in 1:nR])
            if dmax_out > 0 && haskey(outgoing_edges, i)
                out_edges = outgoing_edges[i]
                order_edges = sortperm([u[e] for e in out_edges], rev=true)
                h = 0
                while h < length(out_edges) && sum(u[out_edges[order_edges[j]]] for j in 1:(h+1)) < dmax_out
                    h += 1
                end
                if h > 0
                    @constraint(m1, sum(z[e] for e in out_edges) >= h+1)
                    nci_outgoing_added += 1
                end 
            end

            dmax_in = maximum([sum([-b[i,k,r] for k in 1:nC if b[i,k,r] < 0]) for r in 1:nR])
            if dmax_in > 0 && haskey(incoming_edges, i)
                in_edges = incoming_edges[i]
                order_edges = sortperm([u[e] for e in in_edges], rev=true)
                h = 0
                while h < length(in_edges) && sum(u[in_edges[order_edges[j]]] for j in 1:(h+1)) < dmax_in
                    h += 1
                end
                if h > 0
                    @constraint(m1, sum(z[e] for e in in_edges) >= h+1)
                    nci_incoming_added += 1
                end    
            end
        end
    end
    print_level > 0 && println("[add_vi_constraints] NUMBER OF CONS ADDED:\n\t\t $(vi_1_outgoing_added) outgoing VI 1 constraints\n\t\t $(vi_1_incoming_added) incoming VI 1 constraints\n\t\t $(vi_2_outgoing_added) outgoing VI 2 constraints\n\t\t $(vi_2_incoming_added) incoming VI 2 constraints\n\t\t $(nci_outgoing_added) outgoing NCI constraints\n\t\t $(nci_incoming_added) incoming NCI constraints")
end
