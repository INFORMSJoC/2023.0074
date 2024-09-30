[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# stochastic-unified-network-design

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[A Stochastic Benders Decomposition Scheme for
Large-Scale Stochastic Network Design](https://doi.org/10.1287/ijoc.2023.0074) by Dimitris Bertsimas, Ryan Cory-Wright, Jean Pauphilet, and Periklis Petridis.


# Quick Start

### DataNetworkDesign.jl 
The main julia file for running everything. You can set all arguments from the command line. 

Example usage:
```bash
julia DataNetworkDesign.jl --time 120  --method scp_hybrid --method_kelley scp_fat --nodes 2 --commodities 5 --days 10 --sample 4 --rootCuts 20 --sample_kelley 4 --use_z0
```

### Parameters 

Explanation of parameters in example above:
- `--time 120` is the time limit in seconds. Default is 3600.
- `--method scp_hybrid` is the method to use for the stochastic cutting planes. Default is `scp_slim`.
- `--method_kelley scp_fat` is the method to use for the Kelley (Root Node) cutting planes. Default is `scp_slim`.
- `--nodes 2` is the number of nodes (multiplied by 10) in the network. Default is 1, i.e. 10 nodes.
- `--commodities 5` is the number of commodities. Default is 5.
- `--days 10` is the number of scenarios to generate. Default is 10.
- `--sample 4` is the sampling **ratio** for the stochastic cutting planes. Default is 4, i.e. sample 1/4 (25%) of scenarios/days at each iteration.
- `--rootCuts 20` is the number of root node cuts to add at the start of the algorithm. Default is 0, we usually do 20 when enabled.
- `--sample_kelley 4` is the sampling **ratio** for the Kelley cuts. Default is 4, i.e. sample 1/4 (25%) of scenarios/days at each iteration.
- `--use_z0` is a flag to to use a feasible spanning tree as a starting point. Default is false.


For more information on the parameters, see the function `parse_commandline()` in `src/general_library.jl`:

```julia
# src/general_library.jl
function parse_commandline(args)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--jobid", "--jid" 
            help = "job id corresponding to this run"
            arg_type = String
            default = "00000"
        ...
        ...
    end
    return parse_args(args, s)
end

```


# Test Run
You can also run the test file to check if everything is working properly. This runs a few different configurations of the algorithm and checks if the results are as expected.

```bash
julia tests/TestDataNetworkDesign.jl
```

This should generate something like the following output:
```bash
Compact Summary (values rounded to 3 decimals):
Method                                        | Obj Value    | True Obj     | Obj Bnd      | Gap     | True Gap | Time   | Edges  |conf_bnd_gap | outter_iters
----------------------------------------------|--------------|--------------|--------------|---------|----------|--------|--------|-------------|-------------
naive                                         |    74443.664 |    74443.535 |    74443.535 |     0.0 |      0.0 |   0.27 |   88.0 |         0.0 |         0.0
scp_fat_kelley-scp_fat_cuts-0_det             |    74443.542 |    74443.535 |    74443.542 |     0.0 |      0.0 |   4.71 |   88.0 |       0.011 |         0.0
scp_fat_kelley-scp_fat_cuts-0_stoch           |    74443.542 |    74443.535 |    74443.542 |     0.0 |      0.0 |  6.481 |   88.0 |      -0.009 |         2.0
scp_fat_kelley-scp_fat_cuts-20_det            |    74443.597 |    74443.535 |    74443.597 |     0.0 |      0.0 |  3.891 |   88.0 |       0.011 |         0.0
scp_fat_kelley-scp_fat_cuts-20_stoch          |    74443.542 |    74443.535 |    74443.542 |     0.0 |      0.0 |  9.387 |   88.0 |      -0.034 |         5.0
```

# Code Structure
Main code separated in the following files:
- `DataNetworkDesign.jl`: Main file to run the algorithm.
- `src/general_library.jl`: General functions used in the algorithm.
- `src/scp_library.jl`: Functions related to the stochastic cutting planes.


# Quick start for SLURM Environments
Use commands analogous to the following to run the code in HPC cluster environment with SLURM (interactively).
Below are the commands for MIT Engaging.

```bash
# LOAD MODULES
module load julia/1.7.3
module load sloan/gurobi/9.1.2
module load mosek/9.3

# ACTIVATE INTERACTIVE ENVIRONMENT (SRUN)
srun --pty --partition=sched_mit_sloan_interactive --cpus-per-task=1 --mem=4G bash

# RUN JULIA
julia DataNetworkDesign.jl --rootCuts 20 --sample_kelley 4 --days 10 --nodes 1 --sample 4
```


 
# Citing
To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0074

https://doi.org/10.1287/ijoc.2023.0074.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{scp-network-design,
  author =        {Dimitris Bertsimas, Ryan Cory-Wright, Jean Pauphilet, Periklis Petridis},
  publisher =     {INFORMS Journal on Computing},
  title =         {{A Stochastic Benders Decomposition Scheme for
Large-Scale Stochastic Network Design}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0074.cd},
  url =           {https://github.com/INFORMSJoC/2023.0074},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0074},
}  
```