# using LightGraphs,
import Pkg
# Pkg.activate(".")
using DataFrames, JuMP,Gurobi, StatsBase, LinearAlgebra, Distributions, Random, TimerOutputs
using Plots, ArgParse, JLD2, Clustering
using Suppressor
using MosekTools, Mosek

const TYPES_DET = ["naive", "det", "det_slim"]
const TYPES_SCP = ["scp", "scp_avg", "scp_slim", "scp_slim_feas"]
global TIME_LIMIT = 3600
# global TIME_LIMIT = 72000
# const VER_CON = "LONG_SPARSE_NAIVE"
# const VER_CON = "SMPL"
# const VER_CON = "EDGE"
const VER_CON = "BNDS"

println("Load Functions ...")

include("types.jl")
include("general_library.jl")
include("scp_library.jl")