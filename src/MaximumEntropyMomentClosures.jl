module MaximumEntropyMomentClosures

using LinearAlgebra
using LoopVectorization
using Tullio
using FastGaussQuadrature
using Parameters
using Optimization, OptimizationOptimJL

# export types
export GaussLegendre, GaussHermite
export NewtonParameters

# export functions
export moments!, hessian!, lagrange_multipliers!
export central_normalized_moments!, raw_moments!
export init!

include("base.jl")
include("utils.jl")
include("quadrature.jl")
include("moments.jl")
include("lagrange_multipliers.jl")

end
