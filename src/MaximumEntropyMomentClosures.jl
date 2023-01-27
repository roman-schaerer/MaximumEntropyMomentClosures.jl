module MaximumEntropyMomentClosures

using LinearAlgebra
using LoopVectorization
using Tullio
using FastGaussQuadrature
using Parameters
using Optimization
using OptimizationOptimJL

using Infiltrator #FIXME: Remove


# export types
export Interval
export MomentSystemFullTensor
export MonomialBasisEvaluation, AdaptiveBasisEvaluation
export MaximumEntropyMomentClosures
export GaussLegendre, GaussHermite
export NewtonParameters

# export functions
export moments!, hessian!, lagrange_multipliers!
export central_normalized_moments!, raw_moments!
export init!

include("base.jl")
include("domain.jl")
include("moment_system.jl")
include("moment_evaluation.jl")
include("quadrature.jl")
include("exact_integration.jl")
include("maximum_entropy_moment_closure.jl")
include("utils.jl")
include("moments.jl")
include("lagrange_multipliers.jl")

end
