using MaximumEntropyMomentClosures
using Test, SafeTestsets

@time @safetestset "Quadrature" begin include("quadrature_test.jl") end
@time @safetestset "Moments" begin include("moments_test.jl") end
@time @safetestset "LagrangeMultipliers" begin include("lagrange_multipliers_test.jl") end
