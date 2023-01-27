@enum RETURN_CODES Failure Success


@with_kw mutable struct NewtonParameters{FloatType, IntType}
    τ::FloatType = sqrt(eps(FloatType)) # τ > 0
    ϵ::FloatType = eps(FloatType) # ϵ > 0
    ϵγ::FloatType = 1e-2 # ϵγ > 0
    ξ::FloatType = 1e-3 # ξ ∈ (0, 1/2)
    χ::FloatType = 0.5 # χ ∈ (0, 1)
    max_iter::IntType = 100
    max_iter_line_search::IntType = 40
end

@with_kw mutable struct NewtonStatistics{IntType}
    num_iter_newton::IntType = 0
    num_iter_line_search::IntType = 0
end

abstract type LagrangeParameterEvaluation end

mutable struct MonomialBasisEvaluation <: LagrangeParameterEvaluation
end
