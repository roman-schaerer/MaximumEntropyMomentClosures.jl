abstract type MomentEvaluation end
abstract type ExactIntegration{T} <: MomentEvaluation end
abstract type Quadrature{T} <: MomentEvaluation end
