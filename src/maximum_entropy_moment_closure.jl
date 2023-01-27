

mutable struct AdaptiveBasisEvaluation{T} <: LagrangeParameterEvaluation where T
    hessian::Array{T} # hessian matrix
    tmat::Array{T} # change of basis matrix
    g::Vector{T} # gradient vector
    d::Vector{T} # descent vector
    params::NewtonParameters{T, Int64} # parameters of the Newton search algorithm

    function AdaptiveBasisEvaluation{T}(n=0) where T
        hessian = Array{T,2}(undef, n, n)
        tmat = Array{T,2}(I, n, n) # initialize as identity matrix
        g = Vector{T}(undef, n)
        d = Vector{T}(undef, n)
        params = NewtonParameters{T, Int64}()
        new{T}(hessian, tmat, g, d, params)
    end
end

mutable struct MaximumEntropyMomentClosure{MS<:MomentSystem, 
                                           ME<:MomentEvaluation, 
                                           LE<:LagrangeParameterEvaluation}
    ms::MS
    me::ME
    le::LE

    function MaximumEntropyMomentClosure(ms::MS, me::ME, le::LE) where 
                                         {MS<:MomentSystem, 
                                          ME<:MomentEvaluation,
                                          LE<:LagrangeParameterEvaluation}
        new{MS, ME, LE}(ms, me, le)
    end
end

