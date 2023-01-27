function double_factorial(n::Integer)
    p = 1
    for i in range(1+iseven(n),n,step=2)
        p *= i
    end
    return p
end

# monomial basis
function init_gl_mb_1d(;lower=lower, upper=uppper, degree=degree, nq=nq)
    ms = MomentSystemFullTensor(Interval(lower=lower, upper=upper), degree)
    gl = GaussLegendre(ms, nq)
    le = MonomialBasisEvaluation()
    memc = MaximumEntropyMomentClosure(ms, gl, le)
    return memc
end

# adaptive basis
function init_gl_ab_1d(;lower=lower, upper=uppper, degree=degree, nq=nq)
    ms = MomentSystemFullTensor(Interval(lower=lower, upper=upper), degree)
    gl = GaussLegendre(ms, nq)
    le = AdaptiveBasisEvaluation{Float64}(degree+1)
    memc = MaximumEntropyMomentClosure(ms, gl, le)
    return memc
end