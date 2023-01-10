using MaximumEntropyMomentClosures: init!, 
                                    moments!,
                                    hessian!,
                                    GaussLegendre, 
                                    change_basis!, 
                                    raw_moments!, 
                                    central_normalized_moments!
using Test
using BenchmarkTools
using Random
using LinearAlgebra

Random.seed!(123)

function basis_evaluation_test(;lower=-2.0,upper=3.0,nq=16,m=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)
    
    @test all(gl.ϕ[1,:] .≈ 1)
    @test all(gl.ϕ[3,:] .≈ gl.ϕ[2,:] .* gl.x)
end

function central_normalized_moments_test(;lower=-1.0,upper=1.0,nq=64,m=2)
    ucn_vec = rand(m+1)
    u_vec = similar(ucn_vec)
    raw_moments!(u_vec, ucn_vec)
    ucn_vec_res = similar(ucn_vec)
    central_normalized_moments!(ucn_vec_res, u_vec)    
    @test isapprox(ucn_vec_res, ucn_vec, rtol=1e-6)
end

function moment_vector_test(;lower=-1.0,upper=1.0,nq=64,m=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)
    
    αvec = [1.0, -0.5, -0.25]
    uvec = [5.2133447, -0.8007096, 1.67866228]
    
    uvec_res = similar(uvec)
    moments!(uvec_res, αvec, gl)
    
    @test maximum(abs.(uvec - uvec_res)) < 1e-6
end

function moment_matrix_test(;lower=-1.0,upper=1.0,nq=64,m=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)
 
    n = 2 # number of integration problems
    αmat = randn(m+1, n)
    αmat[m+1,:] .= -abs.(αmat[m+1,:]) # ensure numerical integrability

    umat = Matrix{Float64}(undef, m+1, n)
    moments!(umat, αmat, gl)

    umat2 = similar(umat)
    uvec = Vector{Float64}(undef, m+1)
    for idx_col in 1:n
        moments!(uvec, αmat[:,idx_col], gl)
        umat2[:,idx_col] = uvec
    end
    @test isapprox(umat, umat2)
end

function hessian_test(;lower=-1.0, upper=1.0, m=2, nq=64)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)

    αvec = randn(m)
    hess = Matrix{Float64}(undef,m,m)
    hessian!(hess, αvec, gl)

    # check that hess is symmetric and positive definite
    @test LinearAlgebra.isposdef(hess) && LinearAlgebra.issymmetric(hess)
end

basis_evaluation_test()
moment_vector_test()
moment_matrix_test()
central_normalized_moments_test()
hessian_test()
hessian_test(m=4)
hessian_test(m=6)
hessian_test(m=8)