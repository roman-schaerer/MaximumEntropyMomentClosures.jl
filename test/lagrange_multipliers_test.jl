using MaximumEntropyMomentClosures: 
                                    MaximumEntropyMomentClosures,
                                    init!, 
                                    moments!, 
                                    hessian!,
                                    GaussLegendre, 
                                    change_basis!,
                                    lagrange_multipliers!,
                                    NewtonParameters,
                                    lagrange_multipliers!
using Test
using BenchmarkTools
using Random
using LinearAlgebra

Random.seed!(123)


function lagrange_vec_test(;lower=-1.0, upper=1.0, nq=64, m=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)

    αvec = randn(m+1)
    uvec = similar(αvec)
    moments!(uvec, αvec, gl)
    αvec_res = similar(αvec)
    @test lagrange_multipliers!(αvec_res, uvec, gl)
    @test isapprox(αvec, αvec_res, rtol=1e-6)
end

function lagrange_mat_test(;lower=-1.0, upper=1.0, nq=64, m=2, n=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)

    αmat = randn(m+1, n)
    umat = similar(αmat)
    moments!(umat, αmat, gl)
    αmat_res = similar(αmat)
    @test lagrange_multipliers!(αmat_res, umat, gl)
    @test isapprox(αmat, αmat_res, rtol=1e-6)
end




function change_basis_test(;lower=-1.0, upper=1.0, nq=64, m=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)

    αvec = randn(m+1)
    uvec = similar(αvec)
    moments!(uvec, αvec, gl)

    tmat = Matrix{Float64}(I, m+1, m+1)
    hmat = similar(tmat)
    hessian!(hmat, αvec, gl)

    uvec_orig = copy(uvec)

    change_basis!(αvec, uvec, tmat, hmat, gl)

    uvec_check = similar(αvec)
    moments!(uvec_check, αvec, gl)

    @test isapprox(uvec_check, uvec)
    @test isapprox(uvec_orig, tmat*uvec)

    hmat_check = similar(hmat)
    hessian!(hmat_check, αvec, gl)
    @test isapprox(Matrix{Float64}(I, m+1, m+1), hmat_check, rtol=1e-8)
end

function lagrange_multipliers_test(;lower=-1.0, upper=1.0, nq=64, m=2)
    gl = GaussLegendre(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)

    αvec_ref = randn(m+1)
    αvec_ref[m+1] = -abs(αvec_ref[m+1])
    uvec_ref = similar(αvec_ref)
    moments!(uvec_ref, αvec_ref, gl)

    αvec_pert = randn(m+1)
    αvec_pert[m+1] = -abs(αvec_pert[m+1])

    αvec = αvec_ref + αvec_pert

    uvec = copy(uvec_ref)
    tmat = Matrix{Float64}(I, m+1, m+1)
    gvec = zeros(Float64, m+1)
    dvec = Vector{Float64}(undef, m+1)
    hmat = Matrix{Float64}(undef, m+1, m+1)

    newton_params = NewtonParameters{Float64, Int64}(max_iter=40)
    result = lagrange_multipliers!(αvec, uvec, tmat, gvec, hmat, dvec, gl, newton_params)
    @test result[1] == MaximumEntropyMomentClosures.Success

    uvec_res = similar(uvec_ref)
    moments!(uvec_res, αvec, gl)

    uvec_res .= tmat * uvec_res # transform to standard basis
    αvec_res = copy(transpose(tmat) \ αvec) # transform to standard basis

    @test isapprox(αvec_res, αvec_ref, rtol=1e-6)
    @test isapprox(uvec_res, uvec_ref, rtol=1e-6)

    # initialize αvec_ref as a perturbation of αvec
    αvec_pert = randn(m+1)
    αvec_pert[m+1] = -abs(αvec_pert[m+1])
    αvec_ref = αvec + αvec_pert

    moments!(uvec_ref, αvec_ref, gl)
    uvec .= uvec_ref
    tmat_ref = copy(tmat)

    result = lagrange_multipliers!(αvec, uvec, tmat, gvec, hmat, dvec, gl, newton_params)
    @test result[1] == MaximumEntropyMomentClosures.Success

    uvec_res = similar(uvec)
    moments!(uvec_res, αvec, gl)

    uvec_res = tmat * uvec_res # transform to standard basis
    αvec_res = transpose(tmat) \ αvec # transform to standard basis

    αvec_ref = transpose(tmat_ref) \ αvec_ref # transform to standard basis
    uvec_ref = tmat_ref*uvec_ref # transform to standard basis

    @test isapprox(αvec_res, αvec_ref, rtol=1e-6)
    @test isapprox(uvec_res, uvec_ref, rtol=1e-6)
end

lagrange_vec_test()
lagrange_mat_test()
lagrange_mat_test(lower=-0.5, upper=1.5)

change_basis_test()
change_basis_test(nq=4, m=2)
change_basis_test(nq=8, m=4)
change_basis_test(nq=16, m=4)
change_basis_test(nq=128, m=4)
change_basis_test(nq=1024, m=8)

lagrange_multipliers_test()
lagrange_multipliers_test(m=0)
lagrange_multipliers_test(m=2)
lagrange_multipliers_test(m=4)
lagrange_multipliers_test(m=6)
lagrange_multipliers_test(m=8)

lagrange_multipliers_test(m=2,lower=-0.2,upper=0.8)
lagrange_multipliers_test(m=2,lower=-20.0,upper=20.0)

lagrange_multipliers_test(m=4,lower=-2.0,upper=1.5)
# lagrange_multipliers_test(m=4,lower=-4.0,upper=4.0)