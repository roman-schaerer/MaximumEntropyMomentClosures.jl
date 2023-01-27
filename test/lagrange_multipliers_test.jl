using MaximumEntropyMomentClosures
using MaximumEntropyMomentClosures: MaximumEntropyMomentClosure,
                                    init!, 
                                    moments!, 
                                    hessian!,
                                    GaussLegendre, 
                                    change_basis!,
                                    lagrange_multipliers!,
                                    NewtonParameters,
                                    lagrange_multipliers!
using Test
using Random
using LinearAlgebra

Random.seed!(123)

include("../src/utils.jl")

function lagrange_vec_test(memc::MaximumEntropyMomentClosure{MS, ME, LE}) where {MS, ME, LE}
    αvec = randn(memc.ms.degree+1)
    αvec[end] = -abs(αvec[end])

    uvec = similar(αvec)
    moments!(uvec, αvec, memc.me)
    αvec_test = randn(memc.ms.degree+1)
    αvec_test[end] = -abs(αvec_test[end])

    @test lagrange_multipliers!(αvec_test, uvec, memc)
    @test isapprox(αvec_test, αvec, rtol=1e-6)

    uvec_test = similar(αvec)
    moments!(uvec_test, αvec_test, memc.me)
    @test isapprox(uvec_test, uvec, rtol=1e-6)

    αvec_pert = αvec * (1.0 + 1e-2)
    uvec_pert = similar(αvec_pert)
    moments!(uvec_pert, αvec_pert, memc.me)

    αvec_pert_test = αvec
    @test lagrange_multipliers!(αvec_pert_test, uvec_pert, memc)
    @test isapprox(αvec_pert_test, αvec_pert, rtol=1e-6)
end

function lagrange_mat_test(memc::MaximumEntropyMomentClosure{MS, ME, LE}, n) where {MS, ME, LE}
    αmat = randn(memc.ms.degree+1, n)
    αmat[end,:] .= -abs.(αmat[end,:])

    umat = similar(αmat)
    moments!(umat, αmat, memc.me)
    αmat_test = randn(memc.ms.degree+1, n)
    αmat_test[end,:] .= -abs.(αmat_test[end,:])

    @test lagrange_multipliers!(αmat_test, umat, memc.me)
    @test isapprox(αmat_test, αmat, rtol=1e-6)

    umat_test = similar(αmat)
    moments!(umat_test, αmat_test, memc.me)
    @test isapprox(umat_test, umat, rtol=1e-6)   
end


function change_basis_test(memc::MaximumEntropyMomentClosure{MS, ME, LE}) where {MS, ME, LE}
    degree = memc.ms.degree

    αvec = randn(degree+1)
    uvec = similar(αvec)
    moments!(uvec, αvec, memc.me)

    tmat = Matrix{Float64}(I, degree+1, degree+1)
    hmat = similar(tmat)
    hessian!(hmat, αvec, memc.me)

    uvec_orig = copy(uvec)

    change_basis!(αvec, uvec, memc)

    uvec_check = similar(αvec)
    moments!(uvec_check, αvec, memc.me)

    @test isapprox(uvec_check, uvec)
    @test isapprox(uvec_orig, memc.le.tmat*uvec)

    hmat_check = similar(hmat)
    hessian!(hmat_check, αvec, memc.me)
    @test isapprox(Matrix{Float64}(I, degree+1, degree+1), hmat_check, rtol=1e-8)
end

lagrange_vec_test(init_gl_mb_1d(lower=-1.0, upper=1.0, degree=2, nq=64))
lagrange_vec_test(init_gl_mb_1d(lower=0.0, upper=1.0, degree=4, nq=32))
lagrange_vec_test(init_gl_mb_1d(lower=-1.0, upper=0.5, degree=6, nq=128))

lagrange_vec_test(init_gl_ab_1d(lower=-1.0, upper=3.0, degree=2, nq=64))
lagrange_vec_test(init_gl_ab_1d(lower=0.0, upper=5.0, degree=4, nq=32))
lagrange_vec_test(init_gl_ab_1d(lower=-1.0, upper=0.5, degree=6, nq=128))

lagrange_vec_test(init_gl_ab_1d(lower=-2.0, upper=2.0, degree=2, nq=128))
lagrange_vec_test(init_gl_ab_1d(lower=-0.6, upper=1.2, degree=2, nq=128))
lagrange_vec_test(init_gl_ab_1d(lower=-2.0, upper=2.0, degree=4, nq=128))
lagrange_vec_test(init_gl_ab_1d(lower=-2.0, upper=2.0, degree=6, nq=128))

lagrange_mat_test(init_gl_mb_1d(lower=-1.0, upper=1.0, degree=2, nq=64), 1)
lagrange_mat_test(init_gl_mb_1d(lower=-2.0, upper=1.1, degree=2, nq=32), 32)
lagrange_mat_test(init_gl_mb_1d(lower=-1.0, upper=1.0, degree=4, nq=32), 32)
lagrange_mat_test(init_gl_mb_1d(lower=-1.0, upper=1.0, degree=6, nq=32), 32)
lagrange_mat_test(init_gl_mb_1d(lower=-1.0, upper=1.0, degree=6, nq=256), 32)


change_basis_test(init_gl_ab_1d(lower=-1.0, upper=1.0, degree=2, nq=4))
change_basis_test(init_gl_ab_1d(lower=-1.0, upper=1.0, degree=2, nq=16))
change_basis_test(init_gl_ab_1d(lower=-1.0, upper=1.0, degree=4, nq=32))
change_basis_test(init_gl_ab_1d(lower=-1.0, upper=1.0, degree=6, nq=32))
change_basis_test(init_gl_ab_1d(lower=-1.0, upper=1.0, degree=8, nq=64))
change_basis_test(init_gl_ab_1d(lower=-1.0, upper=1.0, degree=8, nq=1024))
