using MaximumEntropyMomentClosures: init!, 
                                    moments!, 
                                    GaussLegendre, 
                                    GaussHermite, 
                                    change_basis!
using Statistics
using Test

function gauss_legendre_test(lower=0.5, upper=6.0, nq=16, m=2)
    gl = GaussLegendre{Float64}(lower=lower, upper=upper)
    init!(gl, nq=nq, m=m)

    @test gl.lower == lower && gl.upper == upper
    @test all(gl.x .>= gl.lower) && all(gl.x .<= gl.upper)
    @test sum(gl.w) ≈ upper-lower
end

function gauss_hermite_test(μ=2.5,σ=0.25,nq=16,m=2)
    gh = GaussHermite{Float64}(μ=μ, σ=σ)
    init!(gh, nq=nq, m=m)

    @test sum(gh.w) ≈ sqrt(π)*σ
    @test mean(gh.x) ≈ μ
end

gauss_legendre_test()
gauss_hermite_test()