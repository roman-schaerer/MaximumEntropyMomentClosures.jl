abstract type Quadrature{T} end

@with_kw mutable struct GaussLegendre{T} <: Quadrature{T}
    lower::T = -one(T) # lower bound of domain
    upper::T = one(T) # upper bound of domain
    x::Vector{T} = Float64[] # quadrature nodes
    w::Vector{T} = Float64[] # quadrature weights
    ϕ::Array{T} = Float64[] # basis evaluations at quadrature nodes
end

@with_kw mutable struct GaussHermite{T} <: Quadrature{T}
    μ::T = zero(T)
    σ::T = one(T)
    x::Vector{T} = Float64[]
    w::Vector{T} = Float64[]
    ϕ::Array{T} = Float64[]
end

# m = degree of highest polynomial
function init!(qr::Quadrature{T}; nq=16, m=0) where T
    quadrature_rule!(qr, nq)
    compute_monomials!(qr, m)
end

# m = degree of highest polynomial
# nq = number of quadrature points
function compute_monomials!(q::Quadrature{T}, m) where T
    nq = length(q.x)
    q.ϕ = Matrix{T}(undef, m+1, nq)
    q.ϕ[1,:] = ones(nq)
    for j in range(2,m+1,step=1)
        q.ϕ[j,:] = q.ϕ[j-1,:].*q.x
    end
    return Nothing
end


# nq = number of quadrature nodes
# m = number of moments in the function ansatz
# p = number of moments to be computed
# nq = number of quadrature points
function quadrature_rule!(gl::GaussLegendre{T}, nq) where T
    gl.x, gl.w = gausslegendre( nq )
    @. gl.x = (gl.upper + gl.lower)/2 + (gl.x) * (gl.upper - gl.lower)/2
    @. gl.w *= (gl.upper - gl.lower)/2
end

function quadrature_rule!(gh::GaussHermite{T}, nq) where T
    gh.x, gh.w = gausshermite( nq )
    @. gh.x = gh.x * gh.σ + gh.μ
    @. gh.w *= gh.σ # CHECK THIS
end
