"""
    NewtonParameters{FloatType, IntType}

Store parameters of the Newton search algorithm
"""
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


"""
    lagrange_multipliers!

Determine lagrange parameters α using a Newton-based search algorithm with an adaptive change of basis.
"""
function lagrange_multipliers!(α::AbstractVector{T}, u::AbstractVector{T}, tmat::AbstractArray{T}, 
                 g::AbstractVector{T}, h::AbstractArray{T}, d::AbstractVector{T},
                 qr::Quadrature{T}, params::NewtonParameters{T, IntType}) where {T, IntType}

    stats = NewtonStatistics{Int64}()
    αnew = similar(α)

    f = zeroth_moment_cartesian_basis(α, qr) - dot(α, u)

    for iter in 1:params.max_iter
        if change_basis!(α, u, tmat, h, qr) != Success
            println("error: change_basis failed") 
            return (Failure, stats)
        end
        # evaluate gradient
        fill!(g, 0)
        g[1] = zeroth_moment(α, qr)
        g .-= u
        
        # set descent direction
        d .= -g

        # evaluate mismatch in the moments
        u_tmp = similar(u)
        moments!(u_tmp, α, qr)
        err = norm(u_tmp - u)

        # check stopping criteria
        if err < params.τ && exp(5*norm(tmat \ d, 1)) < 1 + params.ϵγ
            return (Success, stats)
        else
            line_search_iter = 0
            ζ = one(T)
            while params.ξ * norm(d) > params.ϵ * norm(α)
                αnew .= α + ζ * d
                fnew = zeroth_moment_cartesian_basis(αnew, qr) - dot(αnew, u)
                stats.num_iter_line_search += 1
                if fnew <= f + params.ξ * ζ * dot(g, d)
                    α .= αnew
                    f = fnew
                    stats.num_iter_newton += 1
                    break
                end
                ζ = params.χ * ζ
                line_search_iter += 1
                if line_search_iter > params.max_iter_line_search
                    return (Failure, stats)
                end
            end
        end
    end
    return (Failure, stats)
end

"""
    change_basis!

Determine orthogonal basis that diagonalizes the hessian matrix. Update all input parameters to represent values in the new basis.
"""
function change_basis!(αvec::AbstractVector{T}, uvec::AbstractVector{T}, 
    tmat::AbstractArray{T}, hmat::AbstractArray{T}, qr::Quadrature{T}) where {T}

    # evaluate hessian in current basis
    hessian!(hmat, αvec, qr)

    chol = cholesky(hmat, check=true) # TODO: should be pivot here?
    # update basis vectors
    qr.ϕ .= chol.L \ qr.ϕ
    # update matrix t (transpose of change-of-basis matrix)
    tmat .= tmat * chol.L
    # update lagrange parameters
    αvec .= transpose(chol.L) * αvec
    # update moment vector
    uvec .= chol.L \ uvec
    return Success
end

function scaling_factor(m)
    double_factorial(2 * m - 1)^(1 / (2 * m))
end

function scale_basis!(u::AbstractVector{T}, α::AbstractVector{T}) where {T}
    m = length(u) - 1
    # define scaling parameter based on the highest moment of the standard normal distribution
    s = scaling_factor(m)

    for idx in range(2:m+1)
        u[idx] /= s^(idx - 1)
    end
    α .*= s
end


"""
    normal_lagrange_multipliers!(α::AbstractVector{T}, uc::AbstractVector{T}) where T

Determine the Lagrange multipliers vector α of a normal distribution over ``\\mathbb{R}`` for a given vector of central moments uc
"""
function normal_lagrange_multipliers!(α::AbstractVector{T}, uc::AbstractVector{T}) where T
    α[1] = -uc[2]^2/(uc[3]*2) + log(uc[1]/(sqrt(2*π*uc[3])))
    α[2] = uc[2]/uc[3]
    α[3] = -1.0/(2*uc[3])
    return α
end

"""
    init_lagrange_multipliers!(α_init::AbstractArray{T}, u::AbstractArray{T}) where T

Initialize Lagrange multipliers α_init for a given array of raw moments u. 
"""
function init_lagrange_multipliers!(α_init::AbstractArray{T}, u::AbstractArray{T}) where T
    uc = Vector{T}(undef, 3)
    α = Vector{T}(undef, 3)
    fill!(α_init,zero(T))
    idx = 1
    for uvec in eachcol(u)
        central_normalized_moments!(uc, uvec[1:3])
        normal_lagrange_multipliers!(α, uc)
        α_init[1:3, idx] = α
        idx += 1
    end
end

function init_lagrange_multipliers!(α_init::AbstractVector{T}, u::AbstractVector{T}) where T
    @assert length(α_init) >= 3 # current implementation is restricted to m = 2
    uc = Vector{T}(undef, 3)
    α = Vector{T}(undef, 3)
    fill!(α_init, zero(T))
    central_normalized_moments!(uc, u[1:3])
    normal_lagrange_multipliers!(α, uc)
    α_init[1:3] = α
end

"""
lagrange_multipliers!(α::AbstractArray{T}, u::AbstractArray{T}, p) where T

Determine the Lagrange multipliers for each moment vector stored as a column vector in the abstract array u.
"""
function lagrange_multipliers!(α::AbstractArray{T}, u::AbstractArray{T}, qr::Quadrature{T}) where T
    αvec = Vector{T}(undef, size(α)[1])
    success = true
    for idx in 1:size(u)[2]
        success &= lagrange_multipliers!(αvec, view(u,:,idx), qr) 
        α[:,idx] = αvec
    end
    return success
end

function lagrange_multipliers!(α::AbstractVector{T}, u::AbstractVector{T}, qr::Quadrature{T}) where T
    optim_fun = OptimizationFunction(dual_function; grad=dual_gradient!, hess=hessian!)
    init_lagrange_multipliers!(α, u)
    params = (qr=qr, u=u)
    prob = OptimizationProblem(optim_fun, α, params)
    sol = solve(prob, Newton())
    α[:] = sol
    return sol.original.ls_success
end

function dual_function(α::AbstractVector{T}, p) where T
    u0 = zeroth_moment(α, p.qr)
    f = u0 - dot(α, p.u)
    return f
end

function dual_function(α::AbstractArray{T}, p) where T
    u0 = zeroth_moment(α, p.qr)
    f = Vector{T}(undef, size(α)[2])
    @tullio f[j] = u0[j] - α[i,j]*p.u[i,j]
    return f
end

function dual_gradient!(g::AbstractArray{T}, α::AbstractArray{T}, p) where T
    moments!(g, α, p.qr)
    g[:] -= p.u[1:size(α)[1]]
    return Nothing
end

