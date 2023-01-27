
# """
#     normal_lagrange_params!(α::AbstractVector{T}, uc::AbstractVector{T}) where T

# Determine the Lagrange parameter vector α of a normal distribution over ``\\mathbb{R}`` for a given vector of central moments uc
# """
# function normal_lagrange_params!(α::AbstractVector{T}, uc::AbstractVector{T}) where T
#     α[1] = -uc[2]^2/(uc[3]*2) + log(uc[1]/(sqrt(2*π*uc[3])))
#     α[2] = uc[2]/uc[3]
#     α[3] = -1.0/(2*uc[3])
#     return α
# end

# """
#     init_lagrange_params!(α_init::AbstractArray{T}, u::AbstractArray{T}) where T

# Initialize Lagrange parameters α_init for a given array of raw moments u. 
# """
# function init_lagrange_params!(α_init::AbstractArray{T}, u::AbstractArray{T}) where T
#     uc = Vector{T}(undef, 3)
#     α = Vector{T}(undef, 3)
#     fill!(α_init,zero(T))
#     idx = 1
#     for uvec in eachcol(u)
#         central_normalized_moments!(uc, uvec[1:3])
#         normal_lagrange_params!(α, uc)
#         α_init[1:3, idx] = α
#         idx += 1
#     end
# end

# function init_lagrange_params!(α_init::AbstractVector{T}, u::AbstractVector{T}) where T
#     @assert length(α_init) >= 3 # current implementation is restricted to m = 2
#     uc = Vector{T}(undef, 3)
#     α = Vector{T}(undef, 3)
#     fill!(α_init, zero(T))
#     central_normalized_moments!(uc, u[1:3])
#     normal_lagrange_params!(α, uc)
#     α_init[1:3] = α
# end

# """
#     lagrange_params!(α::AbstractArray{T}, u::AbstractArray{T}, p) where T

# Determine the Lagrange parameters for each moment vector stored as a column vector in the abstract array u.
# """
# function lagrange_params!(α::AbstractArray{T}, u::AbstractArray{T}, qr::Quadrature{T}) where T
#     αvec = Vector{T}(undef, size(α)[1])
#     success = true
#     for idx in 1:size(u)[2]
#         success &= lagrange_params!(αvec, view(u,:,idx), qr) 
#         α[:,idx] = αvec
#     end
#     return success
# end

# function lagrange_params!(α::AbstractVector{T}, u::AbstractVector{T}, qr::Quadrature{T}) where T
#     optim_fun = OptimizationFunction(dual_function; grad=dual_gradient!, hess=hessian!)
#     init_lagrange_params!(α, u)
#     params = (qr=qr, u=u)
#     prob = OptimizationProblem(optim_fun, α, params)
#     sol = solve(prob, Newton())
#     α[:] = sol
#     return sol.original.ls_success
# end

# function dual_function(α::AbstractVector{T}, p) where T
#     u0 = zeroth_moment(α, p.qr)
#     f = u0 - dot(α, p.u)
#     return f
# end

# function dual_function(α::AbstractArray{T}, p) where T
#     u0 = zeroth_moment(α, p.qr)
#     f = Vector{T}(undef, size(α)[2])
#     @tullio f[j] = u0[j] - α[i,j]*p.u[i,j]
#     return f
# end

# function dual_gradient!(g::AbstractArray{T}, α::AbstractArray{T}, p) where T
#     moments!(g, α, p.qr)
#     g[:] -= p.u[1:size(α)[1]]
#     return Nothing
# end

