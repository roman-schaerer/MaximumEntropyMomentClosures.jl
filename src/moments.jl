function moments!(u::AbstractVector{T}, α::AbstractVector{T}, 
    memc::MaximumEntropyMomentClosure{MomentSystemType, MomentEvaluationType}) where {T, MomentSystemType, MomentEvaluationType}
    moments!(u, α, memc.me)
end

function moments!(u::AbstractVector{T}, α::AbstractVector{T}, qr::Quadrature{T}) where T
    # u[i] = ith moment
    # α[i] = ith lagrange parameter
    # ϕ[i,j] = ith moment of basis function evaluated at node j

    # i = moment order
    # k = quadrature point
    # @infiltrate
    # for i in eachindex(u)
    #     u[i] = 0.0
    #     for j in 1:length(p.ϕ[i,:])
    #         u[i] += p.w[j] * p.ϕ[i,j] * exp(dot(α, p.ϕ[:,j]))
    #     end
    # end

    @tullio e[k] := α[l+0] * qr.ϕ[l+0,k]
    @tullio u[i+0] = qr.w[k] * qr.ϕ[i+0, k] * exp(e[k])
end

function moments!(u::AbstractArray{T}, α::AbstractArray{T}, qr::Quadrature{T}) where T
    # u[i,j] = jth moment in problem i
    # α[i,j] = jth lagrange parameter in problem i
    # ϕ[i,j] = ith moment of basis function evaluated at node 

    # i = moment order
    # j = cell index
    # k = quadrature point

    @tullio e[j,k] := α[l+0, j]*qr.ϕ[l+0, k]
    @tullio u[i+0,j] = qr.w[k] * qr.ϕ[i+0, k] * exp(e[j,k])
end

function hessian!(h::AbstractArray{T}, α::AbstractVector{T}, p) where T
    qr, _ = p
    hessian!(h, α, qr)    
end

function hessian!(h::AbstractArray{T}, α::AbstractVector{T}, qr::Quadrature{T}) where T
    m = length(α)
    fill!(h,0)
    @tullio e[k] := α[l+0]*qr.ϕ[l+0, k]

    ### # the following assumes a cartesian basis
    #@tullio h[i+0,1] = qr.w[k] * qr.ϕ[i+0, k] * exp(e[k])
    #@tullio h[$m, i+1] = qr.w[k] * qr.ϕ[i+$m, k] * exp(e[k]) 
    # for i in 1:m-1
    #     for j in 2:m-i+1
    #         h[i,j] = h[i+j-1,1] # since h[i,j] = h[i+k,j-k]
    #     end
    #     for j in m-i+2:m
    #         h[i,j] = h[m,j-m+i] # since h[i,j] = h[i+k,j-k]
    #     end
    # end
    ###

    # evaluate hessian elements in the uppwer triangular matrix (including the diagonal)
    for j in 1:m
        for i in 1:j
            @tullio h[$i, $j] = qr.w[k] * qr.ϕ[$i, k] * qr.ϕ[$j, k] * exp(e[k])
        end
    end 

    # populate strictly lower triangular matrix from the upper triangular matrix
    for j in 1:m
        for i in 1:j-1
            h[j, i] = h[i, j]
        end
    end 

end


function zeroth_moment(α::AbstractVector{T}, qr::Quadrature{T}) where T
    @tullio e[k] := α[l+0]*qr.ϕ[l+0,k]
    @tullio f := qr.w[k] * qr.ϕ[1,k] * exp(e[k])
    return f
end

function zeroth_moment(α::AbstractArray{T}, qr::Quadrature{T}) where T
    @tullio e[i,k] := α[l+0,i]*qr.ϕ[l+0,k]
    @tullio f[i] := qr.w[k] * qr.ϕ[1,k] * exp(e[i,k])
    return f
end

function zeroth_moment_cartesian_basis(α::AbstractVector{T}, qr::Quadrature{T}) where T
    @tullio e[k] := α[l+0]*qr.ϕ[l+0,k]
    @tullio f := qr.w[k] * exp(e[k])
    return f
end

function central_normalized_moments!(ucn::AbstractVector{T}, u::AbstractVector{T}) where T
    fill!(ucn, zero(T)) # initialize ucn as a zero vector

    ucn[1] = u[1] # 0th moment
    μ = u[2]/ucn[1] # mean (normalized first central moment)
    ucn[2] = μ
    # uc[3] = u[3]/uc[1] - uc[2]^2 # variance - second central moment
    # @assert uc[1] > 0 && uc[3] > 0
    m = length(u)
    for i in 2:m-1
        for j in 0:i
            ucn[i+1] += (-μ)^(i-j)*binomial(i,j)*u[j+1]
        end
        ucn[i+1] /= u[1]
    end
    #check that the central normalized moment of even order are positive
    @assert all(ucn[1:2:end] .> 0)
end

function raw_moments!(u::AbstractVector{T}, ucn::AbstractVector{T}) where T
    fill!(u, zero(T)) # initialize u as a zero vector
    u[1] = ucn[1] # zeroth moment
    μ = ucn[2] # mean = normalized first central moment
    m = length(u) # number of moments
    for i in 1:m-1 # loop over moment orders
        u[i+1] += (μ)^(i)
        for j in 2:i # loop over binomial expansion
            u[i+1] += (μ)^(i-j)*binomial(i,j)*ucn[j+1]
        end
        u[i+1] *= ucn[1]
    end
end
