abstract type MomentSystem end

struct MomentSystemFullTensor{DomainType<:Domain, IntType<:Integer} <: MomentSystem
    domain::DomainType
    degree::IntType
    MomentSystemFullTensor(domain::DomainType, degree::IntType) where 
        {DomainType<:Domain, IntType<:Integer} = 
        new{DomainType, IntType}(domain, degree)
end
