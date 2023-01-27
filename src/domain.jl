abstract type Domain end

@with_kw struct Interval{T} <: Domain
    lower::T # lower bound
    upper::T # upper bound
end