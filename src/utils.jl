function double_factorial(n::Integer)
    p = 1
    for i in range(1+iseven(n),n,step=2)
        p *= i
    end
    return p
end