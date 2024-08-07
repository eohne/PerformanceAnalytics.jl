################################
# Stats:

# mean:
"""
    mean_arith(x)

Calculates the arithmetic mean of a series.
"""
function mean_arith(x)
    return sum(x)/length(x)
end
function mean_arith(x::Asset)
    return sum(x.values)/length(x.values)
end

"""
    mean_geo(returns::Vector{Float64})

Calculates the geometric mean of a series.
"""
function mean_geo(returns::AbstractArray{Float64})
    n = length(returns)
    return prod(1 .+ returns)^(1/n) - 1
end
function mean_geo(returns::Asset)
    n = length(returns.values)
    return prod(1 .+ returns.values)^(1/n) - 1
end

"""
    mean_geo_log(returns::Vector{Float64})

Calculates the geometric mean of a series.
"""
function mean_geo_log(returns::Vector{Float64})
    return exp(mean(log.(1 .+ returns))) - 1
end
function mean_geo_log(returns::Asset)
    return exp(mean(log.(1 .+ returns.values))) - 1
end
#

# Standard Deviation:
"""
    varp(x)

Calculates the population variance of a series.
"""
function varp(x)
    return sum( (x .- mean_arith(x)).^2 )/length(x)
end
function varp(x::Asset)
    return sum( (x.values .- mean_arith(x.values)).^2 )/length(x.values)
end
"""
    stdvp(x)

Calculates the sample variance of a series.
"""
function stdvp(x)
    return sqrt(varp(x))
end
function stdvp(x::Asset)
    return sqrt(varp(x.values))
end
"""
    vars(x)

Calculates the population standard deviation of a series.
"""
function vars(x)
    return sum( (x .- mean_arith(x)).^2 )/(length(x)-1)
end
function vars(x::Asset)
    return sum( (x.values .- mean_arith(x.values)).^2 )/(length(x.values)-1)
end

"""
    stdvs(x)

Calculates the sample standard deviation of a series.
"""
function stdvs(x)
    return sqrt(vars(x))
end
function stdvs(x::Asset)
    return sqrt(vars(x.values))
end
#

# Skewness:
"""
    skew(R, method="simple")

Calculates the skewness of a series.

# Methods: 
   * Population
   * Sample
   * Fisher
"""
function skew(R, method="simple")
    n = length(R)
    if isequal(method,"sample")
        return ( n/((n-1)*(n-2)) )* sum( ((R.-mean_arith(R))/stdvs(R)).^3 )
    elseif isequal(method,"fisher")
        return (sqrt(n*(n-1))/(n-2))*( sum( ((R.-mean_arith(R)).^3)./n ))/( sum( ((R.-mean_arith(R)).^2)./n)^(3/2))
    elseif isequal(method,"simple")
        return (1/n) * sum( ((R.-mean_arith(R))/stdvp(R)).^3 )
    else
        throw(ArgumentError("$method does not exist please choose one from: simple, sample, fisher"))
    end
end
function skew(R::Asset, method="simple")
    n = length(R.values)
    if isequal(method,"sample")
        return ( n/((n-1)*(n-2)) )* sum( ((R.values .- mean_arith(R))/stdvs(R)).^3 )
    elseif isequal(method,"fisher")
        return (sqrt(n*(n-1))/(n-2))*( sum( ((R.values .-mean_arith(R)).^3)./n ))/( sum( ((R.values .- mean_arith(R)).^2)./n)^(3/2))
    elseif isequal(method,"simple")
        return (1/n) * sum( ((R.values .-mean_arith(R))/stdvp(R)).^3 )
    else
        throw(ArgumentError("$method does not exist please choose one from: simple, sample, fisher"))
    end
end
#

# Kurtosis:
"""
    kurt(R, method="simple")

Calculates the kurtosis of a series.

# Methods: 
   * population
   * excess (population - 3)
   * sample
   * sampleexcess (sample - 3)
   * fisher
"""
function kurt(R, method="simple")
    n = length(R)
    if isequal(method,"excess")
        return (1/n * sum( ((R .- mean_arith(R))./stdvp(R)).^4 )) - 3
    elseif isequal(method,"sample")
        return (n*(n+1))/((n-1)*(n-2)*(n-3)) * (sum( ((R .- mean_arith(R))./stdvs(R)).^4 ))
    elseif isequal(method,"sampleexcess")
        return (n*(n+1))/((n-1)*(n-2)*(n-3)) * (sum( ((R .- mean_arith(R))./stdvs(R)).^4 )) - ((3*((n-1)^2))/((n-2)*(n-3)))
    elseif isequal(method,"fisher")
        return ( ((n+1)*(n-1))/((n-2)*(n-3)) ) *(  (sum(((R.-mean_arith(R)).^4)./n)/(sum(((R.-mean_arith(R)).^2)./n)^2)) - ((3*(n-1))/(n+1))) 
    elseif isequal(method,"simple")
        return (1/n * sum( ((R .- mean_arith(R))./stdvp(R)).^4 ))
    else
        throw(ArgumentError("$method does not exist please choose one from: simple, excess, sample, sampleexcess, fisher"))
    end
end
function kurt(R::Asset, method="simple")
    n = length(R.values)
    if isequal(method,"excess")
        return (1/n * sum( ((R.values .- mean_arith(R))./stdvp(R)).^4 )) - 3
    elseif isequal(method,"sample")
        return (n*(n+1))/((n-1)*(n-2)*(n-3)) * (sum( ((R.values .- mean_arith(R))./stdvs(R)).^4 ))
    elseif isequal(method,"sampleexcess")
        return (n*(n+1))/((n-1)*(n-2)*(n-3)) * (sum( ((R.values .- mean_arith(R))./stdvs(R)).^4 )) - ((3*((n-1)^2))/((n-2)*(n-3)))
    elseif isequal(method,"fisher")
        return ( ((n+1)*(n-1))/((n-2)*(n-3)) ) *(  (sum(((R.values .-mean_arith(R)).^4)./n)/(sum(((R.values .-mean_arith(R)).^2)./n)^2)) - ((3*(n-1))/(n+1))) 
    elseif isequal(method,"simple")
        return (1/n * sum( ((R.values .- mean_arith(R))./stdvp(R)).^4 ))
    else
        throw(ArgumentError("$method does not exist please choose one from: simple, excess, sample, sampleexcess, fisher"))
    end
end
#
# End Stats



################################
# Regression Based Metrics:
"""
    factor_regression(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...;intercept::Bool=true)

Returns alpha and beta estimates of a regression.
"""
function factor_regression(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...;intercept::Bool=true)
    if intercept
        Xs = hcat(ones(length(Y)),X...)
    else
        Xs = hcat(X...)
    end
    β = inv(Xs'*Xs)*Xs'*Y
    return β
end

function factor_regression(Y::Asset,X::Asset...;intercept::Bool=true)
    if intercept
        Xs = hcat(ones(length(Y.values)),values.(X)...)
    else
        Xs = hcat(values.(X)...)
    end
    β = inv(Xs'*Xs)*Xs'*Y.values
    return β
end


"""
    factor_alpha(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...)

Returns the α (intercept) of a (factor) regression.
"""
function factor_alpha(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...)
    Xs = hcat(ones(length(Y)),X...)
    α = (inv(Xs'*Xs)*Xs'*Y)[1]
    return α 
end
function factor_alpha(Y::Asset,X::Asset...)
    Xs = hcat(ones(length(Y.values)),values.(X)...)
    α = (inv(Xs'*Xs)*Xs'*Y.values)[1]
    return α 
end

# end regression

"""
    factor_loadings(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...;intercept::Bool=true)

Returns the factor loadings (beta estimates) of a (factor) regression.
"""
function factor_loadings(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...;intercept::Bool=true)
    if intercept
        Xs = hcat(ones(length(Y)),X...)
    else
        Xs = hcat(X...)
    end
    β = (inv(Xs'*Xs)*Xs'*Y)[2:end]
    return β
end
function factor_loadings(Y::Asset,X::Asset...;intercept::Bool=true)
    if intercept
        Xs = hcat(ones(length(Y.values)),values.(X)...)
    else
        Xs = hcat(values.(X)...)
    end
    β = (inv(Xs'*Xs)*Xs'*Y.values)[2:end]
    return β
end
# end regression


"""
    factor_resid(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...;intercept::Bool=true)

Returns the residuals of a (factor) regression.
"""
function factor_resid(Y::AbstractArray{<:Number},X::AbstractArray{<:Number}...;intercept::Bool=true)
    if intercept
        Xs = hcat(ones(length(Y)),X...)
    else
        Xs = hcat(X...)
    end
    β = (inv(Xs'*Xs)*Xs'*Y)
    return Y - (β'*Xs')'
end
function factor_resid(Y::Asset,X::Asset...;intercept::Bool=true)
    if intercept
        Xs = hcat(ones(length(Y.values)),values.(X)...)
    else
        Xs = hcat(values.(X)...)
    end
    β = (inv(Xs'*Xs)*Xs'*Y.values)
    return Y.values - (β'*Xs')'
end
# End regression based metrics