################################
# Pct Change:
"""
    pct_change(x)

Returns the percent change of an array. (Calculates returns from a price series)
"""
function pct_change(x)
    return x[2:end]./x[1:(end-1)] .-1
end
# end pct_change

################################
# Log Change:
"""
    log_diff(x)

Returns the logarithmic difference. (Calculates log returns from a price series)
"""
function log_diff(x)
    return log.(x[2:end]) .- log.(x[1:(end-1)])
end
# end log diff

################################
# Simple Diff:
"""
    simple_diff(x)

Returns the difference. (Calculates the absolut price change)
"""
function simple_diff(x)
    return x[2:end] .- x[1:(end-1)]
end
# end simple diff