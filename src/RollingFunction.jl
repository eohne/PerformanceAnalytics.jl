# Apply function rolling window:
"""
    roll_apply(data::AbstractVector{<:Number}...;fun::Function, window::Int,retain_length::Bool=false,fill_with = NaN,always_return::Bool=false)

Applies a `fun::Function` over a rolling window.

# Arguments:
   * `data::AbstractArray{<:Number}`: one or multiple arrays.
   * `fun::Function`: the function to be applied. Has to return a single value (not a Vector of values).
   * `window::Int`: the window length.
   * `retain_length::Bool`: whether to original length of the data provided should be retained. If true fills the observations that could not be calculated with `fill_with`.
   * `fill_with`: The value that should be returned for values that could not be calculated.
   * `always_return::Bool`: returns a vector with `fill_with` equal to the initial data length if the window is longer than the length of the input data.
"""
function roll_apply(data::AbstractVector{<:Number}...;fun::Function, window::Int,retain_length::Bool=false,fill_with = NaN,always_return::Bool=false)
    if always_return && (window>length(data[1]))
        return [fill_with for i in eachindex(data[1])]
    end
    @assert window <= length(data[1]) "Window length must be shorter or equal to the length of the input data"
    nvals  = length(data[1]) - window +1
    offset = window - 1
    res = zeros(nvals)
    if length(data)>1
        @inbounds for i in eachindex(res)
            res[i] = fun( view.(data, (i:i+offset,))...)[1]
        end
    else
        @inbounds for i in eachindex(res)
            res[i] = fun( view.(data, (i:i+offset,))...)[1]
        end
    end
    if retain_length
        res = vcat( [fill_with for i in 1: (length(data[1]) - length(res))],res )
    end
    return res
end

function roll_apply(TA::AssetReturn;fun::Function, window::Int,retain_length::Bool=false,fill_with = NaN,always_return::Bool=false)
    d = timestamp(TA)
    data = [values(TA[i]) for i in colnames(TA)]
    if always_return && (window>length(data[1]))
        return [fill_with for i in eachindex(data[1])]
    end
    @assert window <= length(data[1]) "Window length must be shorter or equal to the length of the input data"
    nvals  = length(data[1]) - window +1
    offset = window - 1
    res = zeros(nvals)
    i=1
    if length(data)>1
        @inbounds for i in eachindex(res)
            res[i] = fun( view.(data, (i:i+offset,))...)[1]
        end
    else
        @inbounds for i in eachindex(res)
            res[i] = fun( view.(data, (i:i+offset,))...)[1]
        end
    end
    if retain_length
        res = vcat( [fill_with for i in 1: (length(data[1]) - length(res))],res )
    else
        d=d[window:end]
    end
    return AssetReturn(res,d,TA.freq,TA.scale,TA.id,TA.exchange)
end


# end rolling window