# PerformanceAnalytics.jl
[![Build Status](https://github.com/eohne/PerformanceAnalytics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eohne/PerformanceAnalytics.jl/actions/workflows/CI.yml?query=branch%3Amain) [![][docs-latest-img]][docs-latest-url]   

 
Given the lack of packages that deal with financial analysis in Julia I have decided to translate the functionality of the PerformanceAnalytics package in R to Julia.  


Since I took a long break from this and never published it a new package has emerged. This package is registered and implements many but not all the functions available here: [PortfolioAnalytics.jl](https://github.com/doganmehmet/PortfolioAnalytics.jl)

# List of all functions implemented:

## Asset Price and Return Functions

- `AssetPrice(ticker::String, from::TimeType, to::TimeType, freq::String; exchange_local_time=true)`
  - Downloads prices using YFinance.jl
- `AssetReturn(x::AssetPrice)`
- `AssetReturn(ticker::String, from::TimeType, to::TimeType, freq::String)`
  - Downloads prices using YFinance.jl and converts them to returns

### Return Calculation Functions
- `log_diff(x)`
- `pct_change(x)`
- `simple_diff(x)`

## Statistical Measures

### Basic Measures
- `mean_arith(x)`
- `mean_geo(returns::Vector{Float64})`
- `mean_geo_log(returns::Vector{Float64})`
- `stdvp(x)`
- `stdvs(x)`
- `varp(x)`
- `vars(x)`

### Advanced Measures
- `skew(R, method="simple")`
  - Methods: Population, Sample, Fisher
- `kurt(R, method="simple")`
  - Methods: population, excess, sample, sampleexcess, fisher
- `covariance(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}; corrected::Bool=false)`
  - Methods:
    - Normal (corrected=false, default)
    - Corrected (corrected=true)
- `coskew(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number})`
- `cokurt(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number})`
- `hurstindex(R::AbstractArray{<:Number})`

## Risk Measures

### Drawdown Measures
- `drawdown(x::AbstractArray{<:Number})`
- `drawdownpeak(x::AbstractArray{<:Number})`
- `drawdown_table(x::AbstractArray{<:Number}, dates::AbstractArray{<:TimeType})`
- `avgdrawdown(x::AbstractArray{<:Number})`
- `maxdrawdown(x::AbstractArray{<:Number})`

### Downside Risk Measures
- `downsidedeviation(R::AbstractArray{<:Number}, MAR::Number = 0)`
- `downsidepotential(R::AbstractArray{<:Number}, MAR::Number=0)`
- `semideviation(R::AbstractArray{<:Number})`
- `semivariance(R::AbstractArray{<:Number})`

### Value at Risk and Expected Shortfall
- `valueatrisk(R::AbstractArray{<:Number}, p::Number, method::AbstractString="gaussian")`
  - Methods: gaussian, historical, cornish_fisher
- `expectedshortfall(R::AbstractArray{<:Number}, p::Number, method::AbstractString="gaussian")`
  - Methods: gaussian, historical, cornish_fisher

### Other Risk Measures
- `painindex(R::AbstractArray{<:Number})`

## Performance Measures

### Risk-Adjusted Return Measures
- `sharperatio(R::Number, σ::Number, Rf::Number=0.; scale=1)`
- `sharperatio(R::AbstractArray{<:Number}, Rf::Number=0.; scale)`
- `adjustedsharpe(R::AbstractArray{<:Number}, Rf::Number=0.; method::AbstractString="population")`
- `adjustedsharpe(R::Number, σ::Number, S::Number, K::Number, Rf::Number=0.)`
- `treynorratio(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0; modified::Bool=false, scale=1)`
- `informationratio(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, scale::Number=1)`
- `appraisalratio(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0, method::AbstractString = "appraisal", scale::Number=1)`
  - Methods: 
    - "appraisal" (default)
    - "modified"
    - "alternative"

### Drawdown-based Performance Measures
- `calmarratio(R::AbstractArray{<:Number}, scale::Float64=1.)`
- `burkeratio(R::AbstractArray{<:Number}; Rf::Float64 = 0., modified::Bool = false, scale::Number=1)`
- `burkeratio(Ra::AbstractArray{<:Number}, Drawdowns::AbstractArray{<:Number}; Rf::Number=0., scale=1)`

### Relative Performance Measures
- `activepremium(Ra::Number, Rb::Number; scalea=1, scaleb=1)`
- `jensensalpha(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; scale::Number=1)`
- `trackingerror(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, scale::Number=1)`

### Other Performance Measures
- `bernardoledoitratio(R::AbstractVector{<:Number})`
- `kappa(R::AbstractArray{<:Number}, MAR::Number=0., k::Number=2)`
  - Notable cases:
    - k=1: returns the Sharpe-Omega Ratio
    - k=2: returns the Sortino Ratio

## Factor Analysis Functions

- `factor_alpha(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...)`
- `factor_loadings(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; intercept::Bool=true)`
- `factor_regression(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; intercept::Bool=true)`
- `factor_resid(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; intercept::Bool=true)`

## Risk Decomposition

- `specificrisk(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0)`
- `systematicrisk(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0)`
- `totalrisk(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0)`

## Utility Functions

- `roll_apply(data::AbstractVector{<:Number}...; fun::Function, window::Int, retain_length::Bool=false, fill_with = NaN, always_return::Bool=false)`
- `to_annual_return(R::AssetReturn)`
- `to_monthly_return(R::AssetReturn)`
- `annualize(R::Number, scale::Float64=252)`


A list of all functions and their definitions can be found in the documentation.


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://eohne.github.io/PerformanceAnalytics.jl/dev/