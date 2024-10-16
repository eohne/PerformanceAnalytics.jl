# PerformanceAnalytics.jl
[![Build Status](https://github.com/eohne/PerformanceAnalytics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eohne/PerformanceAnalytics.jl/actions/workflows/CI.yml?query=branch%3Amain) [![][docs-latest-img]][docs-latest-url]   

 
Given the lack of packages that deal with financial analysis in Julia I have decided to translate the functionality of the PerformanceAnalytics package in R to Julia.  


Since I took a long break from this and never published it a new package has emerged. This package is registered and implements many but not all the functions available here: [PortfolioAnalytics.jl](https://github.com/doganmehmet/PortfolioAnalytics.jl)

# List of all functions implemented:

## Asset Price and Return Functions

- ```julia
  AssetPrice(ticker::String, from::TimeType, to::TimeType, freq::String; exchange_local_time=true)
  ```
  - Downloads prices using YFinance.jl
- ```julia
  AssetReturn(x::AssetPrice)
  ```
- ```julia
  AssetReturn(ticker::String, from::TimeType, to::TimeType, freq::String)
  ```
  - Downloads prices using YFinance.jl and converts them to returns

### Return Calculation Functions
- ```julia
  log_diff(x)
  ```
- ```julia
  pct_change(x)
  ```
- ```julia
  simple_diff(x)
  ```

## Statistical Measures

### Basic Measures
- ```julia
  mean_arith(x)
  ```
- ```julia
  mean_geo(returns::Vector{Float64})
  ```
- ```julia
  mean_geo_log(returns::Vector{Float64})
  ```
- ```julia
  stdvp(x)
  ```
- ```julia
  stdvs(x)
  ```
- ```julia
  varp(x)
  ```
- ```julia
  vars(x)
  ```

### Advanced Measures
- ```julia
  skew(R, method="simple")
  ```
  - Methods: Population, Sample, Fisher
- ```julia
  kurt(R, method="simple")
  ```
  - Methods: population, excess, sample, sampleexcess, fisher
- ```julia
  covariance(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}; corrected::Bool=false)
  ```
  - Methods:
    - Normal (corrected=false, default)
    - Corrected (corrected=true)
- ```julia
  coskew(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number})
  ```
- ```julia
  cokurt(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number})
  ```
- ```julia
  hurstindex(R::AbstractArray{<:Number})
  ```

## Risk Measures

### Drawdown Measures
- ```julia
  drawdown(x::AbstractArray{<:Number})
  ```
- ```julia
  drawdownpeak(x::AbstractArray{<:Number})
  ```
- ```julia
  drawdown_table(x::AbstractArray{<:Number}, dates::AbstractArray{<:TimeType})
  ```
- ```julia
  avgdrawdown(x::AbstractArray{<:Number})
  ```
- ```julia
  maxdrawdown(x::AbstractArray{<:Number})
  ```

### Downside Risk Measures
- ```julia
  downsidedeviation(R::AbstractArray{<:Number}, MAR::Number = 0)
  ```
- ```julia
  downsidepotential(R::AbstractArray{<:Number}, MAR::Number=0)
  ```
- ```julia
  semideviation(R::AbstractArray{<:Number})
  ```
- ```julia
  semivariance(R::AbstractArray{<:Number})
  ```

### Value at Risk and Expected Shortfall
- ```julia
  valueatrisk(R::AbstractArray{<:Number}, p::Number, method::AbstractString="gaussian")
  ```
  - Methods: gaussian, historical, cornish_fisher
- ```julia
  expectedshortfall(R::AbstractArray{<:Number}, p::Number, method::AbstractString="gaussian")
  ```
  - Methods: gaussian, historical, cornish_fisher

### Other Risk Measures
- ```julia
  painindex(R::AbstractArray{<:Number})
  ```

## Performance Measures

### Risk-Adjusted Return Measures
- ```julia
  sharperatio(R::Number, σ::Number, Rf::Number=0.; scale=1)
  ```
- ```julia
  sharperatio(R::AbstractArray{<:Number}, Rf::Number=0.; scale)
  ```
- ```julia
  adjustedsharpe(R::AbstractArray{<:Number}, Rf::Number=0.; method::AbstractString="population")
  ```
- ```julia
  adjustedsharpe(R::Number, σ::Number, S::Number, K::Number, Rf::Number=0.)
  ```
- ```julia
  treynorratio(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0; modified::Bool=false, scale=1)
  ```
- ```julia
  informationratio(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, scale::Number=1)
  ```
- ```julia
  appraisalratio(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0, method::AbstractString = "appraisal", scale::Number=1)
  ```
  - Methods: 
    - "appraisal" (default)
    - "modified"
    - "alternative"

### Drawdown-based Performance Measures
- ```julia
  calmarratio(R::AbstractArray{<:Number}, scale::Float64=1.)
  ```
- ```julia
  burkeratio(R::AbstractArray{<:Number}; Rf::Float64 = 0., modified::Bool = false, scale::Number=1)
  ```
- ```julia
  burkeratio(Ra::AbstractArray{<:Number}, Drawdowns::AbstractArray{<:Number}; Rf::Number=0., scale=1)
  ```

### Relative Performance Measures
- ```julia
  activepremium(Ra::Number, Rb::Number; scalea=1, scaleb=1)
  ```
- ```julia
  jensensalpha(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; scale::Number=1)
  ```
- ```julia
  trackingerror(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, scale::Number=1)
  ```

### Other Performance Measures
- ```julia
  bernardoledoitratio(R::AbstractVector{<:Number})
  ```
- ```julia
  kappa(R::AbstractArray{<:Number}, MAR::Number=0., k::Number=2)
  ```
  - Notable cases:
    - k=1: returns the Sharpe-Omega Ratio
    - k=2: returns the Sortino Ratio

## Factor Analysis Functions

- ```julia
  factor_alpha(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...)
  ```
- ```julia
  factor_loadings(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; intercept::Bool=true)
  ```
- ```julia
  factor_regression(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; intercept::Bool=true)
  ```
- ```julia
  factor_resid(Y::AbstractArray{<:Number}, X::AbstractArray{<:Number}...; intercept::Bool=true)
  ```

## Risk Decomposition

- ```julia
  specificrisk(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0)
  ```
- ```julia
  systematicrisk(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0)
  ```
- ```julia
  totalrisk(Ra::AbstractArray{<:Number}, Rb::AbstractArray{<:Number}, Rf::Number=0)
  ```

## Utility Functions

- ```julia
  roll_apply(data::AbstractVector{<:Number}...; fun::Function, window::Int, retain_length::Bool=false, fill_with = NaN, always_return::Bool=false)
  ```
- ```julia
  to_annual_return(R::AssetReturn)
  ```
- ```julia
  to_monthly_return(R::AssetReturn)
  ```
- ```julia
  annualize(R::Number, scale::Float64=252)
  ```

A list of all functions and their definitions can be found in the documentation.


[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://eohne.github.io/PerformanceAnalytics.jl/dev/