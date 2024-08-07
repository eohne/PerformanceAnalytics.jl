# Struct


## Asset  
These structs contain the price or return values, the timestamp information, the scaling factor to annualize the returns, the frequency of returns, the id, and the exchange the asset is traded on.  

If AssetReturn(ticker::String,from::TimeType,to::TimeType,freq::String) is called data is automatically downloaded using YFinance.jl and the exchange, freq, and scale are automatically stored/calculated. The package calculates the scale assuming 252 trading days in a year. Furthermore, the exchange is important to accurately aggregate intra-day data (the trading hours in a day are different for exchanges across the world). 

````@docs
AssetPrice
AssetReturn
````

## Accessing Items

````@docs
values
timestamp
scale
id
````



