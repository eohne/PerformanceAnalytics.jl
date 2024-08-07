module PerformanceAnalytics
    # Dependencies:
    using  YFinance, Dates
    import StatsBase: quantile
    import Distributions: Normal,pdf

    # From Structs.jl
    export Asset, AssetReturn, AssetPrice, show, summary,values, timestamp, scale, id

    # From Stats.jl
    export mean_arith,mean_geo,varp,vars,stdvp,stdvs,skew,kurt
    export factor_regression,factor_alpha,factor_loadings,factor_resid

    #From Returns.jl
    export pct_change,log_diff,simple_diff

    #From RollingFunction.jl
    export roll_apply

    #From ReturnAggregation.jl
    export to_monthly_return, to_annual_return

    #From PortfolioStats.jl
    export drawdown,drawdownpeak,maxdrawdown,avgdrawdown,drawdown_table,annualize
    export activepremium,sharperatio,adjustedsharpe,bernardoledoitratio
    export burkeratio,calmarratio,downsidedeviation,downsidepotential,semideviation
    export semivariance,valueatrisk,expectedshortfall,hurstindex,kappa,painindex,covariance
    export coskew,cokurt,specificrisk,systematicrisk,totalrisk,trackingerror,jensensalpha
    export informationratio,treynorratio,appraisalratio

    include("Structs.jl")
    include("Stats.jl");
    include("Returns.jl");
    include("RollingFunction.jl");
    include("ReturnAggregation.jl");
    include("PortfolioStats.jl");
end

