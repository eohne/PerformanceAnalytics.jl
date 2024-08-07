push!(LOAD_PATH,"../src/")
using PerformanceAnalytics
using Documenter
makedocs(
         sitename = "PerformanceAnalytics.jl",
         modules  = [PerformanceAnalytics],
         format = Documenter.HTML(analytics = "G-TTP1GKJ1LF",
         canonical = "https://eohne.github.io/PerformanceAnalytics.jl/stable/"),
         pages=[
                "Home" => "index.md",
                "Function Documentation" =>[
                    "Struct" => "Struct.md",
                    "Statistics" =>"Statistics.md",
                    "Performance Measures" =>"Performance.md",
                    "Risk Measures" =>"RiskMeasures.md",
                    "All Functions" =>"AllFunctions.md",
                ],
                "Example Usage" => [
                    "Placeholder" => "Placeholder.md"
                ],
                "Version Change Log" => "VersionChanges.md"
               ])
deploydocs(;
    repo="github.com/eohne/PerformanceAnalytics.jl",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "dev" => "dev"]
)
