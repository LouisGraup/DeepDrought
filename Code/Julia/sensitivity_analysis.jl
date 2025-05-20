# sensitivity analysis on LWFBrook90.jl

using CSV, DataFrames, Dates
using HypothesisTests
using StatsPlots, Plots; gr()

met = CSV.read("LWFBcal_output/metrics_20250519.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_20250519.csv", DataFrame);

# density plots of metrics
density(met.nse10, label="NSE10")
density!(met.nse40, label="NSE40")
density!(met.nse60, label="NSE60")
density!(met.nse80, label="NSE80")

# behavioral runs
met_good = met[met.nse10 .> 0 .&& 
              met.nse40 .> 0 .&&
              met.nse60 .> 0 .&&
              met.nse80 .> 0, :];

density(met_good.nse10, label="NSE10")
density!(met_good.nse40, label="NSE40")
density!(met_good.nse60, label="NSE60")
density!(met_good.nse80, label="NSE80")

# calculate quadratic mean of NSE values
met_good.nse_com = sqrt.((met_good.nse10 .^ 2 + met_good.nse40 .^ 2 + met_good.nse60 .^ 2 + met_good.nse80 .^ 2) / 4);

describe(met_good)

# compare behavioral and non-behavioral runs
par_good = par[met_good.scen, :];

met_bad = met[met.scen .∉ [met_good.scen], :];
par_bad = par[met_bad.scen, :];

# calculate K-S statistic to determine sensitive parameters

np = size(par, 2);
par_plots = [];
ks_stat = DataFrame(par=names(par), D=zeros(np), pval=zeros(np));
for i in 1:np
    ks_test = ApproximateTwoSampleKSTest(par_good[!, i], par_bad[!, i]);
    ks_stat.D[i] = ks_test.δ; # D statistic
    ks_stat.pval[i] = pvalue(ks_test); # p-value

    p=ecdfplot(par_good[!, i], label="Behavioral")
    ecdfplot!(p, par_bad[!, i], label="Non-behavioral")
    title!(p, "K-S test for $(ks_stat.par[i]),\n D=$(round(ks_stat.D[i], sigdigits=3)), pval=$(round(ks_stat.pval[i], sigdigits=3))")
    display(p)
    push!(par_plots, p)
end

plot(par_plots..., size=(1000,1000), layout=(4,5), legend=false, titlefontsize=8, guidefontsize=6)