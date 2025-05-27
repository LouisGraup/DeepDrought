# sensitivity analysis on LWFBrook90.jl

using CSV, DataFrames, Dates
using HypothesisTests
using StatsPlots, Plots; gr()

met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250523.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250523.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_20250523.csv", DataFrame);

# function for density plot of metrics
function density_plot(met)
    density(met.nse10, label="NSE10")
    density!(met.nse40, label="NSE40")
    density!(met.nse60, label="NSE60")
    density!(met.nse80, label="NSE80")
end

# function to filter metrics for behavioral runs
function behavioral_met(met)
    return met[met.nse10 .> 0 .&& 
               met.nse40 .> 0 .&&
               met.nse60 .> 0 .&&
               met.nse80 .> 0, :]
end

# function to separate parameters into behavioral and non-behavioral runs
function sep_params(par, met)
    met_good = behavioral_met(met);
    par_good = par[met_good.scen, :];
    
    met_bad = met[met.scen .∉ [met_good.scen], :];
    par_bad = par[met_bad.scen, :];

    return par_good, par_bad
end

# function to calculate and plot K-S statistic for parameters
function KS_plot(par, met)
    par_good, par_bad = sep_params(par, met);
    
    np = size(par, 2);
    ks_stat = DataFrame(par=names(par), D=zeros(np), pval=zeros(np));
    par_plots = [];
    
    # loop through each parameter
    for i in 1:np
        ks_test = ApproximateTwoSampleKSTest(par_good[!, i], par_bad[!, i]);
        ks_stat.D[i] = ks_test.δ; # D statistic
        ks_stat.pval[i] = pvalue(ks_test); # p-value

        p=ecdfplot(par_good[!, i], label="Behavioral")
        ecdfplot!(p, par_bad[!, i], label="Non-behavioral")
        title!(p, "K-S test for $(ks_stat.par[i]),\n D=$(round(ks_stat.D[i], sigdigits=3)), pval=$(round(ks_stat.pval[i], sigdigits=3))")
        push!(par_plots, p)
    end
    
    return ks_stat, par_plots
end

# separate parameters into behavioral and non-behavioral runs
function nse_quad(met)
    return sqrt.((met.nse10 .^ 2 + met.nse40 .^ 2 + met.nse60 .^ 2 + met.nse80 .^ 2) / 4);
end

density_plot(met_ctr)
density_plot(met_irr)

# filter metrics for behavioral runs
met_ctr_good = behavioral_met(met_ctr);
met_irr_good = behavioral_met(met_irr);

density_plot(met_ctr_good)
density_plot(met_irr_good)

# calculate quadratic mean of NSE values
met_ctr_good.nse_com = nse_quad(met_ctr_good);
met_irr_good.nse_com = nse_quad(met_irr_good);

describe(met_ctr_good)

# calculate K-S statistic to determine sensitive parameters
ks_stat_ctr, par_plots_ctr = KS_plot(par, met_ctr);
ks_stat_irr, par_plots_irr = KS_plot(par, met_irr);

plot(par_plots_ctr..., size=(1000,1000), layout=(3,5), legend=false, titlefontsize=8, guidefontsize=6)