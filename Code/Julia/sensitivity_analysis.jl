# sensitivity analysis on LWFBrook90.jl

using CSV, DataFrames, Dates
using HypothesisTests
using StatsPlots, Plots; gr()

met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250526.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250526.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_20250526.csv", DataFrame);

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

# function to plot parameter relationships with combined metric
function par_plot(par, met; met_y="rmse_com", behave=true)
    
    # add combined metric if desired
    if met_y == "rmse_com"
        met.rmse_com = rmse_quad(met);
    end

    if behave
        # only plot behavioral parameter sets
        par, _ = sep_params(par, met);
        met = behavioral_met(met);
    end

    np = size(par, 2);
    par_plots = [];
    
    # loop through each parameter
    for i in 1:np
        p = scatter(par[!, i], met[!, met_y], xlabel=names(par)[i], ylabel=met_y, 
        mc="black", msc="black", ms = 1.5, legend=false)
        title!(p, "$(names(par)[i]) vs $(met_y)")
        push!(par_plots, p)
    end
    
    return par_plots
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

# function to calculate combined NSE metric
function nse_quad(met)
    return sqrt.((met.nse10 .^ 2 + met.nse40 .^ 2 + met.nse60 .^ 2 + met.nse80 .^ 2) / 4);
end

# function to calculate combined RMSE metric
function rmse_quad(met)
    return sqrt.((met.rmse10 .^ 2 + met.rmse40 .^ 2 + met.rmse60 .^ 2 + met.rmse80 .^ 2) / 4);
end

density_plot(met_ctr)
density_plot(met_irr)

# filter metrics for behavioral runs
met_ctr_good = behavioral_met(met_ctr);
met_irr_good = behavioral_met(met_irr);

density_plot(met_ctr_good)
density_plot(met_irr_good)

describe(met_ctr_good)

# compare metrics
@df met_ctr_good scatter(:nse10, :nse40, xlabel="NSE10", ylabel="NSE40", legend=false)
@df met_ctr_good scatter(:nse40, :nse60, xlabel="NSE40", ylabel="NSE60", legend=false)
@df met_ctr_good scatter(:nse60, :nse80, xlabel="NSE60", ylabel="NSE80", legend=false)
@df met_ctr_good scatter(:nse10, :nse80, xlabel="NSE10", ylabel="NSE80", legend=false)

# parameter relationships
par_plots_ctr = par_plot(par, met_ctr);
par_plots_irr = par_plot(par, met_irr);

plot(par_plots_ctr..., size=(1000,1000), layout=(3,5), legend=false, titlefontsize=8, guidefontsize=6)

# calculate K-S statistic to determine sensitive parameters
ks_stat_ctr, ks_plots_ctr = KS_plot(par, met_ctr);
ks_stat_irr, ks_plots_irr = KS_plot(par, met_irr);

plot(ks_plots_ctr..., size=(1000,1000), layout=(3,5), legend=false, titlefontsize=8, guidefontsize=6)