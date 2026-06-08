# sensitivity analysis on LWFBrook90.jl

using CSV, DataFrames, Dates
using HypothesisTests
using StatsPlots, Plots; gr()

# function to filter out scenarios which produced error
function filter_error(met)
    println("Filtering out $(sum(met.swp_nse20 .== 0 .&& met.swp_nse80 .== 0 .&& met.swp_nse110 .== 0)) scenarios out of total $(size(met, 1)) which failed to run.")
    return met[met.swp_nse20 .!= 0 .&& met.swp_nse80 .!= 0 .&& met.swp_nse110 .!= 0, :]
end

# function for density plot of metrics
function density_plot(met)
    density(met.swp_nse20, label="NSE20")
    density!(met.swp_nse80, label="NSE80")
    density!(met.swp_nse110, label="NSE110")
    density!(met.swp_nse160, label="NSE160")
end

# function to filter metrics for behavioral runs
function behavioral_met(met)

    return(met[met.swp_nse20 .> 0.62 .&&
               met.swp_nse80 .> 0.7 .&&
               met.swp_nse110 .> 0.8 .&&
               met.swp_nse160 .> 0.7 .&&
               met.trans_cor .> 0.6, :])

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
function par_plot(par, met; met_y="swp_nse_com", behave=true)

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

# function to plot metric comparisons
function met_plot(met, x, y)
    
    if isa(x, Array)
        # create subplots
        n = length(x);
        met_plots = [];
        for i in 1:n
            p = scatter(met[!, x[i]], met[!, y[i]], xlabel=x[i], ylabel=y[i], legend=false)
            push!(met_plots, p)
        end
        return plot(met_plots..., size=(1000, 1000), layout=(ceil(Int, sqrt(n)), ceil(Int, sqrt(n))))
    else
        return scatter(met[!, x], met[!, y], xlabel=x, ylabel=y, legend=false)
    end

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

# function to add combined metrics to DataFrame
function met_comb!(met)
    met.swp_nse_com = sqrt.((met.swp_nse20 .^ 2 + met.swp_nse80 .^ 2 + met.swp_nse110 .^ 2 + met.swp_nse160 .^ 2) / 4);
    return nothing
end

# function to find best scenario
function met_best_scen(met, metric=:swp_nse_com)
    # find the index of the maximum value
    max_idx = argmax(met[!, metric]);
    best_scen = met.scen[max_idx];
    
    return best_scen, met[max_idx, :]
end

# calibration results
met = CSV.read("LWFBcal_output/metrics_vetroz_20260604.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_vetroz_20260604.csv", DataFrame);

# filter out scenarios which produced an error
met = filter_error(met);

# exploratory
density_plot(met)

# trans density plots
density(met.trans_cor, label="Corr Coef")

density(met.max_trans)

# add combined metrics
met_comb!(met);

# filter metrics for behavioral runs
met_good = behavioral_met(met);
println("$(size(met_good, 1)) behavioral parameter sets out of total $(size(met, 1)) in control scenario.")
describe(met_good)

density_plot(met_good)


# compare metrics across depths
met_plot(met_good, [:swp_nse20, :swp_nse80, :swp_nse110, :swp_nse160], [:swp_nse80, :swp_nse110, :swp_nse160, :swp_nse20])

met_plot(met_good, [:swp_nse20, :swp_nse80, :swp_nse110, :swp_nse160], [:trans_cor, :trans_cor, :trans_cor, :trans_cor])

# best control scenario
scen_max, met_max = met_best_scen(met_good);
# parameter values for the best performing scenario
par_best = par[scen_max, :];
par_best

# parameter relationships
par_plots = par_plot(par, met);

plot(par_plots..., size=(1000,1000), layout=(5,6), legend=false, titlefontsize=8, guidefontsize=6)

# calculate K-S statistic to determine sensitive parameters
ks_stat, ks_plots = KS_plot(par, met);

plot(ks_plots..., size=(1200,1200), layout=(5,8), legend=false, titlefontsize=8, guidefontsize=6)
# behavioral is blue, non-behavioral is red