# sensitivity analysis on LWFBrook90.jl

using CSV, DataFrames, Dates
using HypothesisTests
using StatsPlots, Plots; gr()

# function to filter out scenarios which produced error
function filter_error(met)
    println("Filtering out $(sum(met.swc_rmse10 .== 0)) scenarios out of total $(size(met, 1)) which failed to run.")
    return met[met.swc_rmse10 .!= 0, :]
end

# function for density plot of metrics
function density_plot(met)
    density(met.swc_nse10, label="NSE10")
    density!(met.swc_nse40, label="NSE40")
    density!(met.swc_nse60, label="NSE60")
    density!(met.swc_nse80, label="NSE80")
end

# function to filter metrics for behavioral runs
function behavioral_met(met)
    # control metrics
    #= return met[met.swc_nse10 .> 0.75 .&& 
               met.swc_nse40 .> 0.75 .&&
               met.swc_nse60 .> 0.85 .&&
               #met.swc_nse80 .> 0.5, :]
               met.swc_nse80 .> 0.9 .&&
               met.swp_nse10 .> 0.85 .&&
               met.swp_nse80 .> 0.8, :] =#

    # irrigation metrics
    return met[met.swc_nse10 .> 0.55 .&& 
               met.swc_nse40 .> 0.4 .&&
               met.swc_nse60 .> 0.4 .&&
               #met.swc_nse80 .> 0.5, :]
               met.swc_nse80 .> 0.4 .&&
               met.swp_nse10 .> 0.35 .&&
               met.swp_nse80 .> 0.35, :]
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

# function to compare scenarios across parameters
function scen_plot(par_ctr, par_irr, met_ctr, met_irr)
    # retrieve behavioral parameters for both scenarios
    par_ctr_good, = sep_params(par_ctr, met_ctr);
    par_irr_good, = sep_params(par_irr, met_irr);
    
    np = size(par_ctr, 2);
    ks_stat = DataFrame(par=names(par_ctr), D=zeros(np), pval=zeros(np));
    par_plots = [];
    
    # loop through each parameter
    for i in 1:np
        ks_test = ApproximateTwoSampleKSTest(par_ctr_good[!, i], par_irr_good[!, i]);
        ks_stat.D[i] = ks_test.δ; # D statistic
        ks_stat.pval[i] = pvalue(ks_test); # p-value

        p=ecdfplot(par_ctr_good[!, i], label="C")
        ecdfplot!(p, par_irr_good[!, i], label="I")
        title!(p, "K-S test for $(ks_stat.par[i]),\n D=$(round(ks_stat.D[i], sigdigits=3)), pval=$(round(ks_stat.pval[i], sigdigits=3))")
        push!(par_plots, p)
    end
    
    return ks_stat, par_plots
end

# function to add combined metrics to DataFrame
function met_comb!(met)
    met.swc_nse_com = sqrt.((met.swc_nse10 .^ 2 + met.swc_nse40 .^ 2 + met.swc_nse60 .^ 2 + met.swc_nse80 .^ 2) / 4);
    met.swp_nse_com = sqrt.((met.swp_nse10 .^ 2 + met.swp_nse80 .^ 2) / 2);
    met.rmse_com = sqrt.((met.swc_rmse10 .^ 2 + met.swc_rmse40 .^ 2 + met.swc_rmse60 .^ 2 + met.swc_rmse80 .^ 2) / 4);
    met.met_com = sqrt.((met.swc_nse_com .^ 2 + met.swp_nse_com .^ 2) / 2);
    return nothing
end

# function to find best scenario
function met_best_scen(met, metric=:swc_nse_com)
    # find the index of the maximum value
    max_idx = argmax(met[!, metric]);
    best_scen = met.scen[max_idx];
    
    return best_scen, met[max_idx, :]
end

# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250814.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250814.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_ctr_20250814.csv", DataFrame);
par_irr = CSV.read("LWFBcal_output/param_irr_20250814.csv", DataFrame);

# filter out scenarios which produced an error
met_ctr = filter_error(met_ctr);
met_irr = filter_error(met_irr);

density_plot(met_ctr)
density_plot(met_irr)

# add combined metrics
met_comb!(met_ctr);
met_comb!(met_irr);

# filter metrics for behavioral runs
met_ctr_good = behavioral_met(met_ctr);
met_irr_good = behavioral_met(met_irr);

println("$(size(met_ctr_good, 1)) behavioral parameter sets out of total $(size(met_ctr, 1)) in control scenario.")
println("$(size(met_irr_good, 1)) behavioral parameter sets out of total $(size(met_irr, 1)) in irrigation scenario.")

density_plot(met_ctr_good)
density_plot(met_irr_good)

describe(met_ctr_good)

# compare metrics across depths
met_plot(met_ctr_good, [:swc_nse10, :swc_nse40, :swc_nse60, :swc_nse10], [:swc_nse40, :swc_nse60, :swc_nse80, :swc_nse80])
met_plot(met_irr_good, [:swc_nse10, :swc_nse40, :swc_nse60, :swc_nse10], [:swc_nse40, :swc_nse60, :swc_nse80, :swc_nse80])

met_plot(met_ctr_good, :swp_nse10, :swp_nse80)
met_plot(met_irr_good, :swp_nse10, :swp_nse80)

met_plot(met_ctr_good, [:swc_nse10, :swc_nse80], [:swp_nse10, :swp_nse80])
met_plot(met_irr_good, [:swc_nse10, :swc_nse80], [:swp_nse10, :swp_nse80])

# compare metrics across type
met_plot(met_ctr_good, :swc_nse_com, :swp_nse_com)
met_plot(met_irr_good, :swc_nse_com, :swp_nse_com)

# compare metrics across scenarios
scatter(met_ctr.swc_nse10, met_irr.swc_nse10, xlabel="NSE10 Control", ylabel="NSE10 Irrigation", legend=false)

# compare behavioral parameter sets across both scenarios
ks_stat, ks_plots = scen_plot(par_ctr, par_irr, met_ctr, met_irr);
plot(ks_plots..., size=(1000,1000), layout=(4,5), legend=:bottomright, titlefontsize=8, guidefontsize=6)

# find common parameters
common_params = intersect(met_ctr_good.scen, met_irr_good.scen);

met_ctr_best = met_ctr_good[met_ctr_good.scen .∈ [common_params], :];
met_irr_best = met_irr_good[met_irr_good.scen .∈ [common_params], :];

scatter(met_ctr_best.swc_nse_com, met_irr_best.swc_nse_com, 
    xlabel="NSE Combined Control", ylabel="NSE Combined Irrigation", legend=false)

# best control scenario
scen_max_ctr, met_max_ctr = met_best_scen(met_ctr_good, :met_com);
# best irrigation scenario
scen_max_irr, met_max_irr = met_best_scen(met_irr_good, :met_com);

# parameter values for the best performing scenario
par_ctr_best = par_ctr[scen_max_ctr, :];
par_ctr_best

par_irr_best = par_irr[scen_max_irr, :];
par_irr_best


# parameter relationships
par_plots_ctr = par_plot(par_ctr, met_ctr, met_y="swc_nse80");
par_plots_irr = par_plot(par_irr, met_irr, met_y="swc_nse80");

plot(par_plots_ctr..., size=(1000,1000), layout=(4,5), legend=false, titlefontsize=8, guidefontsize=6)

# calculate K-S statistic to determine sensitive parameters
ks_stat_ctr, ks_plots_ctr = KS_plot(par_ctr, met_ctr);
ks_stat_irr, ks_plots_irr = KS_plot(par_irr, met_irr);

plot(ks_plots_irr..., size=(1000,1000), layout=(4,5), legend=false, titlefontsize=8, guidefontsize=6)
# behavioral is blue, non-behavioral is red



# plotting recipes

function ind_par_plot(par_ctr, par_irr, ind)
    p=density(par_ctr[!,ind], label="Control", color=:red, fillrange=0, alpha=.5);
    p=density!(par_irr[!,ind], label="Irrigation", color=:blue, fillrange=0, alpha=.5);
    return p
end

p1 = ind_par_plot(par_ctr_good, par_irr_good, 12);
p1=xlabel!(p1, "Maximum stomatal conductance (cm/s)");
plot(p1, size=(400,500), legend=:top)

p2 = ind_par_plot(par_ctr_good, par_irr_good, 13);
p2=xlabel!(p2,"VPD at 50% stomatal closure (kPa)");
plot(p2, size=(400,500), legend=:top)

p3 = ind_par_plot(par_ctr_good, par_irr_good, 15);
p3=xlabel!(p3,"Critical leaf water potential (MPa)");
plot(p3, size=(400,500), legend=:top)