# sensitivity analysis on LWFBrook90.jl

using CSV, DataFrames, Dates
using HypothesisTests
using StatsPlots, Measures, Plots; gr()

# function to filter out scenarios which produced error
function filter_error(met)
    println("Filtering out $(sum(met.twd_pd_cor .== 0 .&& met.twd_md_cor .== 0)) scenarios out of total $(size(met, 1)) which failed to run.")
    return met[met.twd_pd_cor .!= 0 .&& met.twd_md_cor .!= 0, :]
end

# function for density plot of metrics
function density_plot(met)
    density(met.twd_pd_cor, label="TWD_pd")
    density!(met.twd_md_cor, label="TWD_md")
    #density!(met.lwp_pd_cor, label="LWP_pd")
    #density!(met.lwp_md_cor, label="LWP_md")
end

# function to filter metrics for behavioral runs
function behavioral_met(met)

    # control
    #return(met[met.twd_pd_cor .< -0.5 .&&
    #           met.twd_md_cor .< -0.5, :])
               
    # irr stop
    return(met[met.twd_pd_cor .< -0.8 .&&
               met.twd_md_cor .< -0.8, :])
    
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
function par_plot(par, met; met_y="met_com", behave=true)

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
        mc="black", msc="black", ms = 3, legend=false)
        title!(p, "$(names(par)[i]) vs $(met_y)")
        push!(par_plots, p)
    end
    
    return par_plots
end

# function to plot all parameter interactions with desired metric
function par2_plot(par, met; met_y="met_com", behave=true)

    if behave
        # only plot behavioral parameter sets
        par, _ = sep_params(par, met);
        met = behavioral_met(met);
        marker_size = 5;
    else
        par = par[met.scen, :];
        marker_size = 3;
    end

    np = size(par, 2);

    pair_idx = [(i, j) for i in 1:(np - 1) for j in (i + 1):np];
    n_pairs = length(pair_idx);
    ncols = ceil(Int, sqrt(n_pairs));
    nrows = ceil(Int, n_pairs / ncols);
    metric_vals = met[!, met_y];
    metric_lim = extrema(metric_vals);
    if metric_lim[1] == metric_lim[2]
        metric_lim = (metric_lim[1] - 1, metric_lim[2] + 1);
    end
    par_plots = Any[];

    for (k, (i, j)) in enumerate(pair_idx)
        p = scatter(
            par[!, i],
            par[!, j],
            marker_z=metric_vals,
            ms=marker_size,
            markerstrokewidth=0,
            xlabel=names(par)[i],
            ylabel=names(par)[j],
            clims=metric_lim,
            colorbar=(k == n_pairs),
            legend=false,
            left_margin=4mm
        )
        push!(par_plots, p)
    end

    return plot(
        par_plots...,
        size=(350 * ncols, 300 * nrows),
        layout=(nrows, ncols),
        plot_title="Parameter interactions colored by $(met_y)",
        titlefontsize=8,
        guidefontsize=7
    )
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
    met.twd_com = sqrt.((met.twd_pd_cor .^ 2 + met.twd_md_cor .^ 2) / 2);
    #met.lwp_com = sqrt.((met.lwp_pd_cor .^ 2 + met.lwp_md_cor .^ 2) / 2);
    #met.met_com = sqrt.((met.lwp_com .^ 2 + met.twd_com .^ 2) / 2);
    return nothing
end

# function to find best scenario
function met_best_scen(met, metric=:met_com)
    # find the index of the maximum value
    max_idx = argmax(met[!, metric]);
    best_scen = met.scen[max_idx];
    
    return best_scen, met[max_idx, :]
end

# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_pfyn_ctr_cap_20260904.csv", DataFrame);
#met_irr = CSV.read("LWFBcal_output/metrics_pfyn_irr_cap_20260904.csv", DataFrame);
met_irst = CSV.read("LWFBcal_output/metrics_pfyn_irst_cap_20260904.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_pfyn_ctr_cap_20260904.csv", DataFrame);
#par_irr = CSV.read("LWFBcal_output/param_pfyn_irr_cap_20260904.csv", DataFrame);
par_irst = CSV.read("LWFBcal_output/param_pfyn_irst_cap_20260904.csv", DataFrame);

# filter out scenarios which produced an error
met_ctr = filter_error(met_ctr);
#met_irr = filter_error(met_irr);
met_irst = filter_error(met_irst);


# exploratory
density_plot(met_ctr)
scatter(met_ctr.plfl_mean, met_ctr.plfl_max, xlabel="Mean Daily Plant Discharge > 0 (mm)", 
    ylabel="Max Daily Plant Discharge (mm)", legend=false)

# calculate pFs
met_ctr.min_pd_pF = log10.(-10 * met_ctr.min_pd_psi);
#met_irr.min_md_pF = log10.(-10 * met_irr.min_md_psi);
met_irst.min_md_pF = log10.(-10 * met_irst.min_md_psi);

# parameter interactions
par = par_ctr;
met = met_ctr;
par2_plot(par, met, met_y="min_pd_pF", behave=false)
par2_plot(par, met, met_y="min_md_pF", behave=false)
par2_plot(par, met, met_y="plfl_mean", behave=false)
par2_plot(par, met, met_y="plfl_max", behave=false)
par2_plot(par, met, met_y="twd_pd_cor", behave=false)
par2_plot(par, met, met_y="twd_md_cor", behave=false)

# add combined metrics
met_comb!(met_ctr);
met_comb!(met_irr);
met_comb!(met_irst);

# filter metrics for behavioral runs
met_ctr_good = behavioral_met(met_ctr);
#met_irr_good = behavioral_met(met_irr);
met_irst_good = behavioral_met(met_irst);

println("$(size(met_ctr_good, 1)) behavioral parameter sets out of total $(size(met_ctr, 1)) in control scenario.")
#println("$(size(met_irr_good, 1)) behavioral parameter sets out of total $(size(met_irr, 1)) in irrigation scenario.")
println("$(size(met_irst_good, 1)) behavioral parameter sets out of total $(size(met_irst, 1)) in irrigation stop scenario.")

density_plot(met_ctr_good)

# compare metrics
met_plot(met_irst, :twd_pd_cor, :twd_md_cor)
met_plot(met_ctr, [:twd_pd_cor, :twd_pd_cor, :lwp_pd_cor, :twd_md_cor], [:twd_md_cor, :lwp_pd_cor, :lwp_md_cor, :lwp_md_cor])
met_plot(met_ctr_good, [:twd_pd_cor, :twd_pd_cor, :lwp_pd_cor, :twd_md_cor], [:twd_md_cor, :lwp_pd_cor, :lwp_md_cor, :lwp_md_cor])

# best control scenario
scen_max, met_max = met_best_scen(met_good);
# parameter values for the best performing scenario
par_best = par[scen_max, :];
par_best

# parameter relationships
par_plots = par_plot(par_ctr, met_ctr, met_y="twd_md_cor");
par_plots = par_plot(par_irst, met_irst, met_y="twd_md_cor");

plot(par_plots..., size=(1000,1000), layout=(2,2), legend=false, titlefontsize=12, guidefontsize=12)

# parameter interactions
par2_plot(par, met, met_y="min_md_psi")
par2_plot(par, met, met_y="plfl_mean")
par2_plot(par, met, met_y="plfl_max")
par2_plot(par, met, met_y="twd_pd_cor")
par2_plot(par, met, met_y="twd_md_cor")
par2_plot(par, met, met_y="lwp_pd_cor")
par2_plot(par, met, met_y="lwp_md_cor")
par2_plot(par, met)

# calculate K-S statistic to determine sensitive parameters
ks_stat_ctr, ks_plots_ctr = KS_plot(par_ctr, met_ctr);
#ks_stat_irr, ks_plots_irr = KS_plot(par_irr, met_irr);
ks_stat_irst, ks_plots_irst = KS_plot(par_irst, met_irst);

plot(ks_plots_irst..., size=(1200,1200), layout=(2,2), legend=false, titlefontsize=8, guidefontsize=6)
# behavioral is blue, non-behavioral is red

