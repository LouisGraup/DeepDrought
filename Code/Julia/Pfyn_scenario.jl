## run best performing scenario from LWFBrook90.jl calibration
## for all Pfynwald scenarios
## and examine output against observed data

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions, RCall;
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Plots; gr()
R"""
library(tidyverse)
"""

include("run_LWFB90_param.jl");

# function to retrieve best parameter set for control and irrigation scenarios
function par_best(met, par)
    
    # add combined NSE metrics
    met.swc_nse_com = (met.swc_nse10 + met.swc_nse40 + met.swc_nse60 + met.swc_nse80) / 4;
    met.swp_nse_com = (met.swp_nse10 + met.swp_nse80) / 2;

    met.nse_com = (met.swc_nse_com + met.swp_nse_com) / 2;

    # best overall scenario
    max_idx = argmax(met.nse_com);
    scen_best = met.scen[max_idx];

    par_best = par[scen_best, :];

    return par_best, scen_best
end

function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function get_swc(sim; shape = "long")
    # retrieve soil water potential data from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = days);
    z.date = Date.(dates_out);

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="VWC");
        z_long.depth = parse.(Int, map(x -> x[(end-4):(end-2)], z_long.depth)) .÷ 10; # convert depth to cm from var name
        return z_long
    else
        return z
    end

end

function get_swp(sim; shape="long")
    # retrieve soil water potential data from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

    z = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = days);
    z.date = Date.(dates_out);

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="SWP");
        z_long.depth = parse.(Int, map(x -> x[(end-4):(end-2)], z_long.depth)) .÷ 10; # convert depth to cm from var name
        
        return z_long
    else
        return z
    end

end

function get_sap(sim)
    # retrieve soil water potential data from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

    z = get_fluxes(sim);
    z.date = Date.(dates_out);
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

function ann_trans(sim)
    # calculates annual transpiration
    z = get_sap(sim);

    z.year = year.(z.date);

    z_ann = combine(groupby(z, :year), :trans => sum);

    return z_ann
end

function swc_comp(sim, obs_swc)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water content data

    z_theta = get_swc(sim);

    filter!(:date => >(obs_swc.date[1]), z_theta);

    z_theta.src .= "sim"; # add source column

    select!(obs_swc, :date, :depth, :VWC); # remove extra columns
    obs_swc.src .= "obs"; # add source column

    # combine observed and simulated data
    swc_comp = [obs_swc; z_theta];

    return swc_comp

end

function swp_comp(sim, obs_swp)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_psi = get_swp(sim);

    filter!(:date => >(obs_swp.date[1]), z_psi);

    z_psi.src .= "sim"; # add source column

    obs_swp_long = stack(obs_swp, Not(:date, :meta), variable_name = "depth", value_name="SWP");
    select!(obs_swp_long, :date, :depth, :SWP);
    obs_swp_long.depth = parse.(Int, map(x -> x[(end-3):(end-2)], obs_swp_long.depth)); # convert depth from var name
    obs_swp_long.src .= "obs"; # add source column

    # combine observed and simulated data
    swp_comp = [obs_swp_long; z_psi];

    return swp_comp

end

function sap_comp(sim, obs_sap)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_trans = get_sap(sim);

    filter!(:date => >(obs_sap.date[1]), z_trans);
    filter!(:date => <(Date(2023, 1, 1)), z_trans);

    obs_sap = select(obs_sap, Not(:meta));

    sap_comp = leftjoin(z_trans, obs_sap, on = :date);

    return sap_comp

end

function get_RWU_centroid(sim)
    # borrow code from LWFBrook90 package
    solu = sim.ODESolution;

    days_to_read_out_d = unique(round.(solu.t));

    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2;

    # Compute RWU centroid
    rows_RWU_mmDay  = reduce(hcat, [solu(t).TRANI.mmday   for t in days_to_read_out_d]);

    RWU_percent = rows_RWU_mmDay ./ sum(rows_RWU_mmDay; dims = 1);
    #RWUcentroidLabel = "mean RWU depth"
    if (any(RWU_percent .< 0))
        #@warn "Some root water outfluxes detected. Centroid of RWU is  based only on uptakes."
        rows_RWU_mmDay_onlyUptake = ifelse.(rows_RWU_mmDay.>0,rows_RWU_mmDay, 0);
        RWU_percent_onlyUptake = rows_RWU_mmDay_onlyUptake ./ sum(rows_RWU_mmDay_onlyUptake; dims = 1);
        RWU_percent = RWU_percent_onlyUptake;
        
        #RWUcentroidLabel = "mean RWU depth\n(based on uptake only)"
    end

    row_RWU_centroid_mm = sum(RWU_percent .* y_center; dims=1);

    col_RWU_centroid_mm = reshape(row_RWU_centroid_mm, :);
    
    return col_RWU_centroid_mm
end

function plot_monthly_water_partitioning(df_partitioning_monthly)
    # color palette
    color_palette = reverse([
            "Transpiration deficit" => :red2,
            "Actual transpiration" => :darkolivegreen2,
            "Interception loss" => :forestgreen,
            "Soil evaporation" => :khaki3,
            #"Snow sublimation" => :white,
            "Runoff" => :lightskyblue,
            "Drainage" => :darkblue,
            # "ETa" => :black,
            # "P2" => :darkblue,
            "Precipitation" => :black,
            # "Swat" => :brown
        ]);

    # Preprocess
    df_part_mth = copy(df_partitioning_monthly);
    rename!(df_part_mth, :Td => "Transpiration deficit",
            :Ta => "Actual transpiration",
            :Einterception => "Interception loss",
            :Esoil => "Soil evaporation",
            :Esnow => "Snow sublimation",
            :R => "Runoff",
            :D => "Drainage",
            :Precip => "Precipitation");

    # Preprocess
    df_part_monthly_forMakie = @chain df_part_mth begin
        stack(Not([:year, :month, :nrow, :date]))
        # only keep variables we need
        @subset(:variable .∈ ([first(pair) for pair in color_palette],))
        # make categorical
        @transform :variable = CategoricalArrays.categorical(:variable, levels = [first(pair) for pair in color_palette])
        @transform :variable_code = CategoricalArrays.levelcode.(:variable)
        # Remove fluxes that were not computed (e.g. removes runoff)
        @subset(:value .!= 0.0)
        end

    # Plot
    aog_monthly = mapping(
            :date => "",
            :value => "Water flux per month (mm)",
            stack = :variable => "Water fluxes",
            color = :variable => "") *
        # bar plot of fluxes
        data(@subset(df_part_monthly_forMakie, :variable .!= "Precipitation")) * visual(BarPlot) +
        mapping(
            :date => "",
            :value => "Water flux per month (mm)",
            color = :variable => "") *
        # line plot of precip input
        data(@subset(df_part_monthly_forMakie, :variable .== "Precipitation")) * visual(Lines)
        
    xticks = sort(unique(Dates.floor.(df_part_monthly_forMakie.date, Dates.Month(6))))

    aog_draw = draw(aog_monthly, scales(Color = (; palette = color_palette)),
        axis = (; ygridvisible = true))
                #xticks = AlgebraOfGraphics.datetimeticks((x -> Dates.format(x, "u\nY")), (Date.(xticks)))))
    return aog_draw
end

function plot_yearly_water_partitioning(df_partitioning_yearly)
    color_palette = reverse([
            "Transpiration deficit" => :red2,
            "Actual transpiration" => :darkolivegreen2,
            "Interception loss" => :forestgreen,
            "Soil evaporation" => :khaki3,
            #"Snow sublimation" => :white,
            "Runoff" => :lightskyblue,
            "Drainage" => :darkblue,
            # "ETa" => :black,
            # "P2" => :darkblue,
            "Precipitation" => :black,
            # "Swat" => :brown
        ]);

    # Preprocess
    df_part_yr = copy(df_partitioning_yearly);
    rename!(df_part_yr, :Td => "Transpiration deficit",
            :Ta => "Actual transpiration",
            :Einterception => "Interception loss",
            :Esoil => "Soil evaporation",
            :Esnow => "Snow sublimation",
            :R => "Runoff",
            :D => "Drainage",
            :Precip => "Precipitation");
    
    df_part_yearly_forMakie = @chain df_part_yr begin
        stack(Not([:year, :nrow, :date]))
        # only keep variables we need
        @subset(:variable .∈ ([first(pair) for pair in color_palette],))
        # make categorical
        @transform :variable = categorical(:variable, levels = [first(pair) for pair in color_palette])
        @transform :variable_code = levelcode.(:variable)
        # Remove fluxes that were not computed (e.g. removes runoff)
        @subset(:value .!= 0.0)
        end

    # Plot
    aog_yearly = mapping(
            :year => "",
            :value => "Water flux per year (mm)",
            stack = :variable => "Water fluxes",
            color = :variable => "") *
        # bar plot of fluxes
        data(@subset(df_part_yearly_forMakie, :variable .!= "Precipitation")) * visual(BarPlot) +
        mapping(
            :year => "",
            :value => "Water flux per year (mm)",
            color = :variable => "") *
        # line plot of precip input
        data(@subset(df_part_yearly_forMakie, :variable .== "Precipitation")) * visual(Lines)
        
    aog_draw = draw(aog_yearly, scales(Color = (; palette = color_palette)))
        
    return aog_draw
end

# calculate NSE for soil water potential
function obj_fun_swp(sim, obs)

    z = get_swp(sim, shape="wide");

    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :SWP_10cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :SWP_80cm]]);

    # match simulated data to available dates for each depth
    sim_10cm = z[z.date .∈ [obs_10cm.date], :psi_kPa_100mm];
    sim_80cm = z[z.date .∈ [obs_80cm.date], :psi_kPa_800mm];

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    # calculate NSE
    nse10 = NSE(sim_10cm, obs_10cm.SWP_10cm);
    nse80 = NSE(sim_80cm, obs_80cm.SWP_80cm);

    return nse10, nse80
end

# calculate NSE for soil water content
function obj_fun_swc(sim, obs)

    z = get_swc(sim, shape="wide");

    # process observed data
    obs = unstack(obs, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
    sort!(obs, :date); # sort by date

    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :VWC_10cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :VWC_80cm]]);

    # match simulated data to available dates for each depth
    sim_10cm = z[z.date .∈ [obs_10cm.date], :theta_m3m3_100mm];
    sim_80cm = z[z.date .∈ [obs_80cm.date], :theta_m3m3_800mm];

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    # calculate NSE
    nse10 = NSE(sim_10cm, obs_10cm.VWC_10cm);
    nse80 = NSE(sim_80cm, obs_80cm.VWC_80cm);

    return nse10, nse80
end

# calculate correlation coefficient between sap flow and transpiration data
function obs_fun_sap(sap_comp)

    # remove missing values
    sap_comp = dropmissing(sap_comp);

    cc = cor(sap_comp.trans, sap_comp.sfd);

    return cc

end

# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250814.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250814.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_ctr_20250814.csv", DataFrame);
par_irr = CSV.read("LWFBcal_output/param_irr_20250814.csv", DataFrame);

par_ctr_best, scen_ctr_best = par_best(met_ctr, par_ctr);
par_irr_best, scen_irr_best = par_best(met_irr, par_irr);

met_ctr[scen_ctr_best, :]
met_irr[scen_irr_best, :]

## behavioral data
# soil water content
obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
#filter!(:date => >(Date(2015, 1, 1)), obs_swc); # filter out early dates
filter!(:date => <(Date(2024, 1, 1)), obs_swc); # filter out late dates

# soil water potential
obs_swp = CSV.read("../../Data/Pfyn/PFY_swpc.csv", DataFrame);
filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
filter!(:date => <(Date(2024, 1, 1)), obs_swp); # filter out late dates

# sap flow
obs_sap = CSV.read("../../Data/Pfyn/PFY_sap.csv", DataFrame);

# separate control and irrigation scenarios
obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
obs_swc_irr = obs_swc[obs_swc.meta .== "irrigated", :]; # select irrigation treatment
obs_swc_irst = obs_swc[obs_swc.meta .== "irrigation_stop", :]; # select irrigation stop treatment

obs_swp_ctr = obs_swp[obs_swp.meta .== "control", :]; # select control treatment
obs_swp_irr = obs_swp[obs_swp.meta .== "irrigation", :]; # select irrigation treatment
obs_swp_irst = obs_swp[obs_swp.meta .== "irrigation_stop", :]; # select irrigation stop treatment

obs_sap_ctr = obs_sap[obs_sap.meta .== "Control", :]; # select control treatment
obs_sap_irr = obs_sap[obs_sap.meta .== "Irrigation", :]; # select irrigation treatment
obs_sap_irst = obs_sap[obs_sap.meta .== "Irrigation Stop", :]; # select irrigation stop treatment

# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2014, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_irrigiso_ambient/", "pfynwald", "LWFB_testrun/irrigation/", irrig=true);
sim_irst = run_LWFB90_param(par_ctr_best, Date(2014, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_irr_stop/", "pfynwald", "LWFB_testrun/irr_stop/");

## combine observed and simulated data
# soil water content
swc_comp_ctr = swc_comp(sim_ctr, obs_swc_ctr);
swc_comp_irr = swc_comp(sim_irr, obs_swc_irr);
swc_comp_irst = swc_comp(sim_irst, obs_swc_irst);

# soil water potential
swp_comp_ctr = swp_comp(sim_ctr, obs_swp_ctr);
swp_comp_irr = swp_comp(sim_irr, obs_swp_irr);
swp_comp_irst = swp_comp(sim_irst, obs_swp_irst);

# sap flow
sap_comp_ctr = sap_comp(sim_ctr, obs_sap_ctr);
sap_comp_irr = sap_comp(sim_irr, obs_sap_irr);
sap_comp_irst = sap_comp(sim_irst, obs_sap_irst);

# compare sap flow against modelled transpiration
obs_fun_sap(sap_comp_ctr)
obs_fun_sap(sap_comp_irr)
obs_fun_sap(sap_comp_irst)

# validation metrics for soil water potential
obj_fun_swp(sim_ctr, obs_swp_ctr[obs_swp_ctr.date .>= Date(2021, 1, 1), :])
obj_fun_swp(sim_irr, obs_swp_irr[obs_swp_irr.date .>= Date(2021, 1, 1), :])

# compare irrigation stop scenario against observations
obj_fun_swp(sim_irst, obs_swp_irst)
obj_fun_swc(sim_irst, obs_swc_irst)

# individually by year
obs_swp_irst.year = year.(obs_swp_irst.date);
obs_swc_irst.year = year.(obs_swc_irst.date);
for year in unique(obs_swp_irst.year)
    # loop through each year and calculate NSE
    nse10, nse80 = obj_fun_swp(sim_irst, obs_swp_irst[obs_swp_irst.year .== year, :]);
    #nse10, nse80 = obj_fun_swc(sim_irst, obs_swc_irst[obs_swc_irst.year .== year, :]);
    println("Year: $year, NSE10: $nse10, NSE80: $nse80");
end


## plot observed and simulated data using ggplot

# soil water potential
R"""
rdf = $swp_comp_irst
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) + theme_bw() +
    labs(x="", y="SMP (kPa)", color="Source", title="Soil Water Potential Comparison for Irrigation Stop") +
    theme(plot.title = element_text(hjust=0.5))
"""

R"""
rdf = $swp_comp_ctr
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) + theme_bw() +
    labs(x="", y="SMP (kPa)", color="Source", title="Soil Water Potential Comparison for Control Scenario") +
    theme(plot.title = element_text(hjust=0.5), legend.text=element_text(size=12),legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"),axis.text=element_text(size=12), axis.title=element_text(size=14))
"""

R"""
rdf = $swp_comp_irr
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation")
"""

# soil water content
R"""
rdf = $swc_comp_ctr
ggplot(rdf, aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Content Comparison for Control")
"""

R"""
rdf = $swc_comp_irr
ggplot(rdf, aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Content Comparison for Irrigation")
"""

R"""
rdf = $swc_comp_irst
ggplot(filter(rdf, depth %in% c(10,80)), aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Content Comparison for Irrigation Stop")
"""

# sap flow
R"""
rdf = $sap_comp_ctr
ggplot(rdf, aes(x=date, y=trans, color="trans")) + geom_point() +
    geom_point(aes(date, sfd, color="sap")) +
    labs(title="Sap Flow Comparison for Control")
"""

R"""
rdf = $sap_comp_ctr
ggplot(rdf, aes(x=trans, y=sfd)) + geom_point() + geom_smooth(method='lm') +
    labs(title="Sap Flow Comparison for Control") + theme_bw()
"""

R"""
rdf = $sap_comp_irr
ggplot(rdf, aes(x=date, y=trans, color="trans")) + geom_point() + geom_smooth( aes(x=date, y=trans, color="trans")) +
    geom_point(aes(date, sfd/1.8, color="sap/1.8")) + geom_smooth(aes(date, sfd/1.8, color="sap/1.8")) + theme_bw() +
    labs(x="", title="Sap Flow Comparison for Irrigation")
"""

R"""
rdf = $sap_comp_irr
ggplot(rdf, aes(x=trans, y=sfd)) + geom_point() + geom_smooth(method='lm') +
    labs(title="Sap Flow Comparison for Irrigation") + theme_bw()
"""

R"""
rdf = $sap_comp_irst
ggplot(rdf, aes(x=date, y=trans, color="trans")) + geom_point() +
    geom_point(aes(date, sfd, color="sap")) +
    labs(title="Sap Flow Comparison for Irrigation Stop")
"""

R"""
rdf = $sap_comp_irst
ggplot(rdf, aes(x=trans, y=sfd)) + geom_point() + geom_smooth(method='lm') +
    labs(title="Sap Flow Comparison for Irrigation Stop") + theme_bw()
"""

# compare soil water potential across scenarios
ctr_swp = get_swp(sim_ctr);
ctr_swp.scen .= "Control";

irr_swp = get_swp(sim_irr);
irr_swp.scen .= "Irrigation";

irst_swp = get_swp(sim_irst);
irst_swp.scen .= "Irrigation Stop";

swp_comp_scen = [ctr_swp; irr_swp; irst_swp];

R"""
rdf = $swp_comp_scen
ggplot(filter(rdf, date>="2014-01-01", date<"2020-01-01"), aes(x=date, y=SWP, color=scen)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) + theme_bw() +
    labs(x="", y="SMP (kPa)", color="Scenario", title="Modelled Soil Water Potential Comparison across Scenarios") +
    theme(legend.position="right",legend.text=element_text(size=12),legend.title=element_text(size=14),plot.title = element_text(hjust=0.5), strip.text=element_text(size=12, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=14)) +
    scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))
"""


# compare RWU depth across scenarios

days, dates_out = get_dates(sim_ctr);

rwu_comp = DataFrame(date=dates_out);

rwu_comp.control = runmean(get_RWU_centroid(sim_ctr), 14);
rwu_comp.irrigation = runmean(get_RWU_centroid(sim_irr), 14);
rwu_comp.irrigation_stop = runmean(get_RWU_centroid(sim_irst), 14);

rwu_comp = stack(rwu_comp, Not(:date), variable_name="treatment", value_name="RWU");
rwu_comp.RWU = replace(rwu_comp.RWU, NaN=>missing);

R"""
rdf = $rwu_comp
ggplot(filter(rdf, date<"2020-01-01"), aes(x=date, y=RWU, color=treatment)) + geom_line() +
    labs(x="", title="RWU Depth Comparison across Scenarios") + theme_bw() +
    scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))
"""


# compare annual transpiration across scenarios
ctr_yr = ann_trans(sim_ctr);
ctr_yr.scen .= "Control";

irr_yr = ann_trans(sim_irr);
irr_yr.scen .= "Irrigation";

irst_yr = ann_trans(sim_irst);
irst_yr.scen .= "Irrigation Stop";

yr_comp = [ctr_yr; irr_yr; irst_yr];

R"""
rdf = $yr_comp
ggplot(filter(rdf, year>=2003), aes(x=year, y=trans_sum, group=scen, fill=scen)) + geom_col(position="dodge") +
    labs(title="Transpiration Comparison Across Scenarios")
"""

# compare RWU against transpiration

df_rwu_tran = get_sap(sim_ctr);
df_rwu_tran.RWU = get_RWU_centroid(sim_ctr);
df_rwu_tran.month = month.(df_rwu_tran.date);

draw(
    data(df_rwu_tran[df_rwu_tran.month .> 5 .&& df_rwu_tran.month .< 12, :])*
    mapping(:trans, :RWU, color=:month => nonnumeric)*visual(Scatter),
    scales(Color = (; palette = from_continuous(:seaborn_colorblind6))),
    axis = (; yreversed = true)
)

# water balance modelling

sim_ctr_wb = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/", new_folder=false, watbal=true);
sim_irr_wb = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_irrigation_ambient/", "pfynwald", "LWFB_testrun/irrigation/", new_folder=false, watbal=true);
sim_irst_wb = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_irr_stop/", "pfynwald", "LWFB_testrun/irr_stop/", new_folder=false, watbal=true);

# water partitioning

# control
wb_d_ctr, wb_m_ctr, wb_y_ctr = get_water_partitioning(sim_ctr_wb);
wb_m_ctr = wb_m_ctr[wb_m_ctr.year .== 2014, :]; # filter out single years
wb_y_ctr = wb_y_ctr[wb_y_ctr.year .>= 2003, :]; # filter out early years

wb_yr_ctr_plot = plot_yearly_water_partitioning(wb_y_ctr);
wb_yr_ctr_plot

# irrigation
wb_d_irr, wb_m_irr, wb_y_irr = get_water_partitioning(sim_irr_wb);
wb_mth_irr = wb_m_irr[wb_m_irr.year .== 2014, :]; # filter out single year
wb_yr_irr = wb_y_irr[wb_y_irr.year .>= 2003, :]; # filter out early years

wb_mth_irr_plot = plot_monthly_water_partitioning(wb_mth_irr);

wb_yr_irr_plot = plot_yearly_water_partitioning(wb_yr_irr);
wb_yr_irr_plot

# irrigation stop
wb_d_irst, wb_m_irst, wb_y_irst = get_water_partitioning(sim_irst_wb);
wb_mth_irst = wb_m_irst[wb_m_irst.year .== 2014, :]; # filter out single year
wb_yr_irst = wb_y_irst[wb_y_irst.year .>= 2003, :]; # filter out early years

wb_mth_irst_plot = plot_monthly_water_partitioning(wb_mth_irst);

wb_yr_irst_plot = plot_yearly_water_partitioning(wb_yr_irst);
wb_yr_irst_plot