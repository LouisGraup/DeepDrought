## run best performing scenario from LWFBrook90.jl calibration
## for all Pfynwald scenarios
## and examine output against observed data

using CSV, DataFrames, Dates, Statistics, RCall;
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

function get_swc(sim)
    # retrieve soil water potential data from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = days);
    z.date = Date.(dates_out);

    return z
end

function get_swp(sim)
    # retrieve soil water potential data from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

    z = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = days);
    z.date = Date.(dates_out);

    return z
end

function swc_comp(sim, obs_swc)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water content data

    z_theta = get_swc(sim);

    filter!(:date => >(obs_swc.date[1]), z_theta);
    select!(z_theta, Not(:time));

    # reshape data
    z_long = stack(z_theta, Not(:date), variable_name = "depth", value_name="VWC");
    z_long.depth = parse.(Int, map(x -> x[(end-4):(end-2)], z_long.depth)) .÷ 10; # convert depth to cm from var name
    z_long.src .= "sim"; # add source column

    select!(obs_swc, :date, :depth, :VWC); # remove extra columns
    obs_swc.src .= "obs"; # add source column

    # combine observed and simulated data
    swc_comp = [obs_swc; z_long];

    return swc_comp

end

function swp_comp(sim, obs_swp)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_psi = get_swp(sim);

    filter!(:date => >(obs_swp.date[1]), z_psi);
    select!(z_psi, Not(:time));

    # reshape data
    z_long = stack(z_psi, Not(:date), variable_name = "depth", value_name="SWP");
    z_long.depth = parse.(Int, map(x -> x[(end-4):(end-2)], z_long.depth)) .÷ 10; # convert depth to cm from var name
    z_long.src .= "sim"; # add source column

    obs_swp_long = stack(obs_swp, Not(:date, :meta), variable_name = "depth", value_name="SWP");
    select!(obs_swp_long, :date, :depth, :SWP);
    obs_swp_long.depth = parse.(Int, map(x -> x[(end-3):(end-2)], obs_swp_long.depth)); # convert depth from var name
    obs_swp_long.src .= "obs"; # add source column

    # combine observed and simulated data
    swp_comp = [obs_swp_long; z_long];

    return swp_comp

end

# calculate NSE for soil water potential
function obj_fun_swp(sim, obs)

    z = get_swp(sim);

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

# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250715.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250715.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_20250715.csv", DataFrame);

par_ctr_best, scen_ctr_best = par_best(met_ctr, par);
par_irr_best, scen_irr_best = par_best(met_irr, par);

met_ctr[scen_ctr_best, :]
met_irr[scen_irr_best, :]

## behavioral data
obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
filter!(:date => >(Date(2004, 1, 1)), obs_swc); # filter out early dates

obs_swp = CSV.read("../../Data/Pfyn/PFY_swpc.csv", DataFrame);
filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
filter!(:date => <(Date(2021, 1, 1)), obs_swp); # filter out late dates

# separate control and irrigation scenarios
obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
obs_swc_irr = obs_swc[obs_swc.meta .== "irrigated", :]; # select irrigation treatment

obs_swp_ctr = obs_swp[obs_swp.meta .== "control", :]; # select control treatment
obs_swp_irr = obs_swp[obs_swp.meta .== "irrigated", :]; # select irrigation treatment
obs_swp_irst = obs_swp[obs_swp.meta .== "irrigation_stop", :]; # select irrigation stop treatment


# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/Pfyn_irrigation_ambient/", "pfynwald", "LWFB_testrun/irrigation/");
sim_irst = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/Pfyn_irr_stop/", "pfynwald", "LWFB_testrun/irr_stop/");

# combine observed and simulated data
# soil water content
swc_comp_ctr = swc_comp(sim_ctr, obs_swc_ctr);
swc_comp_irr = swc_comp(sim_irr, obs_swc_irr);

swp_comp_ctr = swp_comp(sim_ctr, obs_swp_ctr);
swp_comp_irr = swp_comp(sim_irr, obs_swp_irr);
swp_comp_irst = swp_comp(sim_irst, obs_swp_irst);

# compare irrigation stop scenario against observations
obj_fun_swp(sim_irst, obs_swp_irst)

# plot observed and simulated data using ggplot
R"""
rdf = $swp_comp_irst
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation Stop")
"""

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
rdf = $swp_comp_ctr
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Control")
"""

R"""
rdf = $swp_comp_irr
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation")
"""