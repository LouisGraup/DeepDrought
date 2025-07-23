## run best performing scenario from LWFBrook90.jl calibration
## for all Pfynwald VPDrought scenarios
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

    return z
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


# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250722.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250722.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_ctr_20250722.csv", DataFrame);
par_irr = CSV.read("LWFBcal_output/param_irr_20250722.csv", DataFrame);

par_ctr_best, scen_ctr_best = par_best(met_ctr, par_ctr);
par_irr_best, scen_irr_best = par_best(met_irr, par_irr);

met_ctr[scen_ctr_best, :]
met_irr[scen_irr_best, :]


## behavioral data
obs = CSV.read("../../Data/Pfyn/PFY_VPD_swp_swc.csv", DataFrame);

obs_ctr = obs[obs.treatment .== "control", :];
obs_irr = obs[obs.treatment .== "irrigation", :];
obs_irr_vpd = obs[obs.treatment .== "irrigation_vpd", :];
obs_drt = obs[obs.treatment .== "roof", :];
obs_drt_vpd = obs[obs.treatment .== "roof_vpd", :];