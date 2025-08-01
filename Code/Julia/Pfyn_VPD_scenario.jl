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

function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function get_swc(sim; shape = "long")
    # retrieve soil water content data from sim
    
    days, dates_out = get_dates(sim);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 800, 1200], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    # rename columns
    rename!(z, ["0.1", "0.8", "1.2", "date"]);

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="VWC");
        z_long.depth = parse.(Float64, z_long.depth);

        return z_long
    else
        return z
    end

end

function get_swp(sim; shape="long")
    # retrieve soil water potential data from sim
    
    days, dates_out = get_dates(sim);

    z = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800, 1200], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));
    
    # rename columns
    rename!(z, ["0.1", "0.8", "1.2", "date"]);

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="SWP");
        z_long.depth = parse.(Float64, z_long.depth);
        
        return z_long
    else
        return z
    end

end

function get_sap(sim)
    # retrieve soil water potential data from sim
    
    days, dates_out = get_dates(sim);

    z = get_fluxes(sim);
    z.date = dates_out;
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

function met_comp(sim, obs)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil moisture data

    z_theta = get_swc(sim);
    z_psi = get_swp(sim);

    z_sim = innerjoin(z_theta, z_psi, on = [:date, :depth]);

    filter!(:date => >(obs.date[1]), z_sim);

    z_sim.src .= "sim"; # add source column

    obs = select(obs, Not(:treatment)); # remove extra columns
    obs.src .= "obs"; # add source column

    # combine observed and simulated data
    met_comp = [obs; z_sim];

    return met_comp

end


# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250724.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250724.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_ctr_20250724.csv", DataFrame);
par_irr = CSV.read("LWFBcal_output/param_irr_20250724.csv", DataFrame);

par_ctr_best, scen_ctr_best = par_best(met_ctr, par_ctr);
par_irr_best, scen_irr_best = par_best(met_irr, par_irr);

met_ctr[scen_ctr_best, :]
met_irr[scen_irr_best, :]


## behavioral data
obs = CSV.read("../../Data/Pfyn/PFY_VPD_swp_swc.csv", DataFrame, missingstring="NA");
rename!(obs, :SWP_corr => :SWP);
obs.depth = -1 * obs.depth;

obs_ctr = obs[obs.treatment .== "control", :];
obs_irr = obs[obs.treatment .== "irrigation", :];
obs_irr_vpd = obs[obs.treatment .== "irrigation_vpd", :];
obs_drt = obs[obs.treatment .== "roof", :];
obs_drt_vpd = obs[obs.treatment .== "roof_vpd", :];


# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2025, 6, 30), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_VPD/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2025, 6, 30), "LWFBinput/Pfyn_irrigation_ambient/", "pfynwald", "LWFB_VPD/irr_amb/");
sim_irr_vpd = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2025, 6, 30), "LWFBinput/Pfyn_irrigation_VPD/", "pfynwald", "LWFB_VPD/irr_vpd/");
sim_drt = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2025, 6, 30), "LWFBinput/Pfyn_drought_ambient/", "pfynwald", "LWFB_VPD/drt_amb/");
sim_drt_vpd = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2025, 6, 30), "LWFBinput/Pfyn_drought_VPD/", "pfynwald", "LWFB_VPD/drt_vpd/");

ctr_comp = met_comp(sim_ctr, obs_ctr);
irr_comp = met_comp(sim_irr, obs_irr);
irr_vpd_comp = met_comp(sim_irr_vpd, obs_irr_vpd);
drt_comp = met_comp(sim_drt, obs_drt);
drt_vpd_comp = met_comp(sim_drt_vpd, obs_drt_vpd);


R"""
rdf = $ctr_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Control")
"""


R"""
rdf = $irr_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation")
"""


R"""
rdf = $irr_vpd_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation VPD")
"""


R"""
rdf = $drt_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Potential Comparison for Roof")
"""


R"""
rdf = $drt_vpd_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Roof VPD")
"""


R"""
rdf1 = $drt_comp
rdf2 = $drt_vpd_comp
ggplot(filter(rdf1, src=="sim"), aes(x=date, y=SWP, color="roof")) + geom_point(size=.5) +
    geom_point(data=filter(rdf2, src=="sim"), aes(x=date, y=SWP, color="roof_vpd"), size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Drought Scenarios")
"""

# compare RWU depth across scenarios

days, dates_out = get_dates(sim_ctr);

rwu_comp = DataFrame(date=dates_out);

rwu_comp.control = get_RWU_centroid(sim_ctr);
rwu_comp.irrigation = get_RWU_centroid(sim_irr);
rwu_comp.irrigation_vpd = get_RWU_centroid(sim_irr_vpd);
rwu_comp.roof = get_RWU_centroid(sim_drt);
rwu_comp.roof_vpd = get_RWU_centroid(sim_drt_vpd);

rwu_comp = stack(rwu_comp, Not(:date), variable_name="treatment", value_name="RWU");
rwu_comp.RWU = replace(rwu_comp.RWU, NaN=>missing);

R"""
rdf = $rwu_comp
ggplot(filter(rdf, date>="2024-01-01"), aes(x=date, y=RWU, color=treatment)) + geom_line() +
    facet_wrap(~treatment)+
    labs(title="RWU Depth Comparison across Scenarios")
"""

rwu_comp = dropmissing(rwu_comp);
rwu_comp.year = year.(rwu_comp.date);
rwu_yr = combine(groupby(rwu_comp, [:year, :treatment]), :RWU .=> [mean maximum]);

R"""
rdf = $rwu_yr
ggplot(filter(rdf, year>2003), aes(x=year, y=RWU_mean, group=treatment, fill=treatment)) + geom_col(position="dodge") +
    labs(title="Mean RWU Depth Comparison Across Scenarios")
"""