# sensitivity analysis of irrigation stop scenario
# comparing LAI trajectories

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions;
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Measures, StatsPlots, Plots; gr()

include("run_LWFB90_param.jl");

## load behavioral data and define objective function

# behavioral data
# soil water content and soil matric potential
obs_swc = CSV.read("../../Data/Pfyn/Pfyn_swat.csv", DataFrame);
filter!(:date => >=(Date(2016, 1, 1)), obs_swc); # filter out early dates
filter!(:date => <(Date(2021, 1, 1)), obs_swc); # filter out late dates

obs_swp = CSV.read("../../Data/Pfyn/Pfyn_swp.csv", DataFrame);
filter!(:date => >=(Date(2016, 1, 1)), obs_swp); # filter out early dates
filter!(:date => <(Date(2021, 1, 1)), obs_swp); # filter out late dates

# up-scaled sap flow data
obs_sap = CSV.read("../../Data/Pfyn/Pfyn_trans_2011_17.csv", DataFrame);

# separate out irrigation scenario
obs_swc_irst = obs_swc[obs_swc.meta .== "stop", :]; # select irrigation stop treatment
filter!(:date => >=(Date(2016, 1, 1)), obs_swc_irst); # filter out early dates
select!(obs_swc_irst, :date, :depth, :VWC); # remove extra columns
obs_swc_irst = unstack(obs_swc_irst, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
sort!(obs_swc_irst, :date); # sort by date

obs_swp_irst = obs_swp[obs_swp.meta .== "stop", :]; # select irrigation stop treatment
select!(obs_swp_irst, :date, :depth, :SWP); # remove extra columns
obs_swp_irst = unstack(obs_swp_irst, :date, :depth, :SWP, renamecols=x->Symbol("SWP_$(x)cm")); # reshape data
sort!(obs_swp_irst, :date); # sort by date

obs_sap_irst = obs_sap[obs_sap.scen .== "Irrigation stop", :]; # select irrigation stop treatment

# objective function to compare model output to observed data
function obj_fun_swc(sim, obs)

    obs = obs[obs.date .∈ [sim.dates], :];

    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :VWC_10cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :VWC_80cm]]);

    # match simulated data to available dates for each depth
    sim_10cm = sim[sim.dates .∈ [obs_10cm.date], :theta_m3m3_100mm];
    sim_80cm = sim[sim.dates .∈ [obs_80cm.date], :theta_m3m3_800mm];

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    function RMSE(sim, obs)
        # calculate Root Mean Square Error
        rmse = sqrt(sum((obs .- sim).^2) / length(obs))
        return rmse
    end

    # calculate NSE
    nse10 = NSE(sim_10cm, obs_10cm.VWC_10cm);
    nse80 = NSE(sim_80cm, obs_80cm.VWC_80cm);

    # calculate RMSE
    rmse10 = RMSE(sim_10cm, obs_10cm.VWC_10cm);
    rmse80 = RMSE(sim_80cm, obs_80cm.VWC_80cm);

    return nse10, rmse10, nse80, rmse80
end

function obj_fun_swp(sim, obs)

    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :SWP_10cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :SWP_80cm]]);

    # match simulated data to available dates for each depth
    sim_10cm = sim[sim.dates .∈ [obs_10cm.date], :psi_kPa_100mm];
    sim_80cm = sim[sim.dates .∈ [obs_80cm.date], :psi_kPa_800mm];

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

function obs_fun_trans(sap_comp)

    # remove missing values
    sap_comp = dropmissing(sap_comp);

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    nse = NSE(sap_comp.trans, sap_comp.Tr);

    return nse

end

function obs_fun_sap(sap_comp)

    # remove missing values
    sap_comp = dropmissing(sap_comp);

    cc = cor(sap_comp.trans, sap_comp.Tr);

    return cc

end

function sap_combine(sim, obs_sap)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_trans = get_sap(sim);

    z_trans = z_trans[z_trans.date .∈ [obs_sap.date], :];
    
    sap_comp = leftjoin(z_trans, obs_sap[:, Not(:scen)], on = :date);

    return sap_comp

end

function get_sap(sim)
    # retrieve transpiration from sim
    
    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

function get_output(sim)
    
    function get_dates(sim)
        days = range(sim.ODESolution.prob.tspan...)[Not(end)];
        dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
        return days, Date.(dates_out)
    end
    
    days, dates_to_read_out = get_dates(sim);

    # soil water content
    z_theta = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = days);
    z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = days);

    # add dates column
    z_theta.dates = dates_to_read_out;
    z_psi.dates = dates_to_read_out;

    swc_met = obj_fun_swc(z_theta, obs_swc_irst)
    swp_met = obj_fun_swp(z_psi, obs_swp_irst)

    # transpiration
    sap_comp = sap_combine(sim, obs_sap_irst);
    sap_cor = obs_fun_sap(sap_comp);
    sap_nse = obs_fun_trans(sap_comp);

    z_trans = get_sap(sim);
    max_trans = maximum(z_trans.trans);

    # annual total
    z_trans.year = year.(z_trans.date);
    z_trans_yr = combine(groupby(z_trans, :year), :trans => sum);
    mean_ann_trans = mean(z_trans_yr.trans_sum);

    tr_met = [sap_cor, sap_nse, max_trans, mean_ann_trans];

    return (swc_met, swp_met, tr_met)
end


### BEGIN USER INPUT ###

## parameter inputs
input_prefix = "pfynwald";

# simulation dates
start_date = Date(2010, 1, 1);
end_date = Date(2024, 12, 31);

# use single parameter set for analysis
par_irst = CSV.read("LWFBcal_output/Pfyn_irst_param_best.csv", DataFrame);
scen_best = 1505;
par_irst_best = par_irst[par_irst.scen .== scen_best, Not(:scen)];

# default run
sim_irst_def = run_LWFB90_param(par_irst_best, start_date, end_date, "LWFBinput/Pfyn_irrigiso_stop/", input_prefix, "LWFB_LAIsens/irr_stop_def/", irrig=true);
out_irst_def = get_output(sim_irst_def);

# positive legacy effect
sim_irst_pos = run_LWFB90_param(par_irst_best, start_date, end_date, "LWFB_LAIsens/Pfyn_irrigiso_stop_pos/", input_prefix, "LWFB_LAIsens/irr_stop_pos/", irrig=true);
out_irst_pos = get_output(sim_irst_pos);

# negative legacy effect
sim_irst_neg = run_LWFB90_param(par_irst_best, start_date, end_date, "LWFB_LAIsens/Pfyn_irrigiso_stop_neg/", input_prefix, "LWFB_LAIsens/irr_stop_neg/", irrig=true);
out_irst_neg = get_output(sim_irst_neg);


# intialize metric dataframes
metrics = DataFrame(scen = Any[],
    swc_nse10 = Float64[],
    swc_rmse10 = Float64[],
    swc_nse80 = Float64[],
    swc_rmse80 = Float64[],
    swp_nse10 = Float64[],
    swp_nse80 = Float64[],
    trans_cor = Float64[], 
    trans_nse = Float64[], 
    max_trans = Float64[],
    ann_trans = Float64[]);

# loop through results

# default run
swc_def, swp_def, tr_def = out_irst_def;
row_def = ["default", swc_def..., swp_def..., tr_def...];
push!(metrics, row_def);

# positive run
swc_pos, swp_pos, tr_pos = out_irst_pos;
row_pos = ["positive", swc_pos..., swp_pos..., tr_pos...];
push!(metrics, row_pos);

# negative run
swc_neg, swp_neg, tr_neg = out_irst_neg;
row_neg = ["negative", swc_neg..., swp_neg..., tr_neg...];
push!(metrics, row_neg);

CSV.write("LWFB_LAIsens/metrics_irst_laisens.csv", metrics);
