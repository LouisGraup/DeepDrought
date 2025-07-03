
using LWFBrook90
using CSV, DataFrames, Dates
using Distributions

cd("DeepDrought/Code/Julia/") # for cluster

include("run_LWFB90_param.jl");

obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
filter!(:date => >=(Date(2004, 1, 1)), obs_swc); # filter out early dates

obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
select!(obs_swc_ctr, :date, :depth, :VWC); # remove extra columns
obs_swc_ctr = unstack(obs_swc_ctr, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
sort!(obs_swc_ctr, :date); # sort by date

# objective function to compare model output to observed data
function obj_fun_swc(sim, obs)

    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :VWC_10cm]]);
    obs_40cm = dropmissing(obs[!, [:date, :VWC_40cm]]);
    obs_60cm = dropmissing(obs[!, [:date, :VWC_60cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :VWC_80cm]]);

    filter!(:date => >(Date(2015, 1, 1)), obs_80cm); # filter out early dates

    # match simulated data to available dates for each depth
    sim_10cm = sim[sim.dates .∈ [obs_10cm.date], :theta_m3m3_100mm];
    sim_40cm = sim[sim.dates .∈ [obs_40cm.date], :theta_m3m3_400mm];
    sim_60cm = sim[sim.dates .∈ [obs_60cm.date], :theta_m3m3_600mm];
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
    nse40 = NSE(sim_40cm, obs_40cm.VWC_40cm);
    nse60 = NSE(sim_60cm, obs_60cm.VWC_60cm);
    nse80 = NSE(sim_80cm, obs_80cm.VWC_80cm);

    # calculate RMSE
    rmse10 = RMSE(sim_10cm, obs_10cm.VWC_10cm);
    rmse40 = RMSE(sim_40cm, obs_40cm.VWC_40cm);
    rmse60 = RMSE(sim_60cm, obs_60cm.VWC_60cm);
    rmse80 = RMSE(sim_80cm, obs_80cm.VWC_80cm);

    return nse10, rmse10, nse40, rmse40, nse60, rmse60, nse80, rmse80
end


param = CSV.read("LWFBcal_output/param.csv", DataFrame);
par = param[1, :];

sim = run_LWFB90_param(par, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFBcalibration/cal_ctr1", new_folder=false);


days = range(sim.ODESolution.prob.tspan...);
dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = days);
z.dates = Date.(dates_out);

met = obj_fun_swc(z, obs_swc_ctr);

println(met)