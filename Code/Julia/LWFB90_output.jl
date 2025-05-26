## run best performing scenario from LWFBrook90.jl calibration
## and examine output against observed data

using CSV, DataFrames, Dates, RCall;
using Plots; gr()
R"""
library(tidyverse)
"""

include("run_LWFB90_param.jl");

## behavioral data
obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
obs_swc = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
select!(obs_swc, :date, :depth, :VWC); # remove extra columns
filter!(:date => >(Date(2003, 4, 1)), obs_swc); # filter out early dates
obs_swc.src .= "obs"; # add source column

## retrieve metrics and parameters from calibration
met = CSV.read("LWFBcal_output/metrics_20250523.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_20250523.csv", DataFrame);

met_good = met[met.nse10 .> 0 .&& 
              met.nse40 .> 0 .&&
              met.nse60 .> 0 .&&
              met.nse80 .> 0, :];

# calculate quadratic mean of NSE values
met_good.nse_com = sqrt.((met_good.nse10 .^ 2 + met_good.nse40 .^ 2 + met_good.nse60 .^ 2 + met_good.nse80 .^ 2) / 4);

# find the index of the maximum value
max_idx = argmax(met_good.nse_com);
met_good[max_idx, :]
scen_max = met_good.scen[max_idx];

# parameter values for the best performing scenario
par_best = par[scen_max, :];
par_best

# run LWFBrook90.jl with single parameter set
sim = run_LWFB90_param(par_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/", "pfynwald", "LWFB_testrun/");
# after first time, run the model with new_folder=false
# sim = run_LWFB90_param(par_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/", "pfynwald", "LWFB_testrun/", new_folder=false);

## retrieve model output

# some plots
plotforcingandstates(sim)
plotamounts(sim, :above_and_belowground, :showRWUcentroid)
plotisotopes(sim, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)


days = range(sim.ODESolution.prob.tspan...);
dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

# soil water content
z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = days);
z.date = Date.(dates_out);
filter!(:date => >(obs_swc.date[1]), z);
select!(z, Not(:time));

# reshape data
z_long = stack(z, Not(:date), variable_name = "depth", value_name="VWC");
z_long.depth = parse.(Int, map(x -> x[(end-4):(end-2)], z_long.depth)) .÷ 10; # convert depth to cm from var name
z_long.src .= "sim"; # add source column

# combine observed and simulated data
swc_comp = [obs_swc; z_long];

# plot observed and simulated data using ggplot
R"""
rdf = $swc_comp
ggplot(rdf, aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Content Comparison")
"""

# plot just observed data
R"""
rdf = $obs_swc
ggplot(rdf, aes(x=date, y=VWC, color=as.factor(depth))) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Observed VWC Depth Comparison")
"""