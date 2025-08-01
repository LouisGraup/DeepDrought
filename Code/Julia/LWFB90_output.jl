## run best performing scenario from LWFBrook90.jl calibration
## and examine output against observed data

using CSV, DataFrames, Dates, RCall;
using Plots; gr()
R"""
library(tidyverse)
"""

include("run_LWFB90_param.jl");

# function to filter metrics for behavioral runs
function behavioral_met(met)
    return met[met.nse10 .> 0 .&& 
               met.nse40 .> 0 .&&
               met.nse60 .> 0 .&&
               met.nse80 .> 0, :]
end

# function to calculate combined NSE metric
function nse_quad(met)
    return sqrt.((met.nse10 .^ 2 + met.nse40 .^ 2 + met.nse60 .^ 2 + met.nse80 .^ 2) / 4);
end

## behavioral data
obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
filter!(:date => >(Date(2004, 1, 1)), obs_swc); # filter out early dates

# separate control and irrigation scenarios
obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
select!(obs_swc_ctr, :date, :depth, :VWC); # remove extra columns
obs_swc_ctr.src .= "obs"; # add source column

obs_swc_irr = obs_swc[obs_swc.meta .== "irrigated", :]; # select irrigation treatment
select!(obs_swc_irr, :date, :depth, :VWC); # remove extra columns
obs_swc_irr.src .= "obs"; # add source column

## retrieve metrics and parameters from calibration
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250529.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250529.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_20250529.csv", DataFrame);

met_ctr_good = behavioral_met(met_ctr);
met_irr_good = behavioral_met(met_irr);

# calculate quadratic mean of NSE values
met_ctr_good.nse_com = nse_quad(met_ctr_good);
met_irr_good.nse_com = nse_quad(met_irr_good);

# find the index of the maximum value
max_idx = argmax(met_ctr_good.nse_com);
met_ctr_good[max_idx, :]
scen_max = met_ctr_good.scen[max_idx];

# parameter values for the best performing scenario
par_best = par[scen_max, :];
par_best

# run LWFBrook90.jl with single parameter set
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2020, 12, 31), "LWFBinput/Pfyn_irrigation_ambient/", "pfynwald", "LWFB_testrun/irrigation/");
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
swc_comp = [obs_swc_ctr; z_long];

# plot observed and simulated data using ggplot
R"""
rdf = $swc_comp
ggplot(rdf, aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Content Comparison")
"""

# plot just observed data
R"""
rdf = $obs_swc_irr
ggplot(rdf, aes(x=date, y=VWC, color=as.factor(depth))) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Observed VWC Depth Comparison")
"""


# water balance modelling

using CairoMakie, AlgebraOfGraphics, DataFramesMeta, CategoricalArrays, Chain;

sim_ctr_wb = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/", new_folder=false, watbal=true);

# water partitioning
sim_water_part_daily, sim_water_part_monthly, wp_y = get_water_partitioning(sim_ctr_wb);
describe(sim_water_part_daily)

function plot_monthly_water_partitioning(df_partitioning_monthly)
    # color palette
    color_palette = reverse([
            #"Transpiration deficit" => :red2,
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
            #"Transpiration deficit" => :red2,
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

wb_mth = plot_monthly_water_partitioning(sim_water_part_monthly);

wb_yr = plot_yearly_water_partitioning(wp_y);

color_palette_Schmidt_Walter2020 = reverse([
            "Td"            => colorant"#d01c8b", #:palevioletred4,
            "Ta"            => colorant"#abdda4", #:darkseagreen,
            "Einterception" => colorant"#fdae61", #:lightyellow,
            "Esoil"         => colorant"#FFD700", #:navajowhite, #:orange2,
            "Esnow"         => :white,
            "R"             => colorant"#91bfdb", #:slategray2,
            "D"             => colorant"#2b83ba", #:skyblue4,
            "Precip"        => :black,
            # "ETa"  => :black,
            # "P2"   => :darkblue,
            # "Swat" => :brown
        ])

color_palette_KM = reverse([
            "Td"            => :red2, # transpiration deficit
            "Ta"            => :darkolivegreen2, # actual transpiration
            "Einterception" => :forestgreen, # evaporation from interception
            "Esoil"         => :khaki3, # soil evaporation
            #"Esnow"         => :white, # snow sublimation
            "R"             => :lightskyblue, # saturation excess runoff + downslope flow + bypass flow
            "D"             => :steelblue4, # vertical flow (drainage)
            # "ETa"   => :black, # actual evapotranspiration 
            # "P2"    => :darkblue,
            "Precip"     => :darkblue,
            # "Swat" => :brown # storage
        ]);