## run best performing scenario from LWFBrook90.jl calibration
## for all Pfynwald scenarios
## and examine output against observed data

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions, RCall;
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Measures, Plots; gr()

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

    filter!(:date => >=(obs_swc.date[1]), z_theta);

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

    filter!(:date => >=(obs_swp.date[1]), z_psi);

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

    filter!(:date => >=(obs_sap.date[1]), z_trans);
    filter!(:date => <(Date(2023, 1, 1)), z_trans);

    obs_sap = select(obs_sap, Not(:meta));

    sap_comp = sort(leftjoin(z_trans, obs_sap, on = :date), :date);

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
    
    return col_RWU_centroid_mm, RWU_percent
end

function get_eff_swp(sim)

    days, dates_out = get_dates(sim);
    swp = DataFrame(date = dates_out);
    
    swp.RWU, rwu_per = get_RWU_centroid(sim); # RWU depth and percent

    swp_all = get_soil_(:psi, sim, days_to_read_out_d=days); # swp

    swp.swp_eff .= sum(rwu_per .* Matrix(swp_all[:, Not(:time)])', dims=1)';

    return swp
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
#filter!(:date => <(Date(2024, 1, 1)), obs_swc); # filter out late dates

# soil water potential
obs_swp = CSV.read("../../Data/Pfyn/PFY_swpc_corr.csv", DataFrame);
#filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
filter!(:date => <(Date(2024, 1, 1)), obs_swp); # filter out late dates

# sap flow
obs_sap = CSV.read("../../Data/Pfyn/PFY_sap.csv", DataFrame);

# separate control and irrigation scenarios
obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
obs_swc_irr = obs_swc[obs_swc.meta .== "irrigation", :]; # select irrigation treatment
obs_swc_irst = obs_swc[obs_swc.meta .== "irrigation_stop", :]; # select irrigation stop treatment

obs_swp_ctr = obs_swp[obs_swp.meta .== "control", :]; # select control treatment
obs_swp_irr = obs_swp[obs_swp.meta .== "irrigation", :]; # select irrigation treatment
obs_swp_irst = obs_swp[obs_swp.meta .== "irrigation_stop", :]; # select irrigation stop treatment

obs_sap_ctr = obs_sap[obs_sap.meta .== "Control", :]; # select control treatment
obs_sap_irr = obs_sap[obs_sap.meta .== "Irrigation", :]; # select irrigation treatment
obs_sap_irst = obs_sap[obs_sap.meta .== "Irrigation Stop", :]; # select irrigation stop treatment

# irrigation dates
on = ["2003-06-19", "2004-05-15", "2005-04-23", "2006-05-06", "2007-05-04", "2008-05-15", 
       "2009-05-14", "2010-06-17", "2011-05-14", "2012-05-11", "2013-05-17", "2014-05-19",
       "2015-05-12", "2016-05-30", "2017-05-16", "2018-05-08", "2019-05-17", "2020-05-26",
       "2021-05-18", "2022-05-11", "2023-05-15", "2024-05-26", "2025-04-30"];

off = ["2003-10-21", "2004-10-26", "2005-10-04", "2006-10-25", "2007-10-02", "2008-10-14",
        "2009-10-12", "2010-10-01", "2011-10-16", "2012-10-02", "2013-09-23", "2014-10-01",
        "2015-10-05", "2016-09-26", "2017-10-09", "2018-09-27", "2019-10-15", "2020-10-19",
        "2021-10-11", "2022-10-21", "2023-10-18", "2024-10-24", "2025-07-31"];

irr = DataFrame(on=Date.(on), off=Date.(off));
irr.year = year.(irr.on);

# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2023, 12, 31), "LWFBinput/Pfyn_irrigiso_ambient/", "pfynwald", "LWFB_testrun/irrigation/"; irrig=true);
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
draw(data(swp_comp_irst)*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Potential Comparison for Control Scenario", titlealign = :center)
)

draw(data(swp_comp_ctr)*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; label="Source")), facet = (; linkyaxes = :none),
    figure = (; size=(1200, 800), title="Soil Water Potential Comparison for Control Scenario", titlealign = :center)
)

draw(data(swp_comp_irr)*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Potential Comparison for Control Scenario", titlealign = :center)
)

# soil water content

draw(data(swc_comp_ctr)*
    mapping(:date, :VWC, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Color = (; label="Source")),
    figure = (; size=(1200, 800), title="Soil Water Content Comparison for Control Scenario", titlealign = :center)
)

draw(data(swc_comp_irr)*
    mapping(:date, :VWC, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Content Comparison for Control Scenario", titlealign = :center)
)

draw(data(swc_comp_irst)*
    mapping(:date, :VWC, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Content Comparison for Control Scenario", titlealign = :center)
)

# sap flow

draw(data(stack(sap_comp_ctr, Not(:date), variable_name=:Source))*
    mapping(:date, :value, color=:Source)*visual(Scatter, markersize=6),
    scales(X = (; label="")),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Control Scenario", titlealign = :center)
)

draw(data(dropmissing(sap_comp_ctr))*
    mapping(:trans, :sfd)*(visual(Scatter)+linear()),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Control Scenario", titlealign = :center)
)

sap_comp_irr.trans_sm = runmean(sap_comp_irr.trans, 14);
sap_comp_irr.sfd_sm = runmean(sap_comp_irr.sfd, 14);

sap_comp_irr_long = stack(sap_comp_irr[:, Not([:trans_sm, :sfd_sm])], Not(:date), variable_name=:Source);
sap_comp_irr_long2 = stack(sap_comp_irr[:, Not([:trans, :sfd])], Not(:date), value_name=:value_sm, variable_name=:Source);

sap_comp_irr_long.value_sm = sap_comp_irr_long2.value_sm;

draw(data(sap_comp_irr_long)*
    (mapping(:date, :value, color=:Source)*visual(Scatter, markersize=6)+
    mapping(:date, :value_sm, color=:Source)*visual(Lines)),
    scales(X = (; label="")),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Irrigation Scenario", titlealign = :center)
)

draw(data(dropmissing(sap_comp_irr))*
    mapping(:trans, :sfd)*(visual(Scatter)+linear()),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Irrigation Scenario", titlealign = :center)
)

draw(data(stack(sap_comp_irst, Not(:date), variable_name=:Source))*
    mapping(:date, :value, color=:Source)*visual(Scatter, markersize=6),
    scales(X = (; label="")),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Irrigation Stop Scenario", titlealign = :center)
)

draw(data(dropmissing(sap_comp_irst))*
    mapping(:trans, :sfd)*(visual(Scatter)+linear()),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Irrigation Stop Scenario", titlealign = :center)
)

# compare transpiration across scenarios
ctr_sap = get_sap(sim_ctr);
ctr_sap.trans_sm = runmean(ctr_sap.trans, 14);
ctr_sap.scen .= "Control";

irr_sap = get_sap(sim_irr);
irr_sap.trans_sm = runmean(irr_sap.trans, 14);
irr_sap.scen .= "Irrigation";

irst_sap = get_sap(sim_irst);
irst_sap.trans_sm = runmean(irst_sap.trans, 14);
irst_sap.scen .= "Irrigation Stop";

sap_comp_scen = [ctr_sap; irr_sap; irst_sap];

# format x-axis ticks
ticks = ctr_sap[month.(ctr_sap.date) .== 1 .&& day.(ctr_sap.date) .== 1 .&& year.(ctr_sap.date) .< 2021, :date];
aogticks = datetimeticks(ticks, string.(ticks));
draw(
    data(irr[irr.year .>= 2014 .&& irr.year .< 2020, :])*
    mapping(:on, :off)*visual(VSpan, color=:lightblue, alpha=0.4)+
    data(sap_comp_scen[sap_comp_scen.date.<Date("2020-01-01"), :])*
    (mapping(:date, :trans, color=:scen)*visual(Scatter, markersize=6)+
    mapping(:date, :trans_sm, color=:scen)*visual(Lines, linewidth=3)),
    scales(Color = (; label="Scenario", palette = ["#E69F00","#56B4E9","#009E73"]),
           X = (; label=""), Y= (; label="Transpiration (mm/day)")),
    figure = (; size=(1200, 600)), axis = (; xticks=aogticks, title="Modelled Transpiration Comparison across Scenarios", titlesize=20)
)

# compare soil water potential across scenarios
ctr_swp = get_swp(sim_ctr);
ctr_swp.scen .= "Control";

irr_swp = get_swp(sim_irr);
irr_swp.scen .= "Irrigation";

irst_swp = get_swp(sim_irst);
irst_swp.scen .= "Irrigation Stop";

swp_comp_scen = [ctr_swp; irr_swp; irst_swp];

draw(data(swp_comp_scen)*
    mapping(:date, :SWP, color=:scen, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; palette = ["#E69F00","#56B4E9","#009E73"], label="Scenario")),
    figure = (; size=(800, 600), title="Modelled Soil Water Potential Comparison across Scenarios", titlealign = :center)
)


# compare RWU depth across scenarios

days, dates_out = get_dates(sim_ctr);

ctr_rwu = DataFrame(date=dates_out);
ctr_rwu.RWU, ctr_rwu_per = get_RWU_centroid(sim_ctr);
ctr_rwu.RWU_sm = runmean(ctr_rwu.RWU, 14);
ctr_rwu.scen .= "Control";

irr_rwu = DataFrame(date=dates_out);
irr_rwu.RWU, = get_RWU_centroid(sim_irr);
irr_rwu.RWU_sm = runmean(irr_rwu.RWU, 14);
irr_rwu.scen .= "Irrigation";

irst_rwu = DataFrame(date=dates_out);
irst_rwu.RWU, = get_RWU_centroid(sim_irst);
irst_rwu.RWU_sm = runmean(irst_rwu.RWU, 14);
irst_rwu.scen .= "Irrigation Stop";

rwu_comp = [ctr_rwu; irr_rwu; irst_rwu];

rwu_comp.RWU = replace(rwu_comp.RWU, NaN=>missing);
rwu_comp.RWU_sm = replace(rwu_comp.RWU_sm, NaN=>missing);

rwu_comp.month = month.(rwu_comp.date);
rwu_comp.year = year.(rwu_comp.date);
rwu_comp_med = dropmissing(rwu_comp[rwu_comp.month .> 3 .&& rwu_comp.month .< 12 .&& rwu_comp.year .> 2002, Not(:RWU_sm)]);
med_rwu_comp = combine(groupby(rwu_comp_med, :scen), :RWU .=> [median mean]);

draw(
    data(irr[irr.year .>= 2014 .&& irr.year .< 2020, :])*
    mapping(:on, :off)*visual(VSpan, color=:lightblue, alpha=0.4)+
    data(rwu_comp[rwu_comp.date.<Date("2020-01-01"), :])*
    (mapping(:date, :RWU, color=:scen)*visual(Scatter, markersize=6)+
    mapping(:date, :RWU_sm, color=:scen)*visual(Lines, linewidth=3)),
    scales(Color = (; label="Scenario", palette = ["#E69F00","#56B4E9","#009E73"]),
           X = (; label=""), Y= (; label="Root Water Uptake Depth (mm)")),
    figure = (; size=(1200, 600)), axis = (; xticks=aogticks, title="RWU Depth Comparison across Scenarios", titlesize=20)
)

# compare annual transpiration across scenarios
ctr_yr = ann_trans(sim_ctr);
ctr_yr.scen .= "Control";

irr_yr = ann_trans(sim_irr);
irr_yr.scen .= "Irrigation";

irst_yr = ann_trans(sim_irst);
irst_yr.scen .= "Irrigation Stop";

yr_comp = [ctr_yr; irr_yr; irst_yr];

draw(data(yr_comp[yr_comp.year .>= 2003, :])*
    mapping(:year, :trans_sum, color=:scen, dodge=:scen)*visual(BarPlot),
    scales(Color = (; label="Scenario", palette = ["#E69F00","#56B4E9","#009E73"]),
           X = (; label="Year"), Y= (; label="Annual Transpiration (mm)")),
    axis = (; title="Annual Transpiration Comparison Across Scenarios")
)


# compare RWU against transpiration

ctr_rwu_tran = get_sap(sim_ctr);
ctr_rwu_tran.RWU, = get_RWU_centroid(sim_ctr);
ctr_rwu_tran.scen .= "Control";

irr_rwu_tran = get_sap(sim_irr);
irr_rwu_tran.RWU, = get_RWU_centroid(sim_irr);
irr_rwu_tran.scen .= "Irrigation";

df_rwu_tran = [ctr_rwu_tran; irr_rwu_tran];
df_rwu_tran.RWU = replace(df_rwu_tran.RWU, NaN=>missing);
df_rwu_tran.month = month.(df_rwu_tran.date);
df_rwu_tran.year = year.(df_rwu_tran.date);
df_rwu_tran = df_rwu_tran[df_rwu_tran.year .> 2002, :];

draw(
    data(df_rwu_tran[df_rwu_tran.month .> 3 .&& df_rwu_tran.month .< 12, :])*
    mapping(:trans, :RWU, color=:month => nonnumeric, layout=:scen)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp[med_rwu_comp.scen .!= "Irrigation Stop",:])*mapping(:RWU_mean, layout=:scen)*visual(HLines, color=:black, linestyle=:dash),
    scales(Color = (; label="Month", palette = from_continuous(:seaborn_bright6)),
           X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), facet = (; linkxaxes = :none), figure = (; size=(800, 400))
)

# compare rwu against swp in each layer

days, dates_out = get_dates(sim_ctr);

# control scenario

ctr_rwu_all = get_soil_(:RWU, sim_ctr, days_to_read_out_d=days);
ctr_rwu_all.date = dates_out;
ctr_rwu_all = stack(ctr_rwu_all[:,Not(:time)], Not([:date]), value_name=:RWU);
ctr_rwu_all.depth = parse.(Int, [match(r"[0-9]+", s).match for s in ctr_rwu_all.variable]);
select!(ctr_rwu_all, Not(:variable));

ctr_swp_all = get_soil_(:psi, sim_ctr, days_to_read_out_d=days);
ctr_swp_all.date = dates_out;
ctr_swp_all = stack(ctr_swp_all[:,Not(:time)], Not([:date]), value_name=:SWP);
ctr_swp_all.depth = parse.(Int, [match(r"[0-9]+", s).match for s in ctr_swp_all.variable]);
select!(ctr_swp_all, Not(:variable));

ctr_rwu_swp = leftjoin(ctr_rwu_all, ctr_swp_all, on=[:date, :depth]);
depth_bins = [0, 70, 100, 200, 250, 400, 500, 600, 800, 1000, 1200, 1500, 1900, 2000];
ctr_rwu_swp.depth_bin = cut(ctr_rwu_swp.depth .- 1, depth_bins, extend=true);

ctr_rwu_swp = ctr_rwu_swp[year.(ctr_rwu_swp.date) .> 2002, :];
ctr_rwu_swp.scenario .= "Control";

draw(data(ctr_rwu_swp[ctr_rwu_swp.RWU .> 0, :])*
    mapping(:SWP, :RWU, color=:depth_bin)*visual(Scatter, alpha=.5, markersize=6),
    scales(Color = (; label="Uptake Depth (mm)", palette = from_continuous(:viridis)),
           X = (; label="Soil Water Potential (kPa)"), Y= (; label="Root Water Uptake (mm/day)")),
    axis = (; title="Daily Root Water Uptake by Depth and Soil Water Potential", titlesize=20),
    figure = (; size=(800, 600))
)

# daily average RWU
ctr_rwu_swp.month = month.(ctr_rwu_swp.date);
ctr_rwu_swp.day = day.(ctr_rwu_swp.date);

ctr_rwu_swp_daily = combine(groupby(ctr_rwu_swp[ctr_rwu_swp.depth .!= 2000, :], [:month, :day, :depth_bin]), :RWU .=> mean);
ctr_rwu_swp_daily.date = Date.(0000, ctr_rwu_swp_daily.month, ctr_rwu_swp_daily.day); # dummy year for plotting

ctr_rwu_swp_wide = unstack(ctr_rwu_swp_daily[:, Not([:month, :day])], :depth_bin, :RWU_mean);
ctr_rwu_swp_run = hcat(map(x -> runmean(x, 7), eachcol(ctr_rwu_swp_wide[:, Not(:date)]))...);

tickdates = Date.(["0000-01-01", "0000-04-01", "0000-07-01", "0000-10-01", "0001-01-01"]);
areaplot(ctr_rwu_swp_wide.date, ctr_rwu_swp_run, label=reshape(unique(ctr_rwu_swp_daily.depth_bin), 1, 12),
    xticks=(tickdates, Dates.format.(tickdates, "mm-dd")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="Day of Year", ylabel="Root Water Uptake (mm/day)", legendtitle="RWU Depth (mm)",
    title="\nMean Daily Root Water Uptake by Depth for Control Scenario", size=(1200, 800))

# irrigation scenario

irr_rwu_all = get_soil_(:RWU, sim_irr, days_to_read_out_d=days);
irr_rwu_all.date = dates_out;
irr_rwu_all = stack(irr_rwu_all[:,Not(:time)], Not([:date]), value_name=:RWU);
irr_rwu_all.depth = parse.(Int, [match(r"[0-9]+", s).match for s in irr_rwu_all.variable]);
select!(irr_rwu_all, Not(:variable));

irr_swp_all = get_soil_(:psi, sim_irr, days_to_read_out_d=days);
irr_swp_all.date = dates_out;
irr_swp_all = stack(irr_swp_all[:,Not(:time)], Not([:date]), value_name=:SWP);
irr_swp_all.depth = parse.(Int, [match(r"[0-9]+", s).match for s in irr_swp_all.variable]);
select!(irr_swp_all, Not(:variable));

irr_rwu_swp = leftjoin(irr_rwu_all, irr_swp_all, on=[:date, :depth]);
depth_bins = [0, 70, 100, 200, 250, 400, 500, 600, 800, 1000, 1200, 1500, 1900, 2000];
irr_rwu_swp.depth_bin = cut(irr_rwu_swp.depth .- 1, depth_bins, extend=true);

irr_rwu_swp = irr_rwu_swp[year.(irr_rwu_swp.date) .> 2002, :];
irr_rwu_swp.scenario .= "Irrigation";

# daily average RWU
irr_rwu_swp.month = month.(irr_rwu_swp.date);
irr_rwu_swp.day = day.(irr_rwu_swp.date);

irr_rwu_swp_daily = combine(groupby(irr_rwu_swp[irr_rwu_swp.depth .!= 2000, :], [:month, :day, :depth_bin]), :RWU .=> mean);
irr_rwu_swp_daily.date = Date.(0000, irr_rwu_swp_daily.month, irr_rwu_swp_daily.day); # dummy year for plotting

irr_rwu_swp_wide = unstack(irr_rwu_swp_daily[:, Not([:month, :day])], :depth_bin, :RWU_mean);
irr_rwu_swp_run = hcat(map(x -> runmean(x, 7), eachcol(irr_rwu_swp_wide[:, Not(:date)]))...);

areaplot(irr_rwu_swp_wide.date, irr_rwu_swp_run, label=reshape(unique(irr_rwu_swp_daily.depth_bin), 1, 12),
    xticks=(tickdates, Dates.format.(tickdates, "mm-dd")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="Day of Year", ylabel="Root Water Uptake (mm/day)", legendtitle="RWU Depth (mm)",
    title="\nMean Daily Root Water Uptake by Depth for Irrigation Scenario", size=(1200, 800))

# scenario comparison

comp_rwu_swp = [ctr_rwu_swp; irr_rwu_swp];

draw(data(comp_rwu_swp[comp_rwu_swp.RWU .> 0, :])*
    mapping(:SWP, :RWU, color=:depth_bin, layout=:scenario)*visual(Scatter, alpha=.5, markersize=6),
    scales(Color = (; label="Uptake Depth (mm)", palette = from_continuous(:viridis)),
           X = (; label="Soil Water Potential (kPa)"), Y= (; label="Root Water Uptake (mm/day)")), facet = (; linkxaxes = :none),
    figure = (; size=(1600, 800), title="Daily Root Water Uptake by Depth and Soil Water Potential", titlesize=20, titlealign = :center)
)

# weighted-average soil water potential

# control
ctr_rwu_swp_med = get_eff_swp(sim_ctr);
ctr_rwu_swp_med = leftjoin(ctr_rwu_swp_med, get_sap(sim_ctr), on=:date);
ctr_rwu_swp_med = ctr_rwu_swp_med[ctr_rwu_swp_med.date .>= Date(2003, 1, 1), :];
ctr_rwu_swp_med.month = month.(ctr_rwu_swp_med.date);
ctr_rwu_swp_med.scenario .= "Control";

draw(data(ctr_rwu_swp_med)*mapping(:date, :swp_eff)*visual(Lines))

draw(data(ctr_rwu_swp_med)*
    mapping(:swp_eff, :trans, color=:RWU)*visual(Scatter, alpha=.5, markersize=6),
    scales(X = (; label="Weighted-Average Soil Water Potential (kPa)"), Y= (; label="Transpiration (mm/day)")),
    axis = (; title="Transpiration vs. Effective SWP", titlesize=16)
)

draw(data(ctr_rwu_swp_med)*
    mapping(:swp_eff, :RWU, color=:month => nonnumeric)*visual(Scatter, alpha=.5, markersize=6),
    scales(Color = (; palette = from_continuous(:seaborn_colorblind6)),
        X = (; label="Weighted-Average Soil Water Potential (kPa)"), Y= (; label="Weighted-Average Root Water Uptake Depth (mm)")),
    axis = (; yreversed = true, title="Weighted-Average RWU Depth vs. SWP", titlesize=16)
)

# irrigation
irr_rwu_swp_med = get_eff_swp(sim_irr);
irr_rwu_swp_med = leftjoin(irr_rwu_swp_med, get_sap(sim_irr), on=:date);
irr_rwu_swp_med = irr_rwu_swp_med[irr_rwu_swp_med.date .>= Date(2003, 1, 1), :];
irr_rwu_swp_med.month = month.(irr_rwu_swp_med.date);
irr_rwu_swp_med.scenario .= "Irrigation";

comp_rwu_swp_eff = [ctr_rwu_swp_med; irr_rwu_swp_med];

draw(data(comp_rwu_swp_eff)*mapping(:date, :swp_eff, color=:scenario)*visual(Lines))

draw(data(comp_rwu_swp_eff)*
    mapping(:swp_eff, :trans, color=:RWU, layout=:scenario)*visual(Scatter, alpha=.5, markersize=6),
    scales(X = (; label="Weighted-Average Soil Water Potential (kPa)"), Y= (; label="Transpiration (mm/day)"), Color = (; label="RWU Depth (mm)")),
    figure = (; size=(1200, 600), title="Transpiration vs. Effective SWP", titlesize=16, titlealign = :center)
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