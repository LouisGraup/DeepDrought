## run best performing scenario from LWFBrook90.jl calibration
## for all Pfynwald scenarios
## and examine output against observed data

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions;
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
    days, dates_out = get_dates(sim);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = days);
    z.date = dates_out;

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
    days, dates_out = get_dates(sim);
    
    z = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = days);
    z.date = dates_out;

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
    days, dates_out = get_dates(sim);
    
    z = get_fluxes(sim);
    z.date = dates_out;
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

function get_soil_iso(sim; shape = "long")
    # retrieve soil water potential data from sim
    days, dates_out = get_dates(sim);

    z = get_soil_([:d18O, :d2H], sim, depths_to_read_out_mm = [50, 200, 400], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "var", value_name="value");
        z_long.iso = map(x -> occursin("d18O", x) ? "d18O" : "d2H", z_long.var);
        z_long.depth = parse.(Int, replace.(z_long.var, r".*_([0-9]+)mm$" => s"\1")) .÷ 10; # convert depth to cm from var name
        select!(z_long, Not(:var));
        z_long = unstack(z_long, [:date, :depth], :iso, :value);
        return z_long
    else
        return z
    end
    
end

function get_xy_iso(sim)
    # retrieve xylem isotope data from sim
    days, dates_out = get_dates(sim);

    z = get_states(sim, days_to_read_out_d = days);
    z.date = Date.(dates_out);
    select!(z, :date, :XYLEM_d18O, :XYLEM_d2H);
    
    return z
end

function get_trans_def(sim)
    days, dates_out = get_dates(sim);

    z = get_fluxes(sim);
    z.date = dates_out;
    z.trans = z.cum_d_tran;
    z.pet = z.cum_d_ptran;
    z.td = z.pet .- z.trans;

    select!(z, :date, :trans, :pet, :td);

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
    # obs is the observed sap flow data

    z_trans = get_sap(sim);

    filter!(:date => >=(obs_sap.date[1]), z_trans);
    filter!(:date => <(Date(2023, 1, 1)), z_trans);

    obs_sap = select(obs_sap, Not(:scen));

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

function get_REW(sim)
    # derive relate extractable soil water based on maximum root depth
    # for now, just use normalized SWAT over all layers

    days, dates_out = get_dates(sim);
    
    #swat = get_states(sim);

    #REW = (swat.SWAT_mm .- minimum(swat.SWAT_mm)) ./ (maximum(swat.SWAT_mm) - minimum(swat.SWAT_mm));

    #return REW

    # get # root zone layers from maximum rooting depth
    max_root_depth = sim.parametrizedSPAC.pars.root_distribution.z_rootMax_m;

    upper_lim = sim.parametrizedSPAC.soil_discretization.df.Upper_m;
    lower_lim = sim.parametrizedSPAC.soil_discretization.df.Lower_m;

    rz_layers = lower_lim[upper_lim .>= max_root_depth];
    nrz = length(rz_layers);

    # soil water content
    swc = get_soil_(:theta, sim, days_to_read_out_d=days); # soil water content
    select!(swc, Not(:time));
    swc = swc[!, 1:nrz]; # only keep soil water content for root zone layers

    # get field capacity and wilting point for available water content
    thf = sim.ODESolution.prob.p.p_soil.p_THETAF[1:nrz]; # field capacity soil water content

    # use swc to derive minimum soil water content for each layer
    # as approximation of wilting point
    # (more correct would be to derive swc from psi_crit using PTF)
    swc_min = minimum.(eachcol(swc));

    awc = thf .- swc_min; # available water capacity per layer
    
    REW_df = deepcopy(swc);
    for i in 1:nrz
        REW_df[!, i] = (swc[!, i] .- swc_min[i]) ./ awc[i];
    end

    REW = mean.(eachrow(REW_df));

    return REW
end

function get_eff_swp(sim)
    # derive effective soil water potential based on RWU
    days, dates_out = get_dates(sim);
    swp = DataFrame(date = dates_out);
    
    swp.RWU, rwu_per = get_RWU_centroid(sim); # RWU depth and percent

    swp_all = get_soil_(:psi, sim, days_to_read_out_d=days); # swp

    swp.swp_eff .= sum(rwu_per .* Matrix(swp_all[:, Not(:time)])', dims=1)';

    return swp
end

function get_eff_swp_rd(sim)
    # derive effective soil water potential based on root distribution
    days, dates_out = get_dates(sim);
    swp = DataFrame(date = dates_out);

    beta = sim.parametrizedSPAC.pars.root_distribution.beta;
    root_density = sim.parametrizedSPAC.soil_discretization.df.Rootden_;

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
            "Irrigation" => :brown
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
            :Precip => "Precipitation",
            :Irrig => "Irrigation");
    
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
        data(@subset(df_part_yearly_forMakie, :variable .∉ (["Precipitation","Irrigation"],))) * visual(BarPlot) +
        mapping(
            :year => "",
            :value => "Water flux per year (mm)",
            color = :variable => "") *
        # line plot of precip input
        data(@subset(df_part_yearly_forMakie, :variable .∈ (["Precipitation","Irrigation"],))) * visual(Lines)
        
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

    obs = obs[obs.date .∈ [z.date], :];

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

    cc = cor(sap_comp.trans, sap_comp.Tr);

    return cc

end

# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250814.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250814.csv", DataFrame);
met_irst = CSV.read("LWFBcal_output/metrics_irst_20251222.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_ctr_20250814.csv", DataFrame);
par_irr = CSV.read("LWFBcal_output/param_irr_20250814.csv", DataFrame);
par_irst = CSV.read("LWFBcal_output/param_irst_20251222.csv", DataFrame);

par_ctr_best, scen_ctr_best = par_best(met_ctr, par_ctr);
par_irr_best, scen_irr_best = par_best(met_irr, par_irr);
par_irst_best = par_irst[met_irst_good.scen[findmax(met_irst_good.trans_nse)[2]], :];

met_ctr[scen_ctr_best, :]
met_irr[scen_irr_best, :]

## behavioral data
# soil water content
obs_swc = CSV.read("../../Data/Pfyn/Pfyn_swat.csv", DataFrame);
#obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
#filter!(:date => >(Date(2015, 1, 1)), obs_swc); # filter out early dates
#filter!(:date => <(Date(2024, 1, 1)), obs_swc); # filter out late dates

# soil water potential
obs_swp = CSV.read("../../Data/Pfyn/PFY_swpc_corr.csv", DataFrame);
#filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
filter!(:date => <(Date(2024, 1, 1)), obs_swp); # filter out late dates

# sap flow
#obs_sap = CSV.read("../../Data/Pfyn/PFY_sap.csv", DataFrame);
obs_sap = CSV.read("../../Data/Pfyn/Pfyn_trans_2011_17.csv", DataFrame);

# leaf water potential data
obs_lwp = CSV.read("../../Data/Pfyn/Insitu_lwp.csv", DataFrame);

# isotope data
obs_soil_iso = CSV.read("../../Data/Pfyn/Pfyn_insitu_soil_iso.csv", DataFrame);
obs_xy_iso = CSV.read("../../Data/Pfyn/Pfyn_insitu_xylem_iso.csv", DataFrame);

soil_iso_depth_dict = Dict("0-10" => 5, "10-30" => 20, "30-50" => 40);
obs_soil_iso.depth = [soil_iso_depth_dict[d] for d ∈ obs_soil_iso.depth];

# separate control and irrigation scenarios
obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
obs_swc_irr = obs_swc[obs_swc.meta .== "irrigation", :]; # select irrigation treatment
obs_swc_irst = obs_swc[obs_swc.meta .== "stop", :]; # select irrigation stop treatment

obs_swp_ctr = obs_swp[obs_swp.meta .== "control", :]; # select control treatment
obs_swp_irr = obs_swp[obs_swp.meta .== "irrigated", :]; # select irrigation treatment
obs_swp_irst = obs_swp[obs_swp.meta .== "irrigation_stop", :]; # select irrigation stop treatment

obs_sap_ctr = obs_sap[obs_sap.scen .== "Control", :]; # select control treatment
obs_sap_irr = obs_sap[obs_sap.scen .== "Irrigation", :]; # select irrigation treatment
obs_sap_irst = obs_sap[obs_sap.scen .== "Irrigation stop", :]; # select irrigation stop treatment

obs_lwp_ctr = obs_lwp[obs_lwp.treatment .== "control", :]; # select control treatment
obs_lwp_irr = obs_lwp[obs_lwp.treatment .== "irrigation", :]; # select irrigation treatment
obs_lwp_irst = obs_lwp[obs_lwp.treatment .== "stop irrigation", :]; # select irrigation stop treatment

obs_soil_iso_ctr = obs_soil_iso[obs_soil_iso.treatment .== "control", :]; # select control treatment
obs_soil_iso_irr = obs_soil_iso[obs_soil_iso.treatment .== "irrigation", :]; # select irrigation treatment
obs_soil_iso_irst = obs_soil_iso[obs_soil_iso.treatment .== "irrigation stop", :]; # select irrigation stop treatment

obs_xy_iso_ctr = obs_xy_iso[obs_xy_iso.treatment .== "control", :]; # select control treatment
obs_xy_iso_irr = obs_xy_iso[obs_xy_iso.treatment .== "irrigation", :]; # select irrigation treatment
obs_xy_iso_irst = obs_xy_iso[obs_xy_iso.treatment .== "irrigation stop", :]; # select irrigation stop treatment

# irrigation dates
irr = CSV.read("../../Data/Pfyn/irrigation.csv", DataFrame);

# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2002, 1, 1), Date(2024, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2002, 1, 1), Date(2024, 12, 31), "LWFBinput/Pfyn_irrigiso_ambient/", "pfynwald", "LWFB_testrun/irrigation/", irrig=true);
sim_irst = run_LWFB90_param(par_irst_best, Date(2002, 1, 1), Date(2024, 12, 31), "LWFBinput/Pfyn_irrigiso_stop/", "pfynwald", "LWFB_testrun/irr_stop/", irrig=true);

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
draw(data(swp_comp_irst[swp_comp_irst.SWP .> -2000, :])*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Potential Comparison for Irrigation Stop Scenario", titlealign = :center)
)

draw(data(swp_comp_ctr)*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; label="Source")), facet = (; linkyaxes = :none),
    figure = (; size=(1200, 800), title="Soil Water Potential Comparison for Control Scenario", titlealign = :center)
)

draw(data(swp_comp_irr)*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Potential Comparison for Irrigation Scenario", titlealign = :center)
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
    figure = (; size=(800, 600), title="Soil Water Content Comparison for Irrigation Scenario", titlealign = :center)
)

draw(data(swc_comp_irst)*
    mapping(:date, :VWC, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Content Comparison for Irrigation Stop Scenario", titlealign = :center)
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


draw(data(stack(sap_comp_irr, Not(:date), variable_name=:Source))*
    mapping(:date, :value, color=:Source)*visual(Scatter, markersize=6),
    scales(X = (; label="")),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Irrigation Scenario", titlealign = :center)
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
    mapping(:trans, :Tr)*(visual(Scatter)+linear()),
    figure = (; size=(800, 600), title="Sap Flow Comparison for Irrigation Stop Scenario", titlealign = :center)
)

# soil isotopes

draw(data(filter!(:date => >(Date(2023,6,1)), get_soil_iso(sim_ctr)))*
    mapping(:date, :d18O, row=:depth => nonnumeric)*visual(Lines, color=:blue)+
    data(obs_soil_iso_ctr)*
    mapping(:date, :d18O, row=:depth => nonnumeric)*visual(Scatter, markersize=8),
    scales(X = (; label=""), Y= (; label="Soil δ18O (‰ VSMOW)")),
    figure = (; size=(800, 600), title="Observed Soil Water δ18O for Control Scenario", titlealign = :center)
)

# xylem isotopes

draw(data(filter!(:date => >(Date(2023,6,1)), get_xy_iso(sim_ctr)))*
    mapping(:date, :XYLEM_d18O)*visual(Lines, color=:green)+
    data(obs_xy_iso_ctr)*
    mapping(:date, :d18O)*visual(Scatter, markersize=8),
    scales(X = (; label=""), Y= (; label="Xylem δ18O (‰ VSMOW)")),
    figure = (; size=(800, 600), title="Observed Xylem δ18O for Control Scenario", titlealign = :center)
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
ctr_swp.treatment .= "control";

irr_swp = get_swp(sim_irr);
irr_swp.treatment .= "irrigation";

irst_swp = get_swp(sim_irst);
irst_swp.treatment .= "stop irrigation";

swp_comp_scen = [ctr_swp; irr_swp; irst_swp];

draw(data(swp_comp_scen)*
    mapping(:date, :SWP, color=:treatment, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SMP (kPa)"), Color = (; palette = ["#E69F00","#56B4E9","#009E73"], label="Scenario")),
    figure = (; size=(800, 600), title="Modelled Soil Water Potential Comparison across Scenarios", titlealign = :center)
)

# compare against leaf water potential
swp_comp_scen = swp_comp_scen[swp_comp_scen.date .>= Date(2022, 4, 1) .&& swp_comp_scen.date .< Date(2025,1,1), :];

draw(
    (data(swp_comp_scen)*
    mapping(:date, :SWP => (x -> x/1000), color=:depth => nonnumeric)*visual(Lines, linewidth=1.5)+
    data(obs_lwp)*
    mapping(:date, :predawn => (x -> -1*x/10))*visual(Scatter, markersize=12))*
    mapping(layout=:treatment),
    scales(X = (; label=""), Y= (; label="SWP, LWP (MPa)")),
    figure = (; size=(1200, 600), title="Comparison between Observed pre-dawn Leaf Water Potential (LWP) and Modelled Soil Water Potential (SWP)")
)

# using effective soil water potential

ctr_eff_swp = get_eff_swp(sim_ctr);
ctr_eff_swp.treatment .= "control";

irr_eff_swp = get_eff_swp(sim_irr);
irr_eff_swp.treatment .= "irrigation";

irst_eff_swp = get_eff_swp(sim_irst);
irst_eff_swp.treatment .= "stop irrigation";

eff_swp_comp = [ctr_eff_swp; irr_eff_swp; irst_eff_swp];
eff_swp_comp = eff_swp_comp[eff_swp_comp.date .>= Date(2022, 4, 1) .&& eff_swp_comp.date .< Date(2025,1,1), :];

draw(
    ((data(eff_swp_comp)*
    mapping(:date, :swp_eff => (x -> x/1000), color=:treatment)*visual(Lines, linewidth=1.5)+
    data(obs_lwp[obs_lwp.date .< Date(2025, 1, 1), :])*
    mapping(:date, :predawn => (x -> -1*x/10))*visual(Scatter, markersize=12))*
    mapping(layout=:treatment)),
    scales(X = (; label=""), Y= (; label="SWP, LWP (MPa)")),
    figure = (; size=(1200, 600), title="Comparison between Observed pre-dawn Leaf Water Potential (LWP) and Modelled Effective Soil Water Potential (SWP)")
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

ctr_rwu_tran = get_trans_def(sim_ctr);
ctr_rwu_tran.RWU, = get_RWU_centroid(sim_ctr);
ctr_rwu_tran.REW = get_REW(sim_ctr);
ctr_rwu_tran_swp = get_eff_swp(sim_ctr);
ctr_rwu_tran.SWP = ctr_rwu_tran_swp.swp_eff;
ctr_rwu_tran.scen .= "Control";

irr_rwu_tran = get_trans_def(sim_irr);
irr_rwu_tran.RWU, = get_RWU_centroid(sim_irr);
irr_rwu_tran.REW = get_REW(sim_irr);
irr_rwu_tran_swp = get_eff_swp(sim_irr);
irr_rwu_tran.SWP = irr_rwu_tran_swp.swp_eff;
irr_rwu_tran.scen .= "Irrigation";

df_rwu_tran = [ctr_rwu_tran; irr_rwu_tran];
df_rwu_tran.RWU = replace(df_rwu_tran.RWU, NaN=>missing);
df_rwu_tran.month = month.(df_rwu_tran.date);
df_rwu_tran.year = year.(df_rwu_tran.date);
df_rwu_tran = df_rwu_tran[df_rwu_tran.year .> 2002 .&& df_rwu_tran.month .> 3 .&& df_rwu_tran.month .< 11, :];

# apply bikini filter from Peters et al. (2018) with forcing data
clim = get_forcing(sim_ctr);
clim.date = Date.(clim.dates);
clim.tmean = (clim.tmax_degC .+ clim.tmin_degC) ./ 2;
clim.vpd = (0.61078 .* exp.(17.26939 .* clim.tmean ./ (clim.tmean .+ 237.3))) - clim.vappres_kPa;
clim_bikini = clim[clim.tmean .> 12 .&& clim.vpd .< 1 .&& clim.prec_mmDay .< 10, :];

df_rwu_tran_filt = dropmissing(df_rwu_tran[df_rwu_tran.date .∈ [clim_bikini.date], :]);
# df_rwu_tran = df_rwu_tran[df_rwu_tran.REW .> 0.5, :]; # filter for water availability

df_rwu_tran_clim = leftjoin(df_rwu_tran_filt, clim_bikini[:, [:date, :tmean, :vpd]], on=:date);

# mean uptake depth
med_rwu_comp = combine(groupby(df_rwu_tran_filt, :scen), :RWU .=> [median mean]);

# RWU vs transpiration
draw(
    data(df_rwu_tran_filt)*
    mapping(:trans, :RWU, color=:month => nonnumeric, layout=:scen)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWU_mean, layout=:scen)*visual(HLines, color=:black, linestyle=:dash),
    scales(Color = (; label="Month", palette = from_continuous(:seaborn_bright6)),
           X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), facet = (; linkxaxes = :none), figure = (; size=(800, 400))
)

# colored by VPD
draw(
    data(df_rwu_tran_clim)*
    mapping(:trans, :RWU, color=:vpd, layout=:scen)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWU_mean, layout=:scen)*visual(HLines, color=:black, linestyle=:dash),
    scales(X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), facet = (; linkxaxes = :none), figure = (; size=(800, 400))
)

# Trans vs VPD colored by RWU
draw(
    data(df_rwu_tran_clim)*
    mapping(:vpd, :trans, color=:SWP, layout=:scen)*visual(Scatter, alpha=.6, markersize=8),
    #data(med_rwu_comp)*mapping(:RWU_mean, layout=:scen)*visual(HLines, color=:black, linestyle=:dash),
    #scales(X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    #axis = (; yreversed = true), facet = (; linkxaxes = :none), figure = (; size=(800, 400))
)

# colored by REW
draw(
    data(df_rwu_tran_clim)*
    mapping(:trans, :RWU, color=:REW, layout=:scen)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWU_mean, layout=:scen)*visual(HLines, color=:black, linestyle=:dash),
    scales(X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), facet = (; linkxaxes = :none), figure = (; size=(800, 400))
)

# RWU vs transpiration deficit
draw(
    data(df_rwu_tran_filt)*
    mapping(:td, :RWU, color=:month => nonnumeric, layout=:scen)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWU_mean, layout=:scen)*visual(HLines, color=:black, linestyle=:dash),
    scales(Color = (; label="Month", palette = from_continuous(:seaborn_bright6)),
           X = (; label="Transpiration Deficit (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
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
depth_bins = [0, 70, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 1900, 2000];
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
    xticks=(tickdates, Dates.format.(tickdates, "u")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="Day of Year", ylabel="Root Water Uptake (mm/day)", legendtitle="RWU Depth (mm)",
    title="\nMean Daily Root Water Uptake by Depth for Control Scenario", size=(1200, 800))

p1=areaplot(ctr_rwu_swp_wide.date, ctr_rwu_swp_run, label=reshape(unique(ctr_rwu_swp_daily.depth_bin), 1, 12),
    xticks=(tickdates, Dates.format.(tickdates, "u")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="", ylabel="Root Water Uptake (mm/day)", legend=false,
    title="\nMean Daily Root Water Uptake by Depth for Control Scenario", size=(1200, 800));

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

p2=areaplot(irr_rwu_swp_wide.date, irr_rwu_swp_run, label=reshape(unique(irr_rwu_swp_daily.depth_bin), 1, 12),
    xticks=(tickdates, Dates.format.(tickdates, "u")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="Day of Year", ylabel="Root Water Uptake (mm/day)", legend=false,
    title="Mean Daily Root Water Uptake by Depth for Irrigation Scenario", size=(1200, 800));
Plots.plot(p1, p2, layout = (2,1))

# scenario comparison

comp_rwu_swp = [ctr_rwu_swp; irr_rwu_swp];

draw(data(comp_rwu_swp[comp_rwu_swp.RWU .> 0, :])*
    mapping(:SWP, :RWU, color=:depth_bin, layout=:scenario)*visual(Scatter, alpha=.5, markersize=6),
    scales(Color = (; label="Uptake Depth (mm)", palette = from_continuous(:viridis)),
           X = (; label="Soil Water Potential (kPa)"), Y= (; label="Root Water Uptake (mm/day)")), facet = (; linkxaxes = :none),
    figure = (; size=(1200, 600), title="Daily Root Water Uptake by Depth and Soil Water Potential", titlesize=20, titlealign = :center)
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

sim_ctr_wb = run_LWFB90_param(par_ctr_best, Date(2000, 1, 1), Date(2024, 12, 31), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_testrun/control/", new_folder=false, watbal=true);
sim_irr_wb = run_LWFB90_param(par_irr_best, Date(2000, 1, 1), Date(2024, 12, 31), "LWFBinput/Pfyn_irrigiso_ambient/", "pfynwald", "LWFB_testrun/irrigation/", new_folder=false, watbal=true, irrig=true);
sim_irst_wb = run_LWFB90_param(par_irst_best, Date(2000, 1, 1), Date(2024, 12, 31), "LWFBinput/Pfyn_irrigiso_stop/", "pfynwald", "LWFB_testrun/irr_stop/", new_folder=false, watbal=true, irrig=true);

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