## run best performing scenario from LWFBrook90.jl calibration for Vetroz site
## and examine output against observed data

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions;
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Measures, Plots; gr()

include("run_LWFB90_param.jl");

# function to retrieve best parameter set
function get_par_best(met, par)
    
    # add combined NSE metric
    met.swp_nse_com = sqrt.((met.swp_nse20 .^ 2 + met.swp_nse80 .^ 2 + met.swp_nse110 .^ 2 + met.swp_nse160 .^ 2) / 4);

    # best overall scenario
    max_idx = argmax(met.swp_nse_com);
    scen_best = met.scen[max_idx];

    par_best = par[scen_best, :];

    return par_best, scen_best
end

function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...)[Not(end)];
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function get_swc(sim; shape = "long")
    # retrieve soil water potential data from sim
    days, dates_out = get_dates(sim);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [200, 800, 1100, 1600], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="VWC");
        z_long.depth = parse.(Int, replace.(z_long.depth, r"theta_m3m3_(\d+)mm" => s"\1")) .÷ 10; # convert depth to cm from var name
        return z_long
    else
        return z
    end

end

function get_swp(sim; shape="long")
    # retrieve soil water potential data from sim
    days, dates_out = get_dates(sim);
    
    z = get_soil_(:psi, sim, depths_to_read_out_mm = [200, 800, 1100, 1600], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="SWP");
        z_long.depth = parse.(Int, replace.(z_long.depth, r"psi_kPa_(\d+)mm" => s"\1")) .÷ 10; # convert depth to cm from var name
        
        return z_long
    else
        return z
    end

end

function get_sap(sim)
    # retrieve transpiration from sim
    
    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

function get_plpsi(sim)
    # retrieve plant water potential from sim

    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    select!(z, :date, :cum_pd_plpsi, :cum_md_plpsi);

    return z
end

function ann_trans(sim)
    # calculates annual transpiration
    z = get_sap(sim);

    z.year = year.(z.date);

    z_ann = combine(groupby(z, :year), :trans => sum);

    return z_ann
end

function f_comp_swp(sim, obs_swp)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_psi = get_swp(sim);

    filter!(:date => >=(obs_swp.date[1]), z_psi);

    z_psi.src .= "sim"; # add source column
    obs_swp.src .= "obs"; # add source column

    # combine observed and simulated data
    swp_comp = [obs_swp; z_psi];

    return swp_comp

end

function f_comp_sap(sim, obs_sap)
    # sim is the LWFBrook90 simulation
    # obs is the observed sap flow data

    z_trans = get_sap(sim);

    filter!(:date => >=(obs_sap.date[1]), z_trans);

    sap_comp = sort(leftjoin(z_trans, obs_sap, on = :date), :date);

    return sap_comp

end

function f_comp_twd(sim, obs_twd)
    # sim is the LWFBrook90 simulation
    # obs is the observed dendrometer data

    z_plpsi = get_plpsi(sim);

    twd_comp = sort(innerjoin(obs_twd, z_plpsi, on = :date), :date);

    return twd_comp
end

function get_RWU_centroid(sim)
    # borrow code from LWFBrook90 package
    solu = sim.ODESolution;
    saved = sim.saved_values;

    days_to_read_out_d = saved.t;

    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2;

    # Compute RWU centroid
    rows_RWU_mmDay  = reduce(hcat, [saved.saveval[t].TRANI for t in 1:(length(days_to_read_out_d)-1)]);

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

function get_PLRF_centroid(sim)
    # borrow code from LWFBrook90 package
    solu = sim.ODESolution;
    saved = sim.saved_values;

    days_to_read_out_d = saved.t;

    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2;

    # Compute PLRF centroid
    rows_PLRF_mmDay  = reduce(hcat, [saved.saveval[t].PLRFI for t in 1:(length(days_to_read_out_d)-1)]);

    PLRF_percent = rows_PLRF_mmDay ./ sum(rows_PLRF_mmDay; dims = 1);

    if (any(PLRF_percent .< 0))
        rows_PLRF_mmDay_onlyUptake = ifelse.(rows_PLRF_mmDay.>0,rows_PLRF_mmDay, 0);
        PLRF_percent_onlyUptake = rows_PLRF_mmDay_onlyUptake ./ sum(rows_PLRF_mmDay_onlyUptake; dims = 1);
        PLRF_percent = PLRF_percent_onlyUptake;
    end

    row_PLRF_centroid_mm = sum(PLRF_percent .* y_center; dims=1);

    col_PLRF_centroid_mm = reshape(row_PLRF_centroid_mm, :);
    
    return col_PLRF_centroid_mm, PLRF_percent
end

function get_REW(sim)
    # derive relate extractable soil water based on maximum root depth
    # for now, just use normalized SWAT over all layers

    days, dates_out = get_dates(sim);
    
    #swat = get_states(sim);

    #REW = (swat.SWAT_mm .- minimum(swat.SWAT_mm)) ./ (maximum(swat.SWAT_mm) - minimum(swat.SWAT_mm));

    #return REW

    # get # root zone layers from maximum rooting depth
    max_root_depth = -1.65;

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

function get_clim(sim)

    met = get_forcing(sim);
    met.date = Date.(met.dates);

    met = select(met, :date, :vappres_kPa, :tmax_degC, :tmin_degC, :prec_mmDay);
    met.month = month.(met.date);
    met.year = year.(met.date);

    met.tmean = (met.tmax_degC .+ met.tmin_degC) ./ 2;

    # calc saturation vapor pressure
    met.Es = 0.61078 .* exp.(17.26939 .* met.tmean ./ (met.tmean .+ 237.3));
    # calc vapor pressure deficit
    met.VPD = met.Es .- met.vappres_kPa;

    return met
end

function combine_fluxes(sim)
    # get all relevant fluxes and combine into single data frame

    df_fluxes = get_fluxes(sim);
    df_fluxes.date = Date.(df_fluxes.dates);
    df_fluxes.trans = df_fluxes.cum_d_tran;
    df_fluxes.pet = df_fluxes.cum_d_ptran;
    df_fluxes.td = df_fluxes.pet .- df_fluxes.trans
    df_fluxes.PLFL = df_fluxes.cum_d_plfl;
    df_fluxes.PLRF = df_fluxes.cum_d_plrf;
    df_fluxes.plpsi_pd = df_fluxes.cum_pd_plpsi;
    df_fluxes.plpsi_md = df_fluxes.cum_md_plpsi;
    df_fluxes.PLFL_Tr_perc = df_fluxes.PLFL ./ df_fluxes.trans * 100;
    df_fluxes.RWUd, = get_RWU_centroid(sim);
    df_fluxes.PLRFd, = get_PLRF_centroid(sim);
    df_fluxes.REW = get_REW(sim);

    eff_swp = get_eff_swp(sim);
    df_fluxes.SWP = eff_swp.swp_eff;

    df_fluxes.RWUd = replace(df_fluxes.RWUd, NaN=>missing);
    df_fluxes.PLRFd = replace(df_fluxes.PLRFd, NaN=>missing);
    df_fluxes.month = month.(df_fluxes.date);
    df_fluxes.year = year.(df_fluxes.date);

    select!(df_fluxes, :date, :trans, :pet, :td, :RWU,
    :PLFL, :PLRF, :plpsi_pd, :plpsi_md, :PLFL_Tr_perc,
    :RWUd, :PLRFd, :REW, :SWP, :month, :year);

    return df_fluxes
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
            "Precipitation" => :black
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

# calibration results
met = CSV.read("LWFBcal_output/metrics_vetroz_20260526.csv", DataFrame);
par = CSV.read("LWFBcal_output/param_vetroz_20260526.csv", DataFrame);

par_best, scen_best = get_par_best(met, par);

met[scen_best, :]

## behavioral data
# soil matric potential
obs_swp = CSV.read("../../Data/Daten_Vetroz_Lorenz/SWP_daily.csv", DataFrame);
filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
obs_swp = obs_swp[obs_swp.depth .!= 1, :];

# sap flow data
obs_sap = CSV.read("../../Data/Daten_Vetroz_Lorenz/sapflow_daily.csv", DataFrame);
obs_sap = obs_sap[obs_sap.sensor_loc .== "stem", :]; # filter to stem sensors
obs_sap.sapflow .= max.(obs_sap.sapflow, 0); # set negative values to zero
select!(obs_sap, Not([:month, :year, :sensor_loc])); # drop unnecessary columns

# leaf water potential data
obs_lwp = CSV.read("../../Data/Daten_Vetroz_Lorenz/LWP_gs.csv", DataFrame);

# dendrometer data
obs_twd = CSV.read("../../Data/Daten_Vetroz_Lorenz/TWD_daily.csv", DataFrame);
obs_twd = obs_twd[obs_twd.sensor_loc .== "stem", :]; # filter to stem sensors

## run LWFBrook90.jl for all scenarios
sim = run_LWFB90_param(par_best, Date(2014, 1, 1), Date(2023, 12, 31), "LWFBinput/Vetroz/", "vetroz", "LWFB_testcap/vetroz/", new_folder=false, iso=false);

# combine observed and simulated data
# soil water potential
swp_comp = f_comp_swp(sim, obs_swp);

# sap flow
sap_comp = f_comp_sap(sim, obs_sap);

# dendro
twd_comp = f_comp_twd(sim, obs_twd);

## plot observed and simulated data

# soil water potential
draw(data(swp_comp)*
    mapping(:date, :SWP, color=:src, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="SWP (kPa)"), Color = (; label="Source")),
    figure = (; size=(800, 600), title="Soil Water Potential Comparison", titlealign = :center)
)

# soil water content
swc_sim = get_swc(sim);
draw(data(swc_sim)*
    mapping(:date, :VWC, row=:depth => nonnumeric)*visual(Scatter, markersize=4),
    scales(X = (; label="")),
    figure = (; size=(1200, 800), title="Soil Water Content", titlealign = :center)
)

# sap flow

draw(data(stack(sap_comp, Not(:date), variable_name=:Source))*
    mapping(:date, :value, color=:Source)*visual(Scatter, markersize=6),
    scales(X = (; label="")),
    figure = (; size=(800, 600), title="Sap Flow Comparison", titlealign = :center)
)

sap_comp.trans_sm = runmean(sap_comp.trans, 14);
sap_comp.sapflow_sm = runmean(sap_comp.sapflow, 14);

sap_comp_long = stack(sap_comp[:, Not([:trans_sm, :sapflow_sm])], Not(:date), variable_name=:Source);
sap_comp_long2 = stack(sap_comp[:, Not([:trans, :sapflow])], Not(:date), value_name=:value_sm, variable_name=:Source);

sap_comp_long.value_sm = sap_comp_long2.value_sm;

draw(data(sap_comp_long)*
    (mapping(:date, :value, color=:Source)*visual(Scatter, markersize=6)+
    mapping(:date, :value_sm, color=:Source)*visual(Lines)),
    scales(X = (; label="")),
    figure = (; size=(800, 600), title="Sap Flow Comparison", titlealign = :center)
)

draw(data(dropmissing(sap_comp))*
    mapping(:trans, :sapflow)*(visual(Scatter)+linear()),
    figure = (; size=(800, 600), title="Sap Flow Comparison", titlealign = :center)
)

# dendro
# predawn correlations
draw(data(twd_comp)*mapping(:TWD_pd, :cum_pd_plpsi)*visual(Scatter, markersize=6),
    scales(X = (; label="Pre-dawn TWD"), Y= (; label="Pre-dawn plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Pre-dawn TWD vs Plant Water Potential", titlealign = :center)
)

draw(data(twd_comp)*mapping(:LWP_pd, :cum_pd_plpsi)*visual(Scatter, markersize=6),
    scales(X = (; label="Pre-dawn TWD-derived LWP (MPa)"), Y= (; label="Pre-dawn plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Pre-dawn LWP vs Plant Water Potential", titlealign = :center)
)

# midday correlations
draw(data(twd_comp)*mapping(:TWD_md, :cum_md_plpsi)*visual(Scatter, markersize=6),
    scales(X = (; label="Midday TWD"), Y= (; label="Midday plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Midday TWD vs Plant Water Potential", titlealign = :center)
)

draw(data(twd_comp)*mapping(:LWP_md, :cum_md_plpsi)*visual(Scatter, markersize=6),
    scales(X = (; label="Midday TWD-derived LWP (MPa)"), Y= (; label="Midday plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Midday LWP vs Plant Water Potential", titlealign = :center)
)

# time series
draw(data(twd_comp)*
    (mapping(:date, :LWP_pd)*visual(Lines, color="black", label="LWP")+
    mapping(:date, :cum_pd_plpsi => (x -> x/1000))*visual(Lines, color="red", label="Model")),
    scales(X = (; label=""), Y= (; label="Water potential (MPa)")),
    legend = (; position = :bottom, framevisible = false),
    figure = (; size=(1200, 600), title="Pre-dawn LWP / Plant Water Potential", titlealign = :center)
)

draw(data(twd_comp)*
    (mapping(:date, :LWP_md)*visual(Lines, color="black", label="LWP")+
    mapping(:date, :cum_md_plpsi => (x -> x/1000))*visual(Lines, color="red", label="Model")),
    scales(X = (; label=""), Y= (; label="Water potential (MPa)")),
    legend = (; position = :bottom, framevisible = false),
    figure = (; size=(1200, 600), title="Midday LWP / Plant Water Potential", titlealign = :center)
)

# compare leaf water potential against effective soil water potential

eff_swp = get_eff_swp(sim);

draw(
    ((data(eff_swp[eff_swp.date .> Date(2021, 5, 1) .&& eff_swp.date .< Date(2022, 12, 31), :])*
    mapping(:date, :swp_eff => (x -> x/1000))*visual(Lines, linewidth=1.5)+
    data(obs_lwp[obs_lwp.pd_md .== "pd", :])*
    mapping(:Date, :LWP_mean => (x -> -1*x/10))*visual(Scatter, color="red", markersize=12))),
    scales(X = (; label=""), Y= (; label="SWP, LWP (MPa)")),
    figure = (; size=(1200, 600), title="Comparison between Observed pre-dawn Leaf Water Potential (LWP) and Modelled Effective Soil Water Potential (SWP)")
)

# compare RWU depth across scenarios

days, dates_out = get_dates(sim);

rwu_sim = DataFrame(date=dates_out);
rwu_sim.RWU, rwu_per = get_RWU_centroid(sim);
rwu_sim.RWU_sm = runmean(rwu_sim.RWU, 14);

rwu_sim.RWU = replace(rwu_sim.RWU, NaN=>missing);
rwu_sim.RWU_sm = replace(rwu_sim.RWU_sm, NaN=>missing);

draw(
    data(rwu_sim)*
    (mapping(:date, :RWU)*visual(Scatter, markersize=6)+
    mapping(:date, :RWU_sm)*visual(Lines, linewidth=3)),
    scales(X = (; label=""), Y= (; label="Root Water Uptake Depth (mm)")),
    figure = (; size=(1200, 600))
)

# compare fluxes

df_fluxes = combine_fluxes(sim);

# restrict to growing season
df_flux_grow = df_fluxes[df_fluxes.month .> 3 .&& df_fluxes.month .< 11, :];

# compare RWU depth and PLRF depth
draw(data(df_flux_grow)*
    mapping(:RWUd, :PLRFd, color=:SWP)*visual(Scatter, markersize=8)
)

# apply bikini filter from Peters et al. (2018) with forcing data
clim = get_clim(sim);
clim_bikini = clim[clim.tmean .> 12 .&& clim.VPD .< 1 .&& clim.prec_mmDay .< 10, :];

df_flux_filt = dropmissing(df_flux_grow[df_flux_grow.date .∈ [clim_bikini.date], :]);

df_flux_clim = leftjoin(df_flux_filt, clim_bikini[:, [:date, :tmean, :VPD]], on=:date);

# mean uptake depth
med_rwu_comp = combine(df_flux_filt, :RWUd .=> [median mean]);

# RWU vs transpiration
draw(
    data(df_flux_filt)*
    mapping(:trans, :RWUd, color=:month => nonnumeric)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWUd_mean)*visual(HLines, color=:black, linestyle=:dash),
    scales(Color = (; label="Month", palette = from_continuous(:seaborn_bright6)),
           X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), figure = (; size=(800, 400))
)

# colored by VPD
draw(
    data(df_flux_clim)*
    mapping(:trans, :RWUd, color=:VPD)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWUd_mean)*visual(HLines, color=:black, linestyle=:dash),
    scales(X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), figure = (; size=(800, 400))
)

# Trans vs VPD colored by RWU
draw(
    data(df_flux_clim)*
    mapping(:VPD, :trans, color=:SWP)*visual(Scatter, alpha=.6, markersize=8),
    #data(med_rwu_comp)*mapping(:RWUd_mean)*visual(HLines, color=:black, linestyle=:dash),
    #scales(X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    #axis = (; yreversed = true), figure = (; size=(800, 400))
)

# colored by REW
draw(
    data(df_flux_clim)*
    mapping(:trans, :RWUd, color=:REW)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWUd_mean)*visual(HLines, color=:black, linestyle=:dash),
    scales(X = (; label="Transpiration (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), figure = (; size=(800, 400))
)

# RWU vs transpiration deficit
draw(
    data(df_flux_filt)*
    mapping(:td, :RWUd, color=:month => nonnumeric)*visual(Scatter, alpha=.6, markersize=8)+
    data(med_rwu_comp)*mapping(:RWUd_mean)*visual(HLines, color=:black, linestyle=:dash),
    scales(Color = (; label="Month", palette = from_continuous(:seaborn_bright6)),
           X = (; label="Transpiration Deficit (mm/day)"), Y= (; label="Weighted-Average\nRoot Water Uptake Depth (mm)")),
    axis = (; yreversed = true), figure = (; size=(800, 400))
)

# plot VPD over time
draw(data(clim[clim.year .> 2014, :])*
    mapping(:date, :VPD)*visual(Scatter, markersize=4),
    scales(X = (; label=""), Y= (; label="VPD (kPa)")),
    figure = (; size=(1200, 600))
)

# compare plant storage contribution to transpiration against effective soil water potential
# only summer months
df_flux_filt = df_flux_filt[df_flux_filt.month .> 5 .&& df_flux_filt.month .< 9, :];

draw(data(df_flux_filt[df_flux_filt.PLFL_Tr_perc .> 10, :])*
    mapping(:SWP, :PLFL_Tr_perc, color=:PLRFd)*(smooth()+visual(Scatter)),
    scales(X = (; label="Effective Soil Water Potential (kPa)"), 
    Y= (; label="Plant Storage Contribution to Transpiration (%)"),
    Color = (; label="Weighted-average RWU depth for plant refill (mm)"))
)

# compare rwu against swp in each layer

days, dates_out = get_dates(sim);

rwu_all = get_soil_(:RWU, sim, days_to_read_out_d=days);
rwu_all.date = dates_out;
rwu_all = stack(rwu_all[:,Not(:time)], Not([:date]), value_name=:RWU);
rwu_all.depth = parse.(Int, [match(r"[0-9]+", s).match for s in rwu_all.variable]);
select!(rwu_all, Not(:variable));

swp_all = get_soil_(:psi, sim, days_to_read_out_d=days);
swp_all.date = dates_out;
swp_all = stack(swp_all[:,Not(:time)], Not([:date]), value_name=:SWP);
swp_all.depth = parse.(Int, [match(r"[0-9]+", s).match for s in swp_all.variable]);
select!(swp_all, Not(:variable));

rwu_swp = leftjoin(rwu_all, swp_all, on=[:date, :depth]);
depth_bins = [0, 70, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1500, 1650];
rwu_swp.depth_bin = cut(rwu_swp.depth .- 1, depth_bins, extend=true);

rwu_swp = rwu_swp[year.(rwu_swp.date) .> 2002, :];
rwu_swp.scenario .= "Control";

draw(data(rwu_swp[rwu_swp.RWU .> 0, :])*
    mapping(:SWP, :RWU, color=:depth_bin)*visual(Scatter, alpha=.5, markersize=6),
    scales(Color = (; label="Uptake Depth (mm)", palette = from_continuous(:viridis)),
           X = (; label="Soil Water Potential (kPa)"), Y= (; label="Root Water Uptake (mm/day)")),
    axis = (; title="Daily Root Water Uptake by Depth and Soil Water Potential", titlesize=20),
    figure = (; size=(800, 600))
)

# daily average RWU
rwu_swp.month = month.(rwu_swp.date);
rwu_swp.day = day.(rwu_swp.date);

rwu_swp_daily = combine(groupby(rwu_swp[rwu_swp.depth .!= 2000, :], [:month, :day, :depth_bin]), :RWU .=> mean);
rwu_swp_daily.date = Date.(0000, rwu_swp_daily.month, rwu_swp_daily.day); # dummy year for plotting

rwu_swp_wide = unstack(rwu_swp_daily[:, Not([:month, :day])], :depth_bin, :RWU_mean);
rwu_swp_run = hcat(map(x -> runmean(x, 7), eachcol(rwu_swp_wide[:, Not(:date)]))...);

tickdates = Date.(["0000-01-01", "0000-04-01", "0000-07-01", "0000-10-01", "0001-01-01"]);
areaplot(rwu_swp_wide.date, rwu_swp_run, label=reshape(unique(rwu_swp_daily.depth_bin), 1, 12),
    xticks=(tickdates, Dates.format.(tickdates, "u")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="Day of Year", ylabel="Root Water Uptake (mm/day)", legendtitle="RWU Depth (mm)",
    title="\nMean Daily Root Water Uptake by Depth for Control Scenario", size=(1200, 800))

p1=areaplot(rwu_swp_wide.date, rwu_swp_run, label=reshape(unique(rwu_swp_daily.depth_bin), 1, 12),
    xticks=(tickdates, Dates.format.(tickdates, "u")), margin=10mm, palette=palette(:viridis, 12),
    xlabel="", ylabel="Root Water Uptake (mm/day)", legend=false,
    title="\nMean Daily Root Water Uptake by Depth for Control Scenario", size=(1200, 800));


# weighted-average soil water potential

# control
rwu_swp_med = get_eff_swp(sim);
rwu_swp_med = leftjoin(rwu_swp_med, get_sap(sim), on=:date);
rwu_swp_med = rwu_swp_med[rwu_swp_med.date .>= Date(2003, 1, 1), :];
rwu_swp_med.month = month.(rwu_swp_med.date);
rwu_swp_med.scenario .= "Control";

draw(data(rwu_swp_med)*mapping(:date, :swp_eff)*visual(Lines))

draw(data(rwu_swp_med)*
    mapping(:swp_eff, :trans, color=:RWU)*visual(Scatter, alpha=.5, markersize=6),
    scales(X = (; label="Weighted-Average Soil Water Potential (kPa)"), Y= (; label="Transpiration (mm/day)")),
    axis = (; title="Transpiration vs. Effective SWP", titlesize=16)
)

draw(data(rwu_swp_med)*
    mapping(:swp_eff, :RWU, color=:month => nonnumeric)*visual(Scatter, alpha=.5, markersize=6),
    scales(Color = (; palette = from_continuous(:seaborn_colorblind6)),
        X = (; label="Weighted-Average Soil Water Potential (kPa)"), Y= (; label="Weighted-Average Root Water Uptake Depth (mm)")),
    axis = (; yreversed = true, title="Weighted-Average RWU Depth vs. SWP", titlesize=16)
)


# water balance modelling

sim_wb = run_LWFB90_param(par_best, Date(2014, 1, 1), Date(2023, 12, 31), "LWFBinput/Vetroz/", "vetroz", "LWFB_testcap/vetroz/", new_folder=false, watbal=true, iso=false);

# water partitioning

# control
wb_d, wb_m, wb_y = get_water_partitioning(sim_wb);
wb_m = wb_m[wb_m.year .== 2014, :]; # filter out single years

wb_yr_plot = plot_yearly_water_partitioning(wb_y);
wb_yr_plot