# run davos simulation

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions;
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Measures, Plots; gr()
using LWFBrook90

function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
end

# input files
input_path = "./LWFBinput/Davos/";
input_prefix = "davos";

model = loadSPAC(input_path, input_prefix; simulate_isotopes = true);

sim = setup(model);
simulate!(sim);

plotforcingandstates(sim)
plotamounts(sim, :above_and_belowground, :showRWUcentroid)
plotisotopes(sim, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)

# compare against flux tower ET
flux = CSV.read("../../Data/Davos/meteo/Davos_meteo.csv", DataFrame);
select!(flux, [:date, :ET]);
flux.meta .= "Flux";

sim_ET = get_fluxes(sim);
sim_ET.date = Date.(sim_ET.dates);
sim_ET.ET = sim_ET.evap;
select!(sim_ET, [:date, :ET]);
sim_ET.meta .= "Model";

comp_ET = [flux; sim_ET];

draw(data(comp_ET)*
    mapping(:date, :ET, color = :meta)*
    visual(Scatter, markersize=6, alpha=0.6),
    scales(X = (; label=""), Color=(; label="Source")),
    figure = (; size=(1200, 800))
)

flux.evap = sim_ET.ET[1:5844];
draw(data(flux)*
    mapping(:ET, :evap)*visual(Scatter)+
    mapping([0], [1]) * visual(ABLines, color=:red, linestyle=:dash)
)

# calculate metrics
cor(flux.ET, flux.evap)
NSE(flux.evap, flux.ET)

# compare to soil water data

# observations
obs_swc = CSV.read("../../Data/Davos/Davos_SWC_cal.csv", DataFrame);
obs_swc = obs_swc[obs_swc.date .<= Date("2023-01-01"), :];
obs_swc.src .= "obs"; # add source column

obs_swp = CSV.read("../../Data/Davos/Davos_SWP_cal.csv", DataFrame);
obs_swp = obs_swp[obs_swp.date .<= Date("2023-01-01"), :];
obs_swp.src .= "obs"; # add source column

swc_d = unique(obs_swc.depth);
swp_d = unique(obs_swp.depth);

# model results
days = range(sim.ODESolution.prob.tspan...);
dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);

# soil water content
z_th = get_soil_(:theta, sim, depths_to_read_out_mm = swc_d, days_to_read_out_d = days);
select!(z_th, Not(:time));
rename!(z_th, string.(swc_d));

z_th.date = Date.(dates_out);
z_th = stack(z_th, Not(:date), variable_name = "depth", value_name="SWC");
z_th.depth = parse.(Int, z_th.depth);

filter!(:date => >=(obs_swc.date[1]), z_th);
z_th.src .= "sim"; # add source column

# combine observed and simulated data
swc_comp = [obs_swc; z_th];

# format x-axis ticks
ticks = unique(z_th[month.(z_th.date) .== 1 .&& day.(z_th.date) .== 1, :date]);
aogticks = datetimeticks(ticks, Dates.format.(ticks, "Y"));

draw(data(swc_comp)*
    mapping(:date, :SWC, color = :src, layout=:depth => nonnumeric)*
    visual(Scatter, markersize=5, alpha=0.6),
    axis = (; xticks = aogticks),
    figure = (; size=(1200, 800)))

# soil water potential
z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = swp_d, days_to_read_out_d = days);
select!(z_psi, Not(:time));
rename!(z_psi, string.(swp_d));

z_psi.date = Date.(dates_out);
z_psi = stack(z_psi, Not(:date), variable_name = "depth", value_name="SWP");
z_psi.depth = parse.(Int, z_psi.depth);

filter!(:date => >=(obs_swp.date[1]), z_psi);
z_psi.src .= "sim"; # add source column

# combine observed and simulated data
swp_comp = [obs_swp; z_psi];

draw(data(swp_comp)*
    mapping(:date, :SWP, color = :src, layout=:depth => nonnumeric)*
    visual(Scatter, markersize=5, alpha=0.6),
    axis = (; xticks = aogticks),
    figure = (; size=(1200, 800)))

# root water uptake depth
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

rwu = DataFrame(date=dates_out);
rwu.RWU, = get_RWU_centroid(sim);
rwu.RWU_sm = runmean(rwu.RWU, 14);
rwu.year = year.(rwu.date);

draw(data(rwu[rwu.year .>= 2020, :])*
    (mapping(:date, :RWU)*visual(Scatter, markersize=6)+
    mapping(:date, :RWU_sm)*visual(Lines, color="red", linewidth=3)),
    scales(X = (; label=""), Y = (; label="Root Water Uptake Depth (mm)")))


# isotopes

# xylem isotopes
xy_iso = get_states(sim);
xy_iso.date = Date.(xy_iso.dates);
select!(xy_iso, [:date, :XYLEM_d18O, :XYLEM_d2H]);

# soil isotopes
soil_iso = get_soil_([:d2H, :d18O], sim, depths_to_read_out_mm = swc_d, days_to_read_out_d = days);
select!(soil_iso, Not(:time));

# combined
sim_iso = [xy_iso soil_iso];

# plot xylem isotopes
xy_iso_long = stack(xy_iso, Not(:date));

draw(data(xy_iso_long[xy_iso_long.date .>= Date("2020-01-01"), :])*
    mapping(:date, :value, color=:variable, row=:variable)*visual(Lines),
    scales(Color=(; legend=false), X=(; label=""), Y = (; label="Xylem isotope value (‰)")),
    facet = (; linkxaxes = :all, linkyaxes = :minimal))

# plot soil isotopes
soil_iso.date = Date.(dates_out);
soil_iso_long = stack(soil_iso, Not(:date), variable_name = "depth");

var = split.(soil_iso_long.depth, "_");
var = reduce(hcat, var);
iso = var[1,:];
depth = var[3,:];

soil_iso_long.iso = iso;
soil_iso_long.depth = parse.(Int, strip.(depth, 'm'));

soil_iso_wide = unstack(soil_iso_long, [:date, :depth], :iso, :value);

# format x-axis ticks
ticks = sim_iso[month.(sim_iso.date) .== 1 .&& day.(sim_iso.date) .== 1 .&& year.(sim_iso.date) .>= 2020, :date];
aogticks = datetimeticks(ticks, Dates.format.(ticks, "Y"));

draw(data(soil_iso_long[soil_iso_long.date .>= Date("2020-01-01"), :])*
    mapping(:date, :value, color=:iso, col=:depth => nonnumeric, row=:iso)*
    visual(Lines),
    scales(X = (; label=""), Y= (; label="Isotope value (‰)")),
    facet = (; linkxaxes = :all, linkyaxes = :minimal),
    axis = (; xticks = aogticks),
    figure = (; size=(1200, 600)))

# water balance

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

simulate!(sim, save_everystep = false, saveat = 0:5844);

wb_d, wb_m, wb_y = get_water_partitioning(sim);

wb_m_plot = plot_monthly_water_partitioning(wb_m);

wb_yr_plot = plot_yearly_water_partitioning(wb_y);
wb_yr_plot
