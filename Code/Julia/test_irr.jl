# test irrigation developments in LWFBrook90
# first add pull request from Github
# add https://github.com/fabern/LWFBrook90.jl#add-below-canopy-irrigation

using CSV, Dates, DataFrames
using LWFBrook90

input_path = "./LWFBinput/Pfyn_irrigation_test/";
input_prefix = "pfynwald";

# no irrigation
model_noirr = loadSPAC(input_path, input_prefix; simulate_isotopes = true,
    root_distribution = (beta = 0.97, z_rootMax_m=-1.9));

sim_noirr = setup(model_noirr);
simulate!(sim_noirr)

# with irrigation
model_irr = loadSPAC(input_path, input_prefix; simulate_isotopes = true, simulate_irrigation = true,
    root_distribution = (beta = 0.97, z_rootMax_m=-1.9));

sim_irr = setup(model_irr, requested_tspan=(start_index, end_index));
simulate!(sim_irr)

# restrict dates for testing
ref_date = Date(model_irr.reference_date);
start_index = Dates.value(Date(2003, 1, 1) - ref_date);
end_index = Dates.value(Date(2003, 12, 31) - ref_date);

simulate!(sim_irr, save_everystep = false, saveat = start_index:end_index); # for watbal comparison

using Plots, Measures; gr();

pl1 = plotamounts(sim_noirr, :above_and_belowground, :showRWUcentroid)
pl1

pl2 = plotisotopes(sim_irr, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
pl2

pl2 = plotisotopes(sim_irr)


# output comparison
days = range(extrema(sim_irr.ODESolution.t)...);
dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(
    days,sim_irr.parametrizedSPAC.reference_date);

# soil water isotopes at shallow depth
iso_no_irr = get_soil_([:δ18O, :δ2H], sim_noirr, depths_to_read_out_mm=[70], days_to_read_out_d=days);
iso_irr = get_soil_([:δ18O, :δ2H], sim_irr, depths_to_read_out_mm=[70], days_to_read_out_d=days);

iso_no_irr.dates = dates_to_read_out;
iso_irr.dates = dates_to_read_out;

# get infiltration
sol_no_irr = sim_noirr.ODESolution;
sol_irr = sim_irr.ODESolution;

slfl_no_irr = [sol_no_irr(t).accum.slfl for t in days];
slfl_irr = [sol_irr(t).accum.slfl for t in days];

# combine outputs
iso_no_irr.slfl = slfl_no_irr;
iso_irr.slfl = slfl_irr;

plot(iso_no_irr.dates, iso_no_irr.slfl, label="no irr")
plot(iso_irr.dates, iso_irr.slfl, label="with irr")

plot(iso_no_irr.dates, iso_no_irr.δ18O_permil_70mm, label="no irr")
plot!(iso_irr.dates, iso_irr.δ18O_permil_70mm, label="with irr")

plot(iso_no_irr.dates, iso_no_irr.δ2H_permil_70mm, label="no irr")
plot!(iso_irr.dates, iso_irr.δ2H_permil_70mm, label="with irr")

# forcing
met = CSV.read(joinpath(input_path, "pfynwald_meteoveg.csv"), DataFrame, skipto=3);
met_irr = get_forcing(sim_irr);
flux_irr = get_fluxes(sim_irr);
met_irr.year = year.(met_irr.dates);
met_irr = met_irr[met_irr.year .== 2003, :];

# differences between meteo input and modelled fluxes
irr_diff = met_irr.irrig_mmDay .- flux_irr.cum_d_irrig;
ind_diff = findall(irr_diff .!= 0);
irr_diff[ind_diff]
met_irr.irrig_mmDay[ind_diff]'
flux_irr.cum_d_irrig[ind_diff]'

# watbal
using CairoMakie, AlgebraOfGraphics, DataFramesMeta, CategoricalArrays, Chain;

wp_d, wp_m, wp_y = get_water_partitioning(sim_irr);

# water partitioning including irrigation as precip
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

    # add irrigation to precipitation
    df_part_yr.Precip = df_part_yr.Precip .+ df_part_yr.Irrig;

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

wb_yr = plot_yearly_water_partitioning(wp_y);