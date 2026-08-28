## run LWFBrook90 for Pfynwald legacy paper

cd("DeepDrought/Code/Julia/")

using Pkg; Pkg.activate("."); Pkg.instantiate()
using Statistics, CategoricalArrays, RollingFunctions

# helper functions
include("run_LWFB90_param.jl");

function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...)[Not(end)];
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function mod_depth_bins(df, depth_bins, varname)
    # group depth variables according to custom depth bins
    
    colnames = names(df)[Not(1)];
    depths = [match(r"_(\d+)mm$", name).captures[1] for name in colnames];
    rename!(df, [:date, Symbol.(depths)...]);
    
    df_long = stack(df, Not(:date), variable_name=:depth, value_name=:value); # make data frame long
    df_long.depth = parse.(Int, df_long.depth); # parse depths
    df_long.depth_bin = cut(df_long.depth .- 1, depth_bins, extend=true); # assign depth bins
    
    # summarize variable by sum for transpiration or mean for soil water vars
    if varname == "TRANS"
        df_long_d = combine(groupby(df_long, [:date, :depth_bin]), :value .=> sum => :value); # group by depth bin
    else
        df_long_d = combine(groupby(df_long, [:date, :depth_bin]), :value .=> mean => :value); # group by depth bin
    end

    df_long_d.depth = depth_bins[levelcode.(df_long_d.depth_bin) .+ 1]; # re-assign depth bins to depths

    # make wide again with new depth bins
    df_wide = unstack(df_long_d[:, Not(:depth_bin)], :depth, :value, renamecols=x->Symbol(varname*"_$(x)mm"));
    
    return df_wide
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

    if (any(RWU_percent .< 0))
        rows_RWU_mmDay_onlyUptake = ifelse.(rows_RWU_mmDay.>0,rows_RWU_mmDay, 0);
        RWU_percent_onlyUptake = rows_RWU_mmDay_onlyUptake ./ sum(rows_RWU_mmDay_onlyUptake; dims = 1);
        RWU_percent = RWU_percent_onlyUptake;
    end

    row_RWU_centroid_mm = sum(RWU_percent .* y_center; dims=1);

    col_RWU_centroid_mm = reshape(row_RWU_centroid_mm, :);
    
    return col_RWU_centroid_mm, RWU_percent
end

function get_eff_swp(sim)
    # derive effective soil water potential based on RWU
    days, dates_out = get_dates(sim);
    swp = DataFrame(date = dates_out);
    
    swp.rwu_depth, rwu_per = get_RWU_centroid(sim); # RWU depth and percent

    swp_all = get_soil_(:psi, sim, days_to_read_out_d=days); # swp

    swp.swp_eff .= sum(rwu_per .* Matrix(swp_all[:, Not(:time)])', dims=1)';

    return swp
end

# behavioral parameters
par_ctr = CSV.read("LWFBcal_output/Pfyn_ctr_param_best.csv", DataFrame);
par_irst = CSV.read("LWFBcal_output/Pfyn_irst_param_best.csv", DataFrame);

# simulation dates
start_date = Date(2010, 1, 1);
end_date = Date(2024, 12, 31);

# common params
input_prefix = "pfynwald";
output_path = "LWFBoutput/";
depth_bins = [0, 70, 100, 200, 250, 400, 500, 600, 800, 1000, 1200, 1500, 1900, 2000];

# initialize scope of output
global flux_out = Any[];
global soil_out = Any[];
global flux_ctr = Any[];
global soil_ctr = Any[];
global flux_irst = Any[];
global soil_irst = Any[];

# loop through parameter sets for each scenario

for s in ["ctr", "irst"]
    
    if s == "ctr"
        # input and output
        par = par_ctr
        input_path = "LWFBinput/Pfyn_control/";
        subdir_name = "ctr";
        irr = false;
    elseif s == "irst"
        # input and output
        par = par_irst
        input_path = "LWFBinput/Pfyn_irrigiso_stop/";
        subdir_name = "irst";
        irr = false;
    end
    
    for i in 1:5
        
        par_id = par.scen[i];

        out_dir = output_path * subdir_name * string(i) * "/";

        # set parameter values
        param_i = par[i, Not(:scen)];
        
        # run the model with the current parameter set
        sim = run_LWFB90_param(param_i, start_date, end_date, input_path, input_prefix, out_dir, irrig=irr, iso=false);
        
        # retrieve output
        days, dates_to_read_out = get_dates(sim);

        # soil water data and fluxes
        z_psi = get_soil_(:psi, sim, days_to_read_out_d = days);
        z_swc = get_soil_(:theta, sim, days_to_read_out_d = days);
        z_flux = get_fluxes(sim);

        # add dates and scen column
        rename!(z_swc, :time => :date);
        rename!(z_psi, :time => :date);
        rename!(z_flux, :dates => :date);
        z_swc.date = dates_to_read_out;
        z_psi.date = dates_to_read_out;
        z_flux.date = Date.(z_flux.date);
        z_flux.param .= par_id;

        # retrieve RWU depth and effective SWP
        z_swp_rwu = get_eff_swp(sim);
        z_swp_rwu.rwu_depth_sm = runmean(z_swp_rwu.rwu_depth, 14);
        z_flux = leftjoin(z_flux, z_swp_rwu, on=:date);

        # additional flux variables
        z_flux.trans_sm = runmean(z_flux.cum_d_tran, 14);
        z_flux.td = z_flux.cum_d_ptran .- z_flux.cum_d_tran;
        z_flux.t_pet = z_flux.cum_d_tran ./  z_flux.cum_d_ptran;

        # need to make depth variables uniform for combining output

        # pull transpiration from fluxes
        z_tran = select(z_flux, (:date, r"TRANI")...);
        z_flux = z_flux[:, Not(r"TRANI")];
        z_flux = z_flux[:, Not(r"PLRFI")]; # remove extra columns from flux
        
        z_tran_mod = mod_depth_bins(z_tran, depth_bins, "TRANS");
        z_flux = leftjoin(z_flux, z_tran_mod, on=:date);
        
        # soil water data
        z_swc_mod = mod_depth_bins(z_swc, depth_bins, "theta");
        z_psi_mod = mod_depth_bins(z_psi, depth_bins, "psi");

        z_soil = leftjoin(z_swc_mod, z_psi_mod, on=:date);
        z_soil.param .= par_id;

        if i==1
            flux_out = z_flux;
            soil_out = z_soil
        else
            flux_out = [flux_out; z_flux];
            soil_out = [soil_out; z_soil];
        end
    end

    if s == "ctr"
        flux_ctr = flux_out;
        flux_ctr.scen .= "ctr";
        soil_ctr = soil_out;
        soil_ctr.scen .= "ctr";
    else
        flux_irst = flux_out;
        flux_irst.scen .= "irst";
        soil_irst = soil_out;
        soil_irst.scen .= "irst";
    end

end

# save output
CSV.write("LWFBoutput/Pfyn_ctr_legacy_flux_output.csv", flux_ctr);
CSV.write("LWFBoutput/Pfyn_ctr_legacy_soil_output.csv", soil_ctr);
CSV.write("LWFBoutput/Pfyn_irst_legacy_flux_output.csv", flux_irst);
CSV.write("LWFBoutput/Pfyn_irst_legacy_soil_output.csv", soil_irst);