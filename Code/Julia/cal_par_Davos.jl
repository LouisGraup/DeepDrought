# parallel calibration of LWFBrook90.jl
# using soil water data as behavioral constraints

using Distributed

addprocs(20)

@everywhere begin
    using Pkg; Pkg.activate("."); Pkg.instantiate()
    using LWFBrook90, DataFrames, Dates, CSV
    using Distributions, QuasiMonteCarlo
end

@everywhere cd("DeepDrought/Code/Julia/")

## load behavioral data and define objective function

@everywhere begin
    # behavioral data
    # soil water content
    obs_swc = CSV.read("../../Data/Davos/Davos_SWC_cal.csv", DataFrame);
    obs_swc = unstack(obs_swc, :date, :depth, :SWC, renamecols=x->Symbol("SWC_$(x)mm")); # reshape data

    # soil matric potential
    obs_swp = CSV.read("../../Data/Davos/Davos_SWP_cal.csv", DataFrame);
    obs_swp = unstack(obs_swp, :date, :depth, :SWP, renamecols=x->Symbol("SWP_$(x)mm")); # reshape data

    # ET data
    obs_ET = CSV.read("../../Data/Davos/Davos_ET.csv", DataFrame);
end

# objective function to compare model output to observed data
@everywhere function obj_fun_swc(sim, obs)

    obs = obs[obs.date .∈ [sim.dates], :];
    
    # separate observed data into different depths and remove missing values
    obs_15cm = dropmissing(obs[!, [:date, :SWC_150mm]]);
    obs_50cm = dropmissing(obs[!, [:date, :SWC_500mm]]);
    obs_80cm = dropmissing(obs[!, [:date, :SWC_800mm]]);

    # match simulated data to available dates for each depth
    sim_15cm = sim[sim.dates .∈ [obs_15cm.date], :theta_m3m3_150mm];
    sim_50cm = sim[sim.dates .∈ [obs_50cm.date], :theta_m3m3_500mm];
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
    nse15 = NSE(sim_15cm, obs_15cm.SWC_150mm);
    nse50 = NSE(sim_50cm, obs_50cm.SWC_500mm);
    nse80 = NSE(sim_80cm, obs_80cm.SWC_800mm);

    # calculate RMSE
    rmse15 = RMSE(sim_15cm, obs_15cm.SWC_150mm);
    rmse50 = RMSE(sim_50cm, obs_50cm.SWC_500mm);
    rmse80 = RMSE(sim_80cm, obs_80cm.SWC_800mm);

    # calculate correlation coefficient
    cor15 = cor(sim_15cm, obs_15cm.SWC_150mm);
    cor50 = cor(sim_50cm, obs_50cm.SWC_500mm);
    cor80 = cor(sim_80cm, obs_80cm.SWC_800mm);

    return nse15, rmse15, cor15, nse50, rmse50, cor50, nse80, rmse80, cor80
end

@everywhere function obj_fun_swp(sim, obs)

    obs = obs[obs.date .∈ [sim.dates], :];

    # separate observed data into different depths and remove missing values
    obs_15cm = dropmissing(obs[!, [:date, :SWP_150mm]]);
    obs_50cm = dropmissing(obs[!, [:date, :SWP_500mm]]);
    obs_80cm = dropmissing(obs[!, [:date, :SWP_800mm]]);

    # match simulated data to available dates for each depth
    sim_15cm = sim[sim.dates .∈ [obs_15cm.date], :psi_kPa_150mm];
    sim_50cm = sim[sim.dates .∈ [obs_50cm.date], :psi_kPa_500mm];
    sim_80cm = sim[sim.dates .∈ [obs_80cm.date], :psi_kPa_800mm];

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    # calculate NSE
    nse15 = NSE(sim_15cm, obs_15cm.SWP_150mm);
    nse50 = NSE(sim_50cm, obs_50cm.SWP_500mm);
    nse80 = NSE(sim_80cm, obs_80cm.SWP_800mm);

    return nse15, nse50, nse80
end

@everywhere function obs_fun_trans(et_comp)

    # remove missing values
    et_comp = dropmissing(et_comp);

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    nse = NSE(et_comp.evap, et_comp.ET);

    return nse

end

@everywhere function obs_fun_sap(et_comp)

    # remove missing values
    et_comp = dropmissing(et_comp);

    cc = cor(et_comp.evap, et_comp.ET);

    return cc

end

@everywhere function ET_combine(sim, obs_ET)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_trans = get_ET(sim);

    z_trans = z_trans[z_trans.date .∈ [obs_ET.date], :];
    
    et_comp = leftjoin(z_trans, obs_ET, on = :date);

    return et_comp

end

@everywhere function get_ET(sim)
    # retrieve transpiration from sim
    
    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    select!(z, :date, :evap);
    
    return z
end

### BEGIN USER INPUT ###

@everywhere begin
    ## parameter input and output paths
    # input
    input_path = "LWFBinput/Davos/";
    input_prefix = "davos";

    # output
    output_path = "LWFBcalibration/";
    subdir_name = "davos";

    ## simulation dates

    start_date = Date(2007, 1, 1);
    end_date = Date(2024, 12, 31);

end


## define calibration parameter sets

n = 1000; # number of parameter sets

# define prior parameter ranges

param = [ # store unmodified parameter ranges here for posterity
    # hydro parameters
    ("DRAIN", 0.1, 0.45), # drainage (0, 1)
    ("INFEXP", 0.1, 0.55), # infiltration exponent (0, 0.9)
    ("IDEPTH_m", 0.4, 0.7), # infiltration depth (m) (0.05, 1.0)
    # meteo parameters
    ("ALB", 0.16, 0.28), # surface albedo (0.1, 0.3)
    ("ALBSN", 0.45, 0.7), # snow surface albedo (0.4, 0.8)
    # soil parameters
    ("RSSA", 100, 600), # soil resistance (20, 1000)
    ("ths1", 0.2, 0.6), # theta_sat (0.15, 0.6)
    ("thr1", 0.01, 0.1), # theta_res (0.0, 0.1)
    ("ksat1", 200, 4200), # k_sat (100, 5000)
    ("alpha1", 7, 20), # alpha (2, 25)
    ("npar1", 1.2, 1.25), # n (1.15, 1.3)
    ("tort1", -1, 8), # tortuosity (−10, 10)
    ("ths2", 0.4, 0.6), # theta_sat (0.15, 0.6)
    ("thr2", 0.01, 0.09), # theta_res (0.0, 0.1)
    ("ksat2", 1000, 4000), # k_sat (100, 5000)
    ("alpha2", 5, 20), # alpha (2, 25)
    ("npar2", 1.15, 1.3), # n (1.15, 1.3)
    ("tort2", -9, 0), # tortuosity (−10, 10)
    ("ths3", 0.45, 0.6), # theta_sat (0.15, 0.6)
    ("thr3", 0.0, 0.08), # theta_res (0.0, 0.1)
    ("ksat3", 1000, 4000), # k_sat (100, 5000)
    ("alpha3", 2, 8), # alpha (2, 25)
    ("npar3", 1.17, 1.28), # n (1.15, 1.3)
    ("tort3", -9, -6), # tortuosity (−10, 10)
    # plant parameters
    ("CINTRL", 0.1, 0.22), # interception storage capacity per unit LAI (0.05, 0.75)
    ("FRINTLAI", 0.02, 0.04), # interception catch fraction per unit LAI (0.02, 0.2)
    ("GLMAX", 0.004, 0.006), # stomatal conductance (0.001, 0.02)
    ("CVPD", 0.5, 1.0), # vpd sensitivity (0.5, 3)
    ("R5", 60, 170), # radiation sensitivity (50, 200)
    ("T1", 6, 12), # low temperature threshold (5, 15)
    ("T2", 20, 27), # high temperature threshold (20, 35)
    ("PSICR", -2.6, -1.0), # critical water potential (-3, -1)
    ("FXYLEM", 0.1, 0.4), # aboveground xylem fraction (0.1, 0.5)
    ("MXKPL", 10.0, 30.0), # maximum plant conductivity (7, 30)
    ("MXRTLN", 2000, 3600), # maximum root length (2000, 4000)
    #("VXYLEM_mm", 5.0, 80.0), # xylem volume (5, 80)
    #("DISPERSIVITY_mm", 30.0, 50.0), # dispersivity coefficient (30, 50)
    ("MAXROOTDEPTH", -1.3, -0.85), # max rooting depth (-1.5, -0.8)
    ("BETAROOT", 0.9, 0.98) # beta root coefficient (0.9, 0.999)
];

### END USER INPUT ###


## process parameter sets
np = size(param, 1); # number of parameters

# separate parameter structure
param_names = [param[i][1] for i in 1:np];
lb = [param[i][2] for i in 1:np];
ub = [param[i][3] for i in 1:np];

nsets = n * np; # total number of samples

# expand parameter sets with LHS
param_sets = QuasiMonteCarlo.sample(nsets, lb, ub, LatinHypercubeSample());

# output parameter sets
param_out = DataFrame(param_sets', param_names); # transpose parameter sets to add column names
curDate = string(Dates.format(today(), "yyyymmdd")); # save date for consistency over multi-day calibrations 
CSV.write("LWFBcal_output/param_davos_" * curDate * ".csv", param_out);


## make output folder structure and create calibration parameter files

# create output directory if it doesn't exist
if !isdir(output_path)
    mkdir(output_path);
end

@everywhere output_prefix = $input_prefix;

# input parameter file

param_file = input_path * input_prefix * "_param.csv";
param0 = CSV.read(param_file, DataFrame);

soil_file = input_path * input_prefix * "_soil_horizons.csv";
soil0 = CSV.read(soil_file, DataFrame, skipto=3);

# function to output soil file
function output_soil_file(df, out_dir)
    # necessary to insert units row into data frame
    units = ["-","m","m","volume fraction (-)","volume fraction (-)","perMeter","-","mm per day","-","volume fraction (-)"]

    df = string.(df); # convert to string
    insert!(df, 1, units); # insert units row

    # create output file name
    output_file = out_dir  * "_soil_horizons.csv";
    # write parameter file
    CSV.write(output_file, df);
end

# soil parameter count
soil_par_count = count(contains.(param_names, "ths"));

# create empty dict for root parameter indices
root_dict = Dict("ROOTS" => 0);
global root_params = false; # boolean to check for root parameters

# loop through parameter sets and create parameter files
for i in 1:nsets
    # create parameter set
    param_set = copy(param0);
    soil_set = copy(soil0);

    for j in 1:np
        # get parameter name and value
        name = param_names[j];
        value = param_sets[j, i];

        if contains(name, "ths")
            if soil_par_count > 1
                # apply ths_volfrac for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.ths_volFrac[k] = round(value, sigdigits=3);
            else
                # apply ths_volfrac for each soil horizon
                soil_set.ths_volFrac = round.(value, sigdigits=3);
            end

        elseif contains(name, "thr")
            if soil_par_count > 1
                # apply thr_volfrac for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.thr_volFrac[k] = round(value, sigdigits=3);
            else
                # apply thr_volfrac for each soil horizon
                soil_set.thr_volFrac = round.(value, sigdigits=3);
            end

        elseif contains(name, "ksat")
            if soil_par_count > 1
                # apply additive factor to log10(ksat) for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.ksat_mmDay[k] = round(value, sigdigits=5);
            else
                # apply additive factor to log10(ksat) for each soil horizon
                soil_set.ksat_mmDay = round(value, sigdigits=5);
            end

        elseif contains(name, "alpha")
            if soil_par_count > 1
                # apply multiplier to alpha_perMeter for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.alpha_perMeter[k] = round(value, sigdigits=4);
            else
                # apply multiplier to alpha_perMeter for each soil horizon
                soil_set.alpha_perMeter = round(value, sigdigits=4);
            end
        
        elseif contains(name, "npar")
            if soil_par_count > 1
                # apply npar_ for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.npar_[k] = round(value, sigdigits=5);
            else
                # apply npar_ for each soil horizon
                soil_set.npar_ = round.(value, sigdigits=5);
            end
        
        elseif contains(name, "tort")
            if soil_par_count > 1
                # apply tort_ for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.tort_[k] = round(value, sigdigits=4);
            else
                # apply tort_ for each soil horizon
                soil_set.tort_ = round.(value, sigdigits=4);
            end

        elseif name ∈ ["BETAROOT", "MAXROOTDEPTH"]
            # save index for later
            push!(root_dict, name => j);
            global root_params = true;

        else
            # get index of parameter name in file
            idx = findall(param_set.param_id .== name)[1];
            param_set.x[idx] = round(value, sigdigits=4);
        end
    end

    # output folders
    out_dir = output_path * subdir_name * string(i) * "/";
    
    # copy folder structure to output folders
    cp(input_path, out_dir, force=true);
    
    # write parameter and soil horizons files
    CSV.write(out_dir * output_prefix * "_param.csv", param_set);
    output_soil_file(soil_set, out_dir * output_prefix);
    
end


## set up calibration runs

# dummy run for reference date
model_temp = loadSPAC(input_path, input_prefix, simulate_isotopes = false);
ref_date = Date(model_temp.reference_date);

start_index = Dates.value(start_date - ref_date);
end_index = Dates.value(end_date - ref_date);

# send variables to workers
@everywhere param_sets = $param_sets
@everywhere root_dict = $root_dict
@everywhere root_params = $root_params
@everywhere start_index = $start_index
@everywhere end_index = $end_index

## run LWFBrook90 for each parameter set
# using parallel processing
@everywhere function run_calibration(i)
    
    cal_dir = output_path * subdir_name * string(i) * "/";
    par_id = i
    
    if root_params
        # retrieve root parameter values
        betaroot = round(param_sets[root_dict["BETAROOT"], par_id], sigdigits=4);
        maxroot = round(param_sets[root_dict["MAXROOTDEPTH"], par_id], sigdigits=4);

        # run model with modified root distribution
        model = loadSPAC(cal_dir, output_prefix, simulate_isotopes = true,
        root_distribution = (beta = betaroot, z_rootMax_m=maxroot));

    else
        # run model with input files
        model = loadSPAC(cal_dir, output_prefix, simulate_isotopes = true)
    end

    # model set up
    sim = setup(model, requested_tspan=(start_index, end_index));

    # error handling
    try
        simulate!(sim);
    catch
        return (par_id, fill(0, 9), fill(0, 3), fill(0, 3)) # skip if simulation fails
    end

    ## retrieve model output

    # soil water content
    z_theta = get_soil_(:theta, sim, depths_to_read_out_mm = [150, 500, 800], days_to_read_out_d = 1:end_index);
    z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = [150, 500, 800], days_to_read_out_d = 1:end_index);

    # add dates column
    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(1:end_index,sim.parametrizedSPAC.reference_date);
    z_theta.dates = Date.(dates_to_read_out);
    z_psi.dates = Date.(dates_to_read_out);

    swc_met = obj_fun_swc(z_theta, obs_swc)
    swp_met = obj_fun_swp(z_psi, obs_swp)
    
    # transpiration
    et_comp = ET_combine(sim, obs_ET);
    et_cor = obs_fun_sap(et_comp);
    et_nse = obs_fun_trans(et_comp);
    
    z_ET = get_ET(sim);
    max_ET = maximum(z_ET.evap);

    et_met = [et_cor, et_nse, max_ET];

    return (par_id, swc_met, swp_met, et_met)

end


# parallel map for calibration runs
results = pmap(i -> run_calibration(i), 1:nsets);

# intialize metric dataframes
metrics = DataFrame(scen = Int[],
    swc_nse15 = Float64[],
    swc_rmse15 = Float64[],
    swc_cor15 = Float64[],
    swc_nse50 = Float64[],
    swc_rmse50 = Float64[],
    swc_cor50 = Float64[],
    swc_nse80 = Float64[],
    swc_rmse80 = Float64[],
    swc_cor80 = Float64[],
    swp_nse15 = Float64[],
    swp_nse50 = Float64[],
    swp_nse80 = Float64[],
    et_cor = Float64[], 
    et_nse = Float64[], 
    max_ET = Float64[]);

# loop through results
for res in results
    # retrieve scenario, parameter id, and metrics
    par_id, swc, swp, et = res;
    row = [par_id, swc..., swp..., et...];
    push!(metrics, row);
end

CSV.write("LWFBcal_output/metrics_davos_" * curDate * ".csv", metrics);
