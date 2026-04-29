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
    # soil matric potential
    obs_swp = CSV.read("../../Data/Daten_Vetroz_Lorenz/SWP_daily.csv", DataFrame);
    filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
    filter!(:date => <(Date(2021, 1, 1)), obs_swp); # filter out late dates
    obs_swp = unstack(obs_swp, :date, :depth, :SWP, renamecols=x->Symbol("SWP_$(x)cm")); # reshape data

    # sap flow data
    obs_sap = CSV.read("../../Data/Daten_Vetroz_Lorenz/sapflow_daily.csv", DataFrame);
    obs_sap = obs_sap[obs_sap.month .∈ [5:11], :]; # filter to growing season
    select!(obs_sap, Not(:month));
end

# objective function to compare model output to observed data
@everywhere function obj_fun_swp(sim, obs)

    # separate observed data into different depths and remove missing values
    obs_20cm = dropmissing(obs[!, [:date, :SWP_20cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :SWP_80cm]]);
    obs_110cm = dropmissing(obs[!, [:date, :SWP_110cm]]);
    obs_160cm = dropmissing(obs[!, [:date, :SWP_160cm]]);

    # match simulated data to available dates for each depth
    sim_20cm = sim[sim.dates .∈ [obs_20cm.date], :psi_kPa_200mm];
    sim_80cm = sim[sim.dates .∈ [obs_80cm.date], :psi_kPa_800mm];
    sim_110cm = sim[sim.dates .∈ [obs_110cm.date], :psi_kPa_1100mm];
    sim_160cm = sim[sim.dates .∈ [obs_160cm.date], :psi_kPa_1600mm];

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    # calculate NSE
    nse20 = NSE(sim_20cm, obs_20cm.SWP_20cm);
    nse80 = NSE(sim_80cm, obs_80cm.SWP_80cm);
    nse110 = NSE(sim_110cm, obs_110cm.SWP_110cm);
    nse160 = NSE(sim_160cm, obs_160cm.SWP_160cm);

    return nse20, nse80, nse110, nse160
end

@everywhere function obs_fun_sap(sap_comp)

    # remove missing values
    sap_comp = dropmissing(sap_comp);

    cc = cor(sap_comp.trans, sap_comp.sapflow);

    return cc

end

@everywhere function sap_combine(sim, obs_sap)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_trans = get_sap(sim);

    z_trans = z_trans[z_trans.date .∈ [obs_sap.date], :];
    
    sap_comp = leftjoin(z_trans, obs_sap, on = :date);

    return sap_comp

end

@everywhere function get_sap(sim)
    # retrieve transpiration from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    
    sol = sim.ODESolution;
    trans = [sol(t).accum.cum_d_tran for t in days];

    z = DataFrame(date = Date.(dates_out), trans = trans);
    
    return z
end

### BEGIN USER INPUT ###

@everywhere begin
    ## parameter input and output paths
    # input
    input_path = "LWFBinput/Vetroz/";
    input_prefix = "vetroz";

    # output
    output_path = "LWFBcalibration/";
    subdir_name = "vetroz";

    ## simulation dates

    start_date = Date(2014, 1, 1);
    end_date = Date(2023, 12, 31);

end


## define calibration parameter sets

n = 1000; # number of parameter sets

# define prior parameter ranges

param = [
    # hydro parameters
    ("DRAIN", 0.0, 1.0), # drainage (0, 1)
    ("INFEXP", 0.0, 0.9), # infiltration exponent (0, 0.9)
    ("IDEPTH_m", 0.05, 1.0), # infiltration depth (m) (0.05, 1.0)
    # meteo parameters
    #("ALB", 0.15, 0.3), # surface albedo (0.1, 0.3)
    #("ALBSN", 0.4, 0.8), # snow surface albedo (0.4, 0.8)
    # soil parameters
    ("RSSA", 20, 1000), # soil resistance (20, 1000)
    ("ths1", 0.15, 0.6), # theta_sat (0.15, 0.6)
    ("thr1", 0.0, 0.1), # theta_res (0.0, 0.1)
    ("ksat1", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("alpha1", 0.5, 1.5), # multiplier on alpha (0.5, 1.5)
    ("npar1", 1.15, 1.3), # n (1.15, 1.3)
    ("ths2", 0.15, 0.6), # theta_sat (0.15, 0.6)
    ("thr2", 0.0, 0.1), # theta_res (0.0, 0.1)
    ("ksat2", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("alpha2", 0.5, 1.5), # multiplier on alpha (0.5, 1.5)
    ("npar2", 1.15, 1.3), # n (1.15, 1.3)
    ("ths3", 0.15, 0.6), # theta_sat (0.15, 0.6)
    ("thr3", 0.0, 0.1), # theta_res (0.0, 0.1)
    ("ksat3", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("alpha3", 0.5, 1.5), # multiplier on alpha (0.5, 1.5)
    ("npar3", 1.15, 1.3), # n (1.15, 1.3)
    # plant parameters
    ("CINTRL", 0.05, 0.75), # interception storage capacity per unit LAI (0.05, 0.75)
    ("FRINTLAI", 0.02, 0.2), # interception catch fraction per unit LAI (0.02, 0.2)
    ("GLMAX", 0.001, 0.02), # stomatal conductance (0.001, 0.02)
    ("CVPD", 0.5, 3.0), # vpd sensitivity (0.5, 3)
    ("R5", 50, 200), # radiation sensitivity (50, 200)
    ("T1", 5, 15), # low temperature threshold (5, 15)
    ("T2", 20, 35), # high temperature threshold (20, 35)
    ("PSICR", -4.0, -1.0), # critical water potential (-3, -1)
    ("FXYLEM", 0.1, 0.5), # aboveground xylem fraction (0.1, 0.5)
    ("MXKPL", 7.0, 30.0), # maximum plant conductivity (7, 30)
    ("MXRTLN", 2000, 4000), # maximum root length (2000, 4000)
    #("MAXROOTDEPTH", -2.0, -0.8), # max rooting depth (-2, -0.8)
    #("BETAROOT", 0.9, 0.999) # beta root coefficient (0.9, 0.999)
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
CSV.write("LWFBcal_output/param_vetroz_" * curDate * ".csv", param_out);


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
                soil_set.ksat_mmDay[k] = round(10 .^ (log10.(soil_set.ksat_mmDay[k]) .+ value), sigdigits=5);
            else
                # apply additive factor to log10(ksat) for each soil horizon
                soil_set.ksat_mmDay = round.(10 .^ (log10.(soil_set.ksat_mmDay) .+ value), sigdigits=5);
            end

        elseif contains(name, "alpha")
            if soil_par_count > 1
                # apply multiplier to alpha_perMeter for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.alpha_perMeter[k] = round(soil_set.alpha_perMeter[k] * value, sigdigits=4);
            else
                # apply multiplier to alpha_perMeter for each soil horizon
                soil_set.alpha_perMeter = round.(soil_set.alpha_perMeter * value, sigdigits=4);
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
        model = loadSPAC(cal_dir, output_prefix, simulate_isotopes = false,
        root_distribution = (beta = betaroot, z_rootMax_m=maxroot));

    else
        # run model with input files
        model = loadSPAC(cal_dir, output_prefix, simulate_isotopes = false)
    end

    # model set up
    sim = setup(model, requested_tspan=(start_index, end_index));

    # error handling
    try
        simulate!(sim);
    catch
        return (par_id, fill(0, 4), fill(0, 2)) # skip if simulation fails
    end

    ## retrieve model output

    # soil water
    z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = [200, 800, 1100, 1600], days_to_read_out_d = start_index:end_index);

    # add dates column
    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(start_index:end_index,sim.parametrizedSPAC.reference_date);
    z_psi.dates = Date.(dates_to_read_out);

    swp_met = obj_fun_swp(z_psi, obs_swp)
    
    # transpiration
    sap_comp = sap_combine(sim, obs_sap);
    sap_cor = obs_fun_sap(sap_comp);
    
    z_trans = get_sap(sim);
    max_trans = maximum(z_trans.trans);

    tr_met = [sap_cor, max_trans];

    return (par_id, swp_met, tr_met)

end


# parallel map for calibration runs
results = pmap(i -> run_calibration(i), 1:nsets);

# intialize metric dataframes
metrics = DataFrame(scen = Int[],
    swp_nse20 = Float64[],
    swp_nse80 = Float64[],
    swp_nse110 = Float64[],
    swp_nse160 = Float64[],
    trans_cor = Float64[], 
    max_trans = Float64[]);

# loop through results
for res in results
    # retrieve scenario, parameter id, and metrics
    par_id, swp, sap = res;
    row = [par_id, swp..., sap...];
    push!(metrics, row);
end

CSV.write("LWFBcal_output/metrics_vetroz_" * curDate * ".csv", metrics);
