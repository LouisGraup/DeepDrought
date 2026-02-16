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
    # soil water content and soil matric potential
    obs_swc = CSV.read("../../Data/Pfyn/Pfyn_swat.csv", DataFrame);
    filter!(:date => >=(Date(2004, 1, 1)), obs_swc); # filter out early dates
    filter!(:date => <(Date(2021, 1, 1)), obs_swc); # filter out late dates

    obs_swp = CSV.read("../../Data/Pfyn/Pfyn_swp.csv", DataFrame);
    filter!(:date => >=(Date(2016, 1, 1)), obs_swp); # filter out early dates
    filter!(:date => <(Date(2021, 1, 1)), obs_swp); # filter out late dates

    # up-scaled sap flow data
    obs_sap = CSV.read("../../Data/Pfyn/Pfyn_trans_2011_17.csv", DataFrame);

    # isotope data
    obs_soil_iso = CSV.read("../../Data/Pfyn/Pfyn_insitu_soil_iso.csv", DataFrame);
    obs_xy_iso = CSV.read("../../Data/Pfyn/Pfyn_insitu_xylem_iso.csv", DataFrame);

    # separate control and irrigation scenarios
    obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
    select!(obs_swc_ctr, :date, :depth, :VWC); # remove extra columns
    obs_swc_ctr = unstack(obs_swc_ctr, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
    sort!(obs_swc_ctr, :date); # sort by date

    obs_swp_ctr = obs_swp[obs_swp.meta .== "control", :]; # select control treatment
    select!(obs_swp_ctr, :date, :depth, :SWP); # remove extra columns
    obs_swp_ctr = unstack(obs_swp_ctr, :date, :depth, :SWP, renamecols=x->Symbol("SWP_$(x)cm")); # reshape data
    sort!(obs_swp_ctr, :date); # sort by date
    
    obs_sap_ctr = obs_sap[obs_sap.scen .== "Control", :]; # select control treatment

    obs_soil_iso_ctr = obs_soil_iso[obs_soil_iso.treatment .== "control", :]; # select control treatment
    obs_soil_iso_ctr = select(obs_soil_iso_ctr, :date, :depth, :d18O); # remove extra columns
    obs_soil_iso_ctr = unstack(obs_soil_iso_ctr, :date, :depth, :d18O); # reshape data
    rename!(obs_soil_iso_ctr, [:date, :d18O_40cm, :d18O_20cm, :d18O_5cm]); # rename columns

    obs_xy_iso_ctr = obs_xy_iso[obs_xy_iso.treatment .== "control", :]; # select control treatment
    select!(obs_xy_iso_ctr, :date, :d18O); # remove extra columns
end

# objective function to compare model output to observed data
@everywhere function obj_fun_swc(sim, obs)

    obs = obs[obs.date .∈ [sim.dates], :];
    
    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :VWC_10cm]]);
    obs_40cm = dropmissing(obs[!, [:date, :VWC_40cm]]);
    obs_60cm = dropmissing(obs[!, [:date, :VWC_60cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :VWC_80cm]]);

    filter!(:date => >(Date(2015, 1, 1)), obs_80cm); # filter out early dates

    # match simulated data to available dates for each depth
    sim_10cm = sim[sim.dates .∈ [obs_10cm.date], :theta_m3m3_100mm];
    sim_40cm = sim[sim.dates .∈ [obs_40cm.date], :theta_m3m3_400mm];
    sim_60cm = sim[sim.dates .∈ [obs_60cm.date], :theta_m3m3_600mm];
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
    nse10 = NSE(sim_10cm, obs_10cm.VWC_10cm);
    nse40 = NSE(sim_40cm, obs_40cm.VWC_40cm);
    nse60 = NSE(sim_60cm, obs_60cm.VWC_60cm);
    nse80 = NSE(sim_80cm, obs_80cm.VWC_80cm);

    # calculate RMSE
    rmse10 = RMSE(sim_10cm, obs_10cm.VWC_10cm);
    rmse40 = RMSE(sim_40cm, obs_40cm.VWC_40cm);
    rmse60 = RMSE(sim_60cm, obs_60cm.VWC_60cm);
    rmse80 = RMSE(sim_80cm, obs_80cm.VWC_80cm);

    return nse10, rmse10, nse40, rmse40, nse60, rmse60, nse80, rmse80
end

@everywhere function obj_fun_swp(sim, obs)

    # separate observed data into different depths and remove missing values
    obs_10cm = dropmissing(obs[!, [:date, :SWP_10cm]]);
    obs_80cm = dropmissing(obs[!, [:date, :SWP_80cm]]);

    # match simulated data to available dates for each depth
    sim_10cm = sim[sim.dates .∈ [obs_10cm.date], :psi_kPa_100mm];
    sim_80cm = sim[sim.dates .∈ [obs_80cm.date], :psi_kPa_800mm];

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

@everywhere function obs_fun_trans(sap_comp)

    # remove missing values
    sap_comp = dropmissing(sap_comp);

    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2))
        return nse
    end

    nse = NSE(sap_comp.trans, sap_comp.Tr);

    return nse

end

@everywhere function obs_fun_sap(sap_comp)

    # remove missing values
    sap_comp = dropmissing(sap_comp);

    cc = cor(sap_comp.trans, sap_comp.Tr);

    return cc

end

@everywhere function obj_fun_soil_iso(sim, obs)

    # separate observed data into different depths and remove missing values
    obs_5cm = dropmissing(obs[!, [:date, :d18O_5cm]]);
    obs_20cm = dropmissing(obs[!, [:date, :d18O_20cm]]);
    obs_40cm = dropmissing(obs[!, [:date, :d18O_40cm]]);
    
    # match simulated data to available dates for each depth
    sim_5cm = sim[sim.dates .∈ [obs_5cm.date], :d18O_permil_50mm];
    sim_20cm = sim[sim.dates .∈ [obs_20cm.date], :d18O_permil_200mm];
    sim_40cm = sim[sim.dates .∈ [obs_40cm.date], :d18O_permil_400mm];

    function RMSE(sim, obs)
        # calculate Root Mean Square Error
        rmse = sqrt(sum((obs .- sim).^2) / length(obs))
        return rmse
    end

    # calculate RMSE
    rmse5 = RMSE(sim_5cm, obs_5cm.d18O_5cm);
    rmse20 = RMSE(sim_20cm, obs_20cm.d18O_20cm);
    rmse40 = RMSE(sim_40cm, obs_40cm.d18O_40cm);

    return rmse5, rmse20, rmse40
end

@everywhere function obj_fun_xy_iso(sim, obs)

    # remove missing values
    obs = dropmissing(obs);

    # match simulated data to available dates
    sim_xy = sim[sim.date .∈ [obs.date], :XYLEM_d18O];

    function RMSE(sim, obs)
        # calculate Root Mean Square Error
        rmse = sqrt(sum((obs .- sim).^2) / length(obs))
        return rmse
    end

    # calculate RMSE
    rmse = RMSE(sim_xy, obs.d18O);

    return rmse
end

@everywhere function sap_combine(sim, obs_sap)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_trans = get_sap(sim);

    z_trans = z_trans[z_trans.date .∈ [obs_sap.date], :];
    
    sap_comp = leftjoin(z_trans, obs_sap[:, Not(:scen)], on = :date);

    return sap_comp

end

@everywhere function get_sap(sim)
    # retrieve transpiration from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    
    z = get_fluxes(sim);
    z.date = Date.(dates_out);
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

@everywhere function get_xy_iso(sim)
    # retrieve xylem isotope data from sim
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    
    z = get_states(sim);
    z.date = Date.(dates_out);
    select!(z, :date, :XYLEM_d18O, :XYLEM_d2H);
    
    return z
end

### BEGIN USER INPUT ###

@everywhere begin
    ## parameter input and output paths
    # input
    input_path = "LWFBinput/Pfyn_control/";
    input_prefix = "pfynwald";

    # output
    output_path = "LWFBcalibration/";
    subdir_name = "cal_ctr";

    ## simulation dates

    start_date = Date(2002, 1, 1);
    end_date = Date(2024, 12, 31);

end


## define calibration parameter sets

n = 500; # number of parameter sets

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
    ("RSSA", 1, 1500), # soil resistance (1, 1500)
    ("ths1", 0.3, 1.5), # multiplier on theta_sat (0.5, 1.5)
    ("ksat1", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("ths2", 0.3, 1.5), # multiplier on theta_sat (0.5, 1.5)
    ("ksat2", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("ths3", 0.3, 1.5), # multiplier on theta_sat (0.5, 1.5)
    ("ksat3", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    # plant parameters
    #("CINTRL", 0.1, 0.75), # interception storage capacity per unit LAI (0.05, 0.75)
    ("FRINTLAI", 0.02, 0.2), # interception catch fraction per unit LAI (0.02, 0.2)
    ("GLMAX", 0.001, 0.03), # stomatal conductance (0.001, 0.03)
    ("CVPD", 1.0, 3.0), # vpd sensitivity (1, 3)
    ("R5", 50, 400), # radiation sensitivity (50, 400)
    #("T1", 6, 12), # low temperature threshold (5, 15)
    #("T2", 20, 35), # high temperature threshold (20, 35)
    ("PSICR", -3, -1.0), # critical water potential (-4, -1)
    ("FXYLEM", 0.2, 0.8), # aboveground xylem fraction (0.2, 0.8)
    ("MXKPL", 1.0, 30.0), # maximum plant conductivity (1, 30)
    ("MXRTLN", 100, 6000), # maximum root length (100, 6000)
    ("VXYLEM_mm", 5.0, 80.0), # xylem volume (5, 80)
    ("DISPERSIVITY_mm", 30.0, 50.0), # dispersivity coefficient (30, 50)
    ("MAXROOTDEPTH", -2.0, -0.5), # max rooting depth (-5, -0.5)
    ("BETAROOT", 0.9, 0.999) # beta root coefficient (0.8, 1.0)
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
CSV.write("LWFBcal_output/param_ctr_" * curDate * ".csv", param_out);


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
                # apply multiplier to ths_volfrac for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.ths_volFrac[k] = soil_set.ths_volFrac[k] * value;
            else
                # apply multiplier to ths_volfrac for each soil horizon
                soil_set.ths_volFrac = soil_set.ths_volFrac * value;
            end

        elseif contains(name, "ksat")
            if soil_par_count > 1
                # apply additive factor to log10(ksat) for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.ksat_mmDay[k] = 10 .^ (log10.(soil_set.ksat_mmDay[k]) .+ value);
            else
                # apply additive factor to log10(ksat) for each soil horizon
                soil_set.ksat_mmDay = 10 .^ (log10.(soil_set.ksat_mmDay) .+ value);
            end

        elseif contains(name, "alpha")
            if soil_par_count > 1
                # apply multiplier to alpha_perMeter for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.alpha_perMeter[k] = soil_set.alpha_perMeter[k] * value;
            else
                # apply multiplier to alpha_perMeter for each soil horizon
                soil_set.alpha_perMeter = soil_set.alpha_perMeter * value;
            end
        
        elseif contains(name, "npar")
            if soil_par_count > 1
                # apply multiplier to npar_ for specific soil horizon
                k = parse(Int, name[end]); # extract horizon number from name
                soil_set.npar_[k] = soil_set.npar_[k] * value;
            else
                # apply multiplier to npar_ for each soil horizon
                soil_set.npar_ = soil_set.npar_ * value;
            end
        
        elseif name ∈ ["BETAROOT", "MAXROOTDEPTH"]
            # save index for later
            push!(root_dict, name => j);
            global root_params = true;

        else
            # get index of parameter name in file
            idx = findall(param_set.param_id .== name)[1];
            param_set.x[idx] = value;
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
model_temp = loadSPAC(input_path, input_prefix);
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
        betaroot = param_sets[root_dict["BETAROOT"], par_id];
        maxroot = param_sets[root_dict["MAXROOTDEPTH"], par_id];

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
        return (par_id, fill(0, 8), fill(0, 2), fill(0, 4), fill(0, 4)) # skip if simulation fails
    end

    ## retrieve model output

    # soil water content
    z_theta = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = 1:end_index);
    z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = 1:end_index);

    # isotopes
    z_soil_iso = get_soil_(:d18O, sim, depths_to_read_out_mm = [50, 200, 400], days_to_read_out_d = 1:end_index);
    z_xy_iso = get_xy_iso(sim);

    # add dates column
    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(1:end_index,sim.parametrizedSPAC.reference_date);
    z_theta.dates = Date.(dates_to_read_out);
    z_psi.dates = Date.(dates_to_read_out);
    z_soil_iso.dates = Date.(dates_to_read_out);

    swc_met = obj_fun_swc(z_theta, obs_swc_ctr)
    swp_met = obj_fun_swp(z_psi, obs_swp_ctr)
    
    # isotopes
    iso_soil_met = obj_fun_soil_iso(z_soil_iso, obs_soil_iso_ctr);
    iso_xy_rmse = obj_fun_xy_iso(z_xy_iso, obs_xy_iso_ctr);

    iso_met = [iso_soil_met..., iso_xy_rmse];

    # transpiration
    sap_comp = sap_combine(sim, obs_sap_ctr);
    sap_cor = obs_fun_sap(sap_comp);
    sap_nse = obs_fun_trans(sap_comp);
    
    z_trans = get_sap(sim);
    max_trans = maximum(z_trans.trans);

    # annual total
    z_trans.year = year.(z_trans.date);
    z_trans_yr = combine(groupby(z_trans, :year), :trans => sum);
    mean_ann_trans = mean(z_trans_yr.trans_sum);

    tr_met = [sap_cor, sap_nse, max_trans, mean_ann_trans];

    return (par_id, swc_met, swp_met, iso_met, tr_met)

end


# parallel map for calibration runs
results = pmap(i -> run_calibration(i), 1:nsets);

# intialize metric dataframes
metrics = DataFrame(scen = Int[],
    swc_nse10 = Float64[],
    swc_rmse10 = Float64[],
    swc_nse40 = Float64[],
    swc_rmse40 = Float64[],
    swc_nse60 = Float64[],
    swc_rmse60 = Float64[],
    swc_nse80 = Float64[],
    swc_rmse80 = Float64[],
    swp_nse10 = Float64[],
    swp_nse80 = Float64[],
    iso_rmse5 = Float64[],
    iso_rmse20 = Float64[],
    iso_rmse40 = Float64[],
    iso_rmse_xy = Float64[],
    trans_cor = Float64[], 
    trans_nse = Float64[], 
    max_trans = Float64[], 
    ann_trans = Float64[])

# loop through results
for res in results
    # retrieve scenario, parameter id, and metrics
    par_id, swc, swp, iso, sap = res;
    row = [par_id, swc..., swp..., iso..., sap...];
    push!(metrics, row);
end

CSV.write("LWFBcal_output/metrics_ctr_" * curDate * ".csv", metrics);
