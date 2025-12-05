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

    obs_swp = CSV.read("../../Data/Pfyn/PFY_swpc_corr.csv", DataFrame);
    filter!(:date => <(Date(2021, 1, 1)), obs_swp); # filter out late dates

    # up-scaled sap flow data
    obs_sap = CSV.read("../../Data/Pfyn/Pfyn_trans11_18.csv", DataFrame);

    # separate out irrigation scenario
    obs_swc_irr = obs_swc[obs_swc.meta .== "irrigation", :]; # select irrigation treatment
    select!(obs_swc_irr, :date, :depth, :VWC); # remove extra columns
    obs_swc_irr = unstack(obs_swc_irr, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
    sort!(obs_swc_irr, :date); # sort by date

    obs_swp_irr = obs_swp[obs_swp.meta .== "irrigation", :]; # select irrigation treatment

    obs_sap_irr = obs_sap[obs_sap.scen .== "Irrigation", :]; # select irrigation treatment
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

@everywhere function sap_comp(sim, obs_sap)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil water potential data

    z_trans = get_sap(sim);

    filter!(:date => >=(obs_sap.date[1]), z_trans);
    filter!(:date => <(Date(2018, 1, 1)), z_trans);

    obs_sap = select(obs_sap, Not(:scen));

    sap_comp = sort(leftjoin(z_trans, obs_sap, on = :date), :date);

    return sap_comp

end

@everywhere function get_sap(sim)
    # retrieve soil water potential data from sim
    days, dates_out = get_dates(sim);
    
    z = get_fluxes(sim);
    z.date = dates_out;
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

### BEGIN USER INPUT ###

@everywhere begin
    ## parameter input and output paths
    # input
    input_path = "LWFBinput/Pfyn_irrigiso_ambient/";
    input_prefix = "pfynwald";

    # output
    output_path = "LWFBcalibration/";
    subdir_name = "cal_irr";

    ## simulation dates

    start_date = Date(2000, 1, 1);
    end_date = Date(2020, 12, 31);

end


## define calibration parameter sets

n = 1000; # number of parameter sets

# define prior parameter ranges

param = [
    # hydro parameters
    ("DRAIN", 0.0, 1.0), # drainage (0, 1)
    ("INFEXP", 0.0, 0.9), # infiltration exponent (0, 0.9)
    ("IDEPTH_m", 0.05, 1.0), # infiltration depth (m) (0.05, 0.5)
    # meteo parameters
    #("ALB", 0.15, 0.3), # surface albedo (0.1, 0.3)
    #("ALBSN", 0.4, 0.8), # snow surface albedo (0.4, 0.8)
    # soil parameters
    ("RSSA", 1, 1500), # soil resistance (1, 1500)
    ("ths1", 0.5, 1.5), # multiplier on theta_sat (0.5, 1.5)
    ("ksat1", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("ths2", 0.5, 1.5), # multiplier on theta_sat (0.5, 1.5)
    ("ksat2", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    ("ths3", 0.5, 1.5), # multiplier on theta_sat (0.5, 1.5)
    ("ksat3", -0.5, 0.5), # additive factor on log10(k_sat) (-0.5, 0.5)
    # plant parameters
    #("CINTRL", 0.1, 0.75), # interception storage capacity per unit LAI (0.05, 0.75)
    ("FRINTLAI", 0.02, 0.2), # interception catch fraction per unit LAI (0.02, 0.2)
    ("GLMAX", 0.001, 0.03), # stomatal conductance (0.001, 0.03)
    ("CVPD", 1.0, 3), # vpd sensitivity (1, 3)
    ("R5", 50, 500), # radiation sensitivity (50, 400)
    #("T1", 6, 12), # low temperature threshold (5, 15)
    #("T2", 20, 35), # high temperature threshold (20, 35)
    ("PSICR", -3.0, -1.0), # critical water potential (-4, -1)
    ("FXYLEM", 0.2, 0.8), # aboveground xylem fraction (0.2, 0.8)
    ("MXKPL", 1.0, 30.0), # maximum plant conductivity (1, 30)
    ("MXRTLN", 100, 6000), # maximum root length (100, 6000)
    #("VXYLEM_mm", 1.0, 100.0), # xylem volume (1, 100)
    #("DISPERSIVITY_mm", 1.0, 100.0), # dispersivity coefficient (1, 100)
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
CSV.write("LWFBcal_output/param_irr_" * curDate * ".csv", param_out);


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
        model = loadSPAC(cal_dir, output_prefix, 
        simulate_isotopes = true, simulate_irrigation = true,
        root_distribution = (beta = betaroot, z_rootMax_m=maxroot));

    else
        # run model with input files
        model = loadSPAC(cal_dir, output_prefix, 
        simulate_isotopes = true, simulate_irrigation = true)
    end

    # model set up
    sim = setup(model, requested_tspan=(start_index, end_index));

    # error handling
    try
        simulate!(sim);
    catch
        return (par_id, fill(0, 14)) # skip if simulation fails
    end

    ## retrieve model output

    # soil water content
    z_theta = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = 1:end_index);
    z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = 1:end_index);

    # add dates column
    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(1:end_index,sim.parametrizedSPAC.reference_date);
    z_theta.dates = Date.(dates_to_read_out);
    z_psi.dates = Date.(dates_to_read_out);

    swc_met = obj_fun_swc(z_theta, obs_swc_irr)
    swp_met = obj_fun_swp(z_psi, obs_swp_irr)
    
    # transpiration
    sap_comp = sap_comp(sim, obs_sap_irr);
    sap_cor = obs_fun_sap(sap_comp);
    sap_nse = obs_fun_trans(sap_comp);

    z_trans = get_sap(sim);
    max_trans = maximum(z_trans.trans);

    # annual total
    z_trans.year = year.(z_trans.Date);
    z_trans_yr = combine(groupby(z_trans, :year), :trans => sum);
    mean_ann_trans = mean(z_trans_yr.trans_sum);

    tr_met = [sap_cor, sap_nse, max_trans, mean_ann_trans]

    return (par_id, swc_met, swp_met, tr_met)

end


# parallel map for calibration runs
results = pmap(i -> run_calibration(i), 1:nsets);

# intialize metric dataframes
metrics_irr = DataFrame(scen = Int[],
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
    trans_cor = Float64[], 
    trans_nse = Float64[], 
    max_trans = Float64[],
    ann_trans = Float64[]);

# loop through results
for res in results
    # retrieve parameter id, and metrics
    par_id, swc, swp, sap = res;
    row = [par_id, swc..., swp..., sap...];
    push!(metrics_irr, row);
end

CSV.write("LWFBcal_output/metrics_irr_" * curDate * ".csv", metrics_irr);
