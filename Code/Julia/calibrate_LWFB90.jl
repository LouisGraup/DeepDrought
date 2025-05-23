# calibration and sensitivity analysis of LWFBrook90.jl
# using soil water data as behavioral constraints

using LWFBrook90
using CSV, DataFrames, Dates
using Distributions, QuasiMonteCarlo
using Distributed
using Plots; gr()

## load behavioral data and define objective function

# behavioral data
obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
obs_swc = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
select!(obs_swc, :date, :depth, :VWC); # remove extra columns
filter!(:date => >(Date(2003, 4, 1)), obs_swc); # filter out early dates

obs_swc = unstack(obs_swc, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
sort!(obs_swc, :date); # sort by date

# objective function to compare model output to observed data
function obj_fun(sim, obs)

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

### BEGIN USER INPUT ###

## parameter input and output paths
# input
input_path = "LWFBinput/";
input_prefix = "pfynwald";

# output
output_path = "LWFBcalibration/";
subdir_name = "cal";

## simulation dates
start_date = Date(2000, 1, 1);
end_date = Date(2020, 12, 31);

## define calibration parameter sets

n = 200; # number of parameter sets

# define prior parameter ranges

param = [
    # hydro parameters
    ("DRAIN", 0.0, 1.0), # drainage
    ("BYPAR", 0.0, 1.0), # bypass flow
    # meteo parameters
    ("ALB", 0.1, 0.3), # surface albedo
    # soil parameters
    ("RSSA", 1.0, 1500.0), # soil resistance
    ("ths", 0.5, 2.0), # multiplier on theta_sat
    ("ksat", -0.5, 0.5), # additive factor on log10(k_sat)
    # plant parameters
    ("CINTRL", 0.05, 0.75), # interception storage capacity per unit LAI
    ("FRINTLAI", 0.02, 0.2), # interception catch fraction per unit LAI
    ("GLMAX", 0.001, 0.03), # stomatal conductance
    ("CVPD", 1.0, 3.0), # vpd sensitivity
    ("R5", 50, 400), # radiation sensitivity
    ("T1", 5, 15), # low temperature threshold
    ("T2", 20, 35), # high temperature threshold
    ("PSICR", -4.0, -1.0), # critical water potential
    ("FXYLEM", 0.2, 0.8), # aboveground xylem fraction
    ("MXKPL", 1.0, 30.0), # maximum plant conductivity
    ("VXYLEM_mm", 1.0, 100.0), # xylem volume
    ("DISPERSIVITY_mm", 1.0, 100.0), # dispersivity coefficient
    ("MAXROOTDEPTH", -5.0, -0.5), # max rooting depth
    ("BETAROOT", 0.8, 1.0) # beta root coefficient
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
param_out = DataFrame(param_sets', param_names);
CSV.write("LWFBcal_output/param_" * string(Dates.format(now(), "yyyymmdd")) * ".csv", param_out);

## make output folder structure and create calibration parameter files

# create output directory if it doesn't exist
if !isdir(output_path)
    mkdir(output_path);
end

output_prefix = input_prefix;

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

# create empty dict for root parameter indices
root_dict = Dict("ROOTS" => 0);

# loop through parameter sets and create parameter files
for i in 1:nsets
    # create parameter set
    param_set = copy(param0);
    soil_set = copy(soil0);

    for j in 1:np
        # get parameter name and value
        name = param_names[j];
        value = param_sets[j, i];

        if name == "ths"
            # apply multiplier to ths_volfrac for each soil horizon
            soil_set.ths_volFrac = soil_set.ths_volFrac * value;

        elseif name == "ksat"
            # apply additive factor to log10(ksat) for each soil horizon
            soil_set.ksat_mmDay = 10 .^ (log10.(soil_set.ksat_mmDay) .+ value);

        elseif name ∈ ["BETAROOT", "MAXROOTDEPTH"]
            # save index for later
            push!(root_dict, name => j);

        else
            # get index of parameter name in file
            idx = findall(param_set.param_id .== name)[1];
            param_set.x[idx] = value;
        end
    end

    # output folder
    out_dir = output_path * subdir_name * string(i) * "/";
    # copy folder structure to output folder
    cp(input_path, out_dir, force=true);

    # create output file name
    output_param_file = out_dir * output_prefix * "_param.csv";

    # write parameter and soil horizons file
    CSV.write(output_param_file, param_set);
    output_soil_file(soil_set, out_dir * output_prefix);
end

## set up calibration runs

# dummy run for reference date
model = loadSPAC(input_path, input_prefix);
ref_date = Date(model.reference_date);

start_index = Dates.value(start_date - ref_date);
end_index = Dates.value(end_date - ref_date);

# initialize metrics data frame
metrics = DataFrame(
    scen = 1:nsets,
    nse10 = fill(0.0, nsets),
    rmse10 = fill(0.0, nsets),
    nse40 = fill(0.0, nsets),
    rmse40 = fill(0.0, nsets),
    nse60 = fill(0.0, nsets),
    rmse60 = fill(0.0, nsets),
    nse80 = fill(0.0, nsets),
    rmse80 = fill(0.0, nsets)
);

## run LWFBrook90 for each parameter set
# using parallel processing

Threads.@threads for i in 1:nsets
#for i in 1:nsets
    cal_dir = output_path * subdir_name * string(i) * "/";

    if "BETAROOT" ∈ keys(root_dict) || "MAXROOTDEPTH" ∈ keys(root_dict)
        # retrieve root parameter values
        betaroot = param_sets[root_dict["BETAROOT"], i];
        maxroot= param_sets[root_dict["MAXROOTDEPTH"], i];

        # run model with modified root distribution
        model = loadSPAC(cal_dir, output_prefix, simulate_isotopes = true,
        Δz_thickness_m = "soil_discretization.csv",
        root_distribution = (beta = betaroot, z_rootMax_m=maxroot),
        IC_soil = (PSIM_init_kPa = -6.5,
        delta18O_init_permil = -13.0,
        delta2H_init_permil = -95.0));

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
        continue; # skip if simulation fails
    end

    ## retrieve model output

    # soil water content
    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = 1:end_index);

    # add dates column
    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(1:end_index,sim.parametrizedSPAC.reference_date);
    z.dates = Date.(dates_to_read_out);

    # calculate goodness of fit to observed data
    metrics[i, 2:9] = obj_fun(z, obs_swc);

end

CSV.write("LWFBcal_output/metrics_" * string(Dates.format(now(), "yyyymmdd")) * ".csv", metrics);
