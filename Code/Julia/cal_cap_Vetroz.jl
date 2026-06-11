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
    # dendrometer data
    obs_twd = CSV.read("../../Data/Daten_Vetroz_Lorenz/TWD_daily.csv", DataFrame);
    obs_twd = obs_twd[obs_twd.sensor_loc .== "stem", :]; # filter to stem sensors
    obs_twd = obs_twd[obs_twd.year .< 2023, :]; # filter out last years
    select!(obs_twd, Not([:month, :year, :sensor_loc, :MDS])); # drop unnecessary columns
end

# objective function to compare model output to observed data
@everywhere function obs_fun_twd(twd_comp)

    # split into pd and md
    twd_comp_pd = twd_comp[:, [:date, :TWD_pd, :cum_pd_plpsi]];
    lwp_comp_pd = twd_comp[:, [:date, :LWP_pd, :cum_pd_plpsi]];
    twd_comp_md = twd_comp[:, [:date, :TWD_md, :cum_md_plpsi]];
    lwp_comp_md = twd_comp[:, [:date, :LWP_md, :cum_md_plpsi]];

    # remove missing values
    twd_comp_pd = dropmissing(twd_comp_pd);
    lwp_comp_pd = dropmissing(lwp_comp_pd);
    twd_comp_md = dropmissing(twd_comp_md);
    lwp_comp_md = dropmissing(lwp_comp_md);

    cc_twd_pd = cor(twd_comp_pd.cum_pd_plpsi, twd_comp_pd.TWD_pd);
    cc_lwp_pd = cor(lwp_comp_pd.cum_pd_plpsi, lwp_comp_pd.LWP_pd);
    cc_twd_md = cor(twd_comp_md.cum_md_plpsi, twd_comp_md.TWD_md);
    cc_lwp_md = cor(lwp_comp_md.cum_md_plpsi, lwp_comp_md.LWP_md);

    return cc_twd_pd, cc_lwp_pd, cc_twd_md, cc_lwp_md

end

@everywhere function twd_combine(sim, obs_twd)
    # sim is the LWFBrook90 simulation
    # obs is the observed dendrometer data

    z_plpsi = get_plpsi(sim);

    twd_comp = sort(innerjoin(obs_twd, z_plpsi, on = :date), :date);

    return twd_comp
end

@everywhere function get_plpsi(sim)
    # retrieve plant water potential from sim

    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    select!(z, :date, :cum_pd_plpsi, :cum_md_plpsi);

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
    # plant parameters
    ("CAPACITANCE", 0.1, 10.0), # capacitance (0.1, 10)
    ("STORAGEK", 0.1, 10.0), # storage conductance (0.1, 10)
    ("VSTORAGE", 1.0, 20.0), # stem storage volume (1, 20)
    ("PSICR", -4.0, -1.7) # critical water potential (-4, -1)
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
        return (par_id, fill(0, 4)) # skip if simulation fails
    end

    ## retrieve model output

    # plant potential
    twd_comp = twd_combine(sim, obs_twd);
    twd_cor = obs_fun_twd(twd_comp);
    
    return (par_id, twd_cor)

end


# parallel map for calibration runs
results = pmap(i -> run_calibration(i), 1:nsets);

# intialize metric dataframes
metrics = DataFrame(scen = Int[],
    twd_pd_cor = Float64[],
    lwp_pd_cor = Float64[],
    twd_md_cor = Float64[],
    lwp_md_cor = Float64[]);

# loop through results
for res in results
    # retrieve scenario, parameter id, and metrics
    par_id, twd = res;
    row = [par_id, twd...];
    push!(metrics, row);
end

CSV.write("LWFBcal_output/metrics_vetroz_" * curDate * ".csv", metrics);
