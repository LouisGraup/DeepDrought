using Distributed

addprocs(5)

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
obs_swc = CSV.read("../../Data/Pfyn/PFY_swat.csv", DataFrame);
obs_swc.VWC = obs_swc.VWC / 100; # convert to decimal
filter!(:date => >=(Date(2004, 1, 1)), obs_swc); # filter out early dates

obs_swp = CSV.read("../../Data/Pfyn/PFY_swpc.csv", DataFrame);
filter!(:date => >=(Date(2015, 1, 1)), obs_swp); # filter out early dates
filter!(:date => <(Date(2021, 1, 1)), obs_swp); # filter out late dates

# separate control and irrigation scenarios
obs_swc_ctr = obs_swc[obs_swc.meta .== "control", :]; # select control treatment
select!(obs_swc_ctr, :date, :depth, :VWC); # remove extra columns
obs_swc_ctr = unstack(obs_swc_ctr, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
sort!(obs_swc_ctr, :date); # sort by date

obs_swc_irr = obs_swc[obs_swc.meta .== "irrigated", :]; # select irrigation treatment
select!(obs_swc_irr, :date, :depth, :VWC); # remove extra columns
obs_swc_irr = unstack(obs_swc_irr, :date, :depth, :VWC, renamecols=x->Symbol("VWC_$(x)cm")); # reshape data
sort!(obs_swc_irr, :date); # sort by date

obs_swp_ctr = obs_swp[obs_swp.meta .== "control", :]; # select control treatment
obs_swp_irr = obs_swp[obs_swp.meta .== "irrigated", :]; # select control treatment

end

# objective function to compare model output to observed data
@everywhere function obj_fun_swc(sim, obs)

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

### BEGIN USER INPUT ###

@everywhere begin
    ## parameter input and output paths
    # input
    input_path_ctr = "LWFBinput/Pfyn_control/";
    input_path_irr = "LWFBinput/Pfyn_irrigation_ambient/";
    input_prefix = "pfynwald";

    # output
    output_path = "LWFBcalibration/";
    subdir_name_ctr = "cal_ctr";
    subdir_name_irr = "cal_irr";

    ## simulation dates

    start_date = Date(2000, 1, 1);
    end_date = Date(2020, 12, 31);


    nsets = 20;

    output_prefix = input_prefix;

    param_in = CSV.read("LWFBcal_output/param.csv", DataFrame);
    param_sets = Matrix(param_in);
    param_sets = param_sets';

    root_dict = Dict("ROOTS" => 0);
    root_params = true; # boolean to check for root parameters
    push!(root_dict, "MAXROOTDEPTH" => 19);
    push!(root_dict, "BETAROOT" => 20);

end

model = loadSPAC(input_path_ctr, input_prefix);
ref_date = Date(model.reference_date);

start_index = Dates.value(start_date - ref_date);
end_index = Dates.value(end_date - ref_date);

# metrics_ctr = DataFrame(
#     scen = 1:nsets,
#     swc_nse10 = fill(0.0, nsets),
#     swc_rmse10 = fill(0.0, nsets),
#     swc_nse40 = fill(0.0, nsets),
#     swc_rmse40 = fill(0.0, nsets),
#     swc_nse60 = fill(0.0, nsets),
#     swc_rmse60 = fill(0.0, nsets),
#     swc_nse80 = fill(0.0, nsets),
#     swc_rmse80 = fill(0.0, nsets),
#     swp_nse10 = fill(0.0, nsets),
#     swp_nse80 = fill(0.0, nsets)
# );


#@everywhere obs_swc_ctr = deepcopy(Main.obs_swc_ctr)
#@everywhere obs_swc_irr = deepcopy(Main.obs_swc_irr)
#@everywhere obs_swp_ctr = deepcopy(Main.obs_swp_ctr)
#@everywhere obs_swp_irr = deepcopy(Main.obs_swp_irr)
#@everywhere param_sets = deepcopy(Main.param_sets)
#@everywhere root_dict = deepcopy(Main.root_dict)
#@everywhere output_path = Main.output_path
#@everywhere subdir_name_ctr = Main.subdir_name_ctr
#@everywhere subdir_name_irr = Main.subdir_name_irr
#@everywhere output_prefix = Main.output_prefix
#@everywhere root_params = Main.root_params
@everywhere start_index = $start_index
@everywhere end_index = $end_index
#@everywhere nsets = Main.nsets


@everywhere function run_calibration(i)
    if i <= nsets
        cal_dir = output_path * subdir_name_ctr * string(i) * "/";
        par_id = i
        obs_swc = obs_swc_ctr
        obs_swp = obs_swp_ctr
    else
        cal_dir = output_path * subdir_name_irr * string(i - nsets) * "/";
        par_id = i - nsets
        obs_swc = obs_swc_irr
        obs_swp = obs_swp_irr
    end

    if root_params
        # retrieve root parameter values
        betaroot = param_sets[root_dict["BETAROOT"], par_id];
        maxroot = param_sets[root_dict["MAXROOTDEPTH"], par_id];

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
        return (i <= nsets ? "ctr" : "irr", par_id, fill(0, 8), fill(0,2)) # skip if simulation fails
    end

    ## retrieve model output

    # soil water content
    z_theta = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600, 800], days_to_read_out_d = 1:end_index);
    z_psi = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800], days_to_read_out_d = 1:end_index);

    # add dates column
    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(1:end_index,sim.parametrizedSPAC.reference_date);
    z_theta.dates = Date.(dates_to_read_out);
    z_psi.dates = Date.(dates_to_read_out);

    swc_met = obj_fun_swc(z_theta, obs_swc)
    swp_met = obj_fun_swp(z_psi, obs_swp)

    return (i <= nsets ? "ctr" : "irr", par_id, swc_met, swp_met)

end


# parallel map for calibration runs
results = pmap(i -> run_calibration(i), 1:(2*nsets));

# intialize metric dataframes
metrics_ctr = DataFrame(scen = Int[],
    swc_nse10 = Float64[],
    swc_rmse10 = Float64[],
    swc_nse40 = Float64[],
    swc_rmse40 = Float64[],
    swc_nse60 = Float64[],
    swc_rmse60 = Float64[],
    swc_nse80 = Float64[],
    swc_rmse80 = Float64[],
    swp_nse10 = Float64[],
    swp_nse80 = Float64[])
metrics_irr = copy(metrics_ctr);

# loop through results
for res in results
    # retrieve scenario, parameter id, and metrics
    scenario, par_id, swc, swp = res;
    row = [par_id, swc..., swp...];
    if scenario == "ctr"
        push!(metrics_ctr, row);
    else
        push!(metrics_irr, row);
    end
end

CSV.write("LWFBcal_output/metrics_ctr_par.csv", metrics_ctr);
CSV.write("LWFBcal_output/metrics_irr_par.csv", metrics_ctr);
