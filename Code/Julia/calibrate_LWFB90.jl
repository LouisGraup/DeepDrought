# calibration and sensitivity analysis of LWFBrook90.jl
# using soil water data as behavioral constraints

using LWFBrook90
using CSV, DataFrames
using Random, Distributions
using Distributed
using Plots; gr()

# objective function to compare model output to observed data
function obj_fun(sim, obs)
    
    function NSE(sim, obs)
        # calculate Nash-Sutcliffe Efficiency
        nse = 1 - (sum((obs - sim)^2) / sum((obs - mean(obs))^2))
        return nse
    end
    
    function RMSE(sim, obs)
        # calculate Root Mean Square Error
        rmse = sqrt(sum((obs - sim)^2) / length(obs))
        return rmse
    end

    # calculate NSE
    nse = NSE(sim, obs)
    # calculate RMSE
    rmse = RMSE(sim, obs)
    

    return nse
end

# behavioral data
obs_swc = CSV.read("./Data/Pfyn/PFY_swat.csv", DataFrame);

## parameter input and output paths
# input
input_path = "./LWFBinput/";
input_prefix = "pfynwald";

# output
output_path = "./LWFBcalibration/";
subdir_name = "cal";

## define calibration parameter sets

n = 100; # number of parameter sets

# define prior parameter ranges

param = [
    # hydro parameters
    ("DRAIN", rand(Uniform(0,1), n)), # drainage
    ("BYPAR", rand(Uniform(0,1), n)), # bypass flow
    # soil parameters
    ("RSSA", rand(Uniform(1,1000), n)), # soil resistance
    #("ths", rand(Uniform(.5,2), n)), # multiplier on theta_sat
    #("ksat", rand(Uniform(.5,2), n)), # multiplier on k_sat
    # plant parameters
    ("GLMAX", rand(Uniform(0.001,0.03), n)), # stomatal conductance
    ("CVPD", rand(Uniform(1,3), n)), # vpd sensitivity
    ("PSICR", rand(Uniform(-4,-1), n)), # critical water potential
    ("FXYLEM", rand(Uniform(0.2,0.8), n)), # aboveground xylem fraction
    ("MXKPL", rand(Uniform(5,30), n)), # maximum plant conductivity
    ("VXYLEM_mm", rand(Uniform(1,100), n)), # xylem volume
    ("DISPERSIVITY_mm", rand(Uniform(1,100), n)), # dispersivity coefficient
    #("MAXROOTDEPTH", rand(Uniform(-5,.5), n)), # max rooting depth
    #("BETAROOT", rand(Uniform(.8,1), n)) # beta root coefficient
];

# expand parameter sets to LHS or Sobol


## make output folder structure and create calibration parameter files

# create output directory if it doesn't exist
if !isdir(output_path)
    mkdir(output_path);
end

output_prefix = input_prefix;

# input parameter file

param_file = input_path * input_prefix * "_param.csv";
param0 = CSV.read(param_file, DataFrame);

# loop through parameter sets and create parameter files
for i in 1:n
    # create parameter set
    param_set = param0;
    for (name, value) in param
        # get index of parameter name in file
        idx = findall(param0.param_id .== name)[1];
        param_set.x[idx] = value[i];
    end

    # output folder
    out_dir = output_path * subdir_name * string(i) * "/";
    # copy folder structure to output folder
    cp(input_path, out_dir, force=true);

    # create output file name
    output_file = out_dir * output_prefix * "_param" * ".csv";

    # write parameter file
    CSV.write(output_file, param_set);
end

# for checking individual parameter values
# for i in 1:length(param)
#     println("Parameter $(param[i][1]) = $(param[i][2][1])")
# end

# run LWFBrook90 for each parameter set
# using parallel processing

#@threads for i in 1:n
for i in 1:n
    cal_dir = output_path * subdir_name * string(i) * "/";

    model = loadSPAC(cal_dir, output_prefix, simulate_isotopes = true,
    Δz_thickness_m = "soil_discretization.csv",
    root_distribution = (beta = 0.95, z_rootMax_m=-1.0),
    IC_soil = (PSIM_init_kPa = -6.5,
    delta18O_init_permil = -13.0,
    delta2H_init_permil = -95.0));

    sim = setup(model, requested_tspan=(0, 365));
    simulate!(sim)

    # retrieve model output

    # time span
    days = range(1, max(sim.ODESolution.t));

    # soil water content
    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 400, 600], days_to_read_out_d = days)

    dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim_root_test.parametrizedSPAC.reference_date)

    z.dates = dates_to_read_out

    # calculate goodness of fit to observed data

end