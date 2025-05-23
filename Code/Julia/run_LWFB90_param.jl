## function to create input files for single LWFBrook90.jl run
## for a given set of parameters

using LWFBrook90, CSV, DataFrames, Dates;

function run_LWFB90_param(par, start_date, end_date, input_path, input_prefix, output_path; new_folder=true)
    # new_folder tells the function whether it needs to create a new folder
    # or if it can skip this step and just run the model with par
    
    if new_folder
        ## make output folder structure and create parameter files

        # check for soil parameters
        if "ths" ∈ names(par) || "ksat" ∈ names(par)
            soil_pars = true;
        end

        # initially assume root parameters are not present
        root_pars = false;

        # input parameter file

        param_file = input_path * input_prefix * "_param.csv";
        param0 = CSV.read(param_file, DataFrame);

        if soil_pars
            soil_file = input_path * input_prefix * "_soil_horizons.csv";
            soil0 = CSV.read(soil_file, DataFrame, skipto=3);
        end

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

        # loop through parameters
        for name in names(par)
            # get parameter value
            value = par[name];

            if name == "ths"
                # apply multiplier to ths_volfrac for each soil horizon
                soil0.ths_volFrac = soil0.ths_volFrac * value;

            elseif name == "ksat"
                # apply additive factor to log10(ksat) for each soil horizon
                soil0.ksat_mmDay = 10 .^ (log10.(soil0.ksat_mmDay) .+ value);

            elseif name ∈ ["BETAROOT", "MAXROOTDEPTH"]
                # root parameters are present
                root_pars = true;

            else
                # get index of parameter name in file
                idx = findall(param0.param_id .== name)[1];
                param0.x[idx] = value;
            end
        end

        # copy folder structure to output folder
        cp(input_path, output_path, force=true);

        # create output file name
        output_prefix = input_prefix;
        output_param_file = output_path * output_prefix * "_param.csv";

        # write parameter and soil horizons file
        CSV.write(output_param_file, param0);
        if soil_pars
            output_soil_file(soil0, output_path * output_prefix);
        end

    else
        # if new_folder is false, just define the output prefix
        output_prefix = input_prefix;

        # and check for root parameters
        if "BETAROOT" ∈ names(par) || "MAXROOTDEPTH" ∈ names(par)
            root_pars = true;
        else
            root_pars = false;
        end

    end

    ## set up model run

    # dummy run for reference date
    model = loadSPAC(output_path, output_prefix);
    ref_date = Date(model.reference_date);

    start_index = Dates.value(start_date - ref_date);
    end_index = Dates.value(end_date - ref_date);

    if root_pars
        # retrieve root parameter values
        betaroot = par["BETAROOT"];
        maxroot= par["MAXROOTDEPTH"];

        # run model with modified root distribution
        model = loadSPAC(output_path, output_prefix, simulate_isotopes = true,
        Δz_thickness_m = "soil_discretization.csv",
        root_distribution = (beta = betaroot, z_rootMax_m=maxroot),
        IC_soil = (PSIM_init_kPa = -6.5,
        delta18O_init_permil = -13.0,
        delta2H_init_permil = -95.0));

    else
        # run model with input files
        model = loadSPAC(output_path, output_prefix, simulate_isotopes = true)
    end

    # model set up
    sim = setup(model, requested_tspan=(start_index, end_index));
    # run model
    simulate!(sim);

    return sim

end