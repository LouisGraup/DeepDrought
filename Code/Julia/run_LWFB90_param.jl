## function to create input files for single LWFBrook90.jl run
## for a given set of parameters

using LWFBrook90, CSV, DataFrames, Dates;

function run_LWFB90_param(par, start_date, end_date, input_path, input_prefix, output_path; new_folder=true, watbal=false, irrig=false)
    # new_folder tells the function whether it needs to create a new folder
    # or if it can skip this step and just run the model with par
    
    if new_folder
        ## make output folder structure and create parameter files

        # create output directory if it doesn't exist
        if !ispath(output_path)
            mkpath(output_path);
        end

        # initially assume root and soil parameters are not present
        global root_pars = false;
        global soil_pars = false;

        # check for soil parameters
        if any(contains.(names(par), "ths")) || any(contains.(names(par), "ksat"))
            soil_pars = true;
            soil_par_count = count(contains.(names(par), "ths"));
        end

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

            if contains(name, "ths")
                if soil_par_count > 1
                    # apply multiplier to ths_volfrac for specific soil horizon
                    k = parse(Int, name[end]); # extract horizon number from name
                    soil0.ths_volFrac[k] = soil0.ths_volFrac[k] * value;
                else
                    # apply multiplier to ths_volfrac for each soil horizon
                    soil0.ths_volFrac = soil0.ths_volFrac * value;
                end

            elseif contains(name, "ksat")
                if soil_par_count > 1
                    # apply additive factor to log10(ksat) for specific soil horizon
                    k = parse(Int, name[end]); # extract horizon number from name
                    soil0.ksat_mmDay[k] = 10 .^ (log10.(soil0.ksat_mmDay[k]) .+ value);
                else
                    # apply additive factor to log10(ksat) for each soil horizon
                    soil0.ksat_mmDay = 10 .^ (log10.(soil0.ksat_mmDay) .+ value);
                end
            elseif contains(name, "alpha")
                if soil_par_count > 1
                    # apply multiplier to alpha_perMeter for specific soil horizon
                    k = parse(Int, name[end]); # extract horizon number from name
                    soil0.alpha_perMeter[k] = soil0.alpha_perMeter[k] * value;
                else
                    # apply multiplier to alpha_perMeter for each soil horizon
                    soil0.alpha_perMeter = soil0.alpha_perMeter * value;
                end
        
            elseif contains(name, "npar")
                if soil_par_count > 1
                    # apply multiplier to npar_ for specific soil horizon
                    k = parse(Int, name[end]); # extract horizon number from name
                    soil0.npar_[k] = soil0.npar_[k] * value;
                else
                    # apply multiplier to npar_ for each soil horizon
                    soil0.npar_ = soil0.npar_ * value;
                end
        
            elseif name ∈ ["BETAROOT", "MAXROOTDEPTH"]
                # root parameters are present
                global root_pars = true;

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
            global root_pars = true;
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
        betaroot = par.BETAROOT[1];
        maxroot= par.MAXROOTDEPTH[1];

        # run model with modified root distribution
        model = loadSPAC(output_path, output_prefix, 
        simulate_isotopes = true, simulate_irrigation = irrig,
        root_distribution = (beta = betaroot, z_rootMax_m=maxroot));

    else
        # run model with input files
        model = loadSPAC(output_path, output_prefix, 
        simulate_isotopes = true, simulate_irrigation = irrig)
    end

    # model set up
    sim = setup(model, requested_tspan=(start_index, end_index));

    if watbal
        # run model with modified output for partitioning
        simulate!(sim, save_everystep = false, saveat = start_index:end_index);
    else
        # run model as normal
        simulate!(sim);
    end

    return sim

end