# solve water balance error

using LWFBrook90
using CSV, DataFrames
using Plots, Measures; gr();

input_path = "./LWFBinput2/";
input_prefix = "pfynwald";

model_root_test = loadSPAC(input_path, input_prefix, simulate_isotopes = true,
    Δz_thickness_m = "soil_discretization.csv",
    root_distribution = (beta = 0.95, z_rootMax_m=-1.0),
    IC_soil = (PSIM_init_kPa = -6.5,
    delta18O_init_permil = -13.0,
    delta2H_init_permil = -95.0));

sim_root_test = setup(model_root_test, requested_tspan=(0, 365));
simulate!(sim_root_test, save_everystep = false, saveat = 0:365)

# water partitioning
sim_water_part_daily, sim_water_part_monthly, wp_y = get_water_partitioning(sim_root_test)
describe(sim_water_part_daily)


sim_states = get_states(sim_root_test)
describe(sim_states)
@show propertynames(sim_states)

sim_fluxes = get_fluxes(sim_root_test)
describe(sim_fluxes)
names(sim_fluxes)
@show propertynames(sim_fluxes)

pl2 = plotisotopes(sim_root_test, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
pl2

pl1 = plotamounts(sim_root_test, :above_and_belowground, :showRWUcentroid)
pl1

pl3 = plotforcingandstates(sim_root_test)
pl3


# For water balance errors
# WB error defined by Ireson et al. 2023
solution = sim_root_test.ODESolution
solu = solution

# Compute In/Out Fluxes and Storages
days = solution.t
dat_daily_cumulativeFluxes_and_Storage = DataFrame(
    # NOTE: we need not to take differences as these cumulative fluxes were set to 0
    #       every day by a solver callback
    slfl  = [solu(t).accum.slfl for t in days],
    byfl  = [solu(t).accum.byfl for t in days],
    vrfln = [solu(t).accum.vrfln for t in days],
    dsfl  = [solu(t).accum.dsfl for t in days],
    swat  = vcat([sum(solu(t).SWATI.mm, dims=1) for t in days]...), # storage
    tran  = [solu(t).accum.cum_d_tran for t in days],
    slvp  = [solu(t).accum.cum_d_slvp for t in days],
    prec  = [solu(t).accum.cum_d_prec for t in days],
    evap  = [solu(t).accum.evap for t in days], # evap is sum of irvp, isvp, snvp, slvp, sum(trani)
    flow  = [solu(t).accum.flow for t in days], # flow is sum of byfli, dsfli, gwfl
    seep  = [solu(t).accum.seep for t in days],
    gwat  = [solu(t).GWAT.mm for t in days], # storage
    snow  = [solu(t).SNOW.mm for t in days], # storage
    intr  = [solu(t).INTR.mm for t in days], # storage
    ints  = [solu(t).INTS.mm for t in days], # storage
)

describe(dat_daily_cumulativeFluxes_and_Storage)
df = dat_daily_cumulativeFluxes_and_Storage

function compute_error_SWATI(df) #
    cumInflow_mm  = cumsum(df.slfl - df.byfl) # =INFL,               corresponds to q(t,0)  in Ireson 2023 eq 16 with additionally sources and sinks
    cumOutflow_mm = cumsum(df.vrfln + df.dsfl + df.tran + df.slvp) # corresponds to q(t,zN) in Ireson 2023 eq 16 with additionally sources and sinks
    error_mm = (cumInflow_mm .- cumOutflow_mm) .- (df.swat .- df.swat[1])
    return cumInflow_mm, cumOutflow_mm, df.swat, error_mm
end
        
function compute_error_ModelDomain(df) #
    cumInflow_mm  = cumsum(df.prec)                     # corresponds to q(t,0)  in Ireson 2023 eq 16 with additionally sources and sinks
    cumOutflow_mm = cumsum(df.evap + df.flow + df.seep) # corresponds to q(t,zN) in Ireson 2023 eq 16 with additionally sources and sinks
    storage = df.swat .+ df.gwat + df.snow + df.intr + df.ints
    error_mm = (cumInflow_mm .- cumOutflow_mm) .- (storage .- storage[1])
    return cumInflow_mm, cumOutflow_mm, storage, error_mm
end

# Compute Water balance error
# a) Cumulative error metric (time evolution of water balance error)
cum_qIn_SWATI, cum_qOut_SWATI, storage_SWATI, εB_cumulative_SWATI = compute_error_SWATI(dat_daily_cumulativeFluxes_and_Storage)
cum_qIn_ALL, cum_qOut_ALL, storage_ALL, εB_cumulative_ALL = compute_error_ModelDomain(dat_daily_cumulativeFluxes_and_Storage)

# b) Scalar error metric
εB_SWATI = εB_cumulative_SWATI[end]
εB_ALL   = εB_cumulative_ALL[end]
εR_SWATI = (diff(cum_qIn_SWATI) .- diff(cum_qOut_SWATI)) .- (diff(storage_SWATI)) # mm, Bias error for intervals t = (0,tM), Ireson 2023 eq. 15
εR_SWATI = sqrt(sum(εR_SWATI.^2)/length(εR_SWATI))
εR_ALL   = (diff(cum_qIn_ALL)   .- diff(cum_qOut_ALL))   .- (diff(storage_ALL))   # mm, Bias error for intervals t = (0,tM), Ireson 2023 eq. 15
εR_ALL = sqrt(sum(εR_ALL.^2)/length(εR_ALL))
