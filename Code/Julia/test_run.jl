# test script to run LWFBrook90

using LWFBrook90

input_path = "./LWFBinput/";
input_prefix = "pfynwald";

model = loadSPAC(input_path, input_prefix; simulate_isotopes = true);
model

simulation = setup(model)
#simulate!(simulation)

sim_test = remakeSPAC(simulation, requested_tspan=(0, 365));
simulate!(sim_test)

model_root_test = loadSPAC(input_path, input_prefix, simulate_isotopes = true,
    Δz_thickness_m = "soil_discretization.csv",
    root_distribution = (beta = 0.95, z_rootMax_m=-1.0),
    IC_soil = (PSIM_init_kPa = -6.5,
    delta18O_init_permil = -13.0,
    delta2H_init_permil = -95.0));

sim_root_test = setup(model_root_test, requested_tspan=(0, 365));
simulate!(sim_root_test)

using Plots, Measures; gr();
pl2 = plotisotopes(sim_root_test, :d18O, (d18O = :auto, d2H = :auto), :showRWUcentroid)
pl2

pl1 = plotamounts(sim_root_test, :above_and_belowground, :showRWUcentroid)
pl1

pl3 = plotforcingandstates(sim_root_test)
pl3

# ODE manipulation

solution = sim_root_test.ODESolution;
t = solution.t; # this is an array
length(t)
j = findall(t.==1); # element-wise comparison with .==
solution[j].t == t[j]
# but cannot do
t(1)

# array interface is easier than working with vectors
# because vector notation accesses index 1 for all components (t, u, etc.)
solution[1]
# whereas array notation accesses u at t = 1
solution(1)
# so these return the same elements but the latter is cleaner
solution[j].u == solution(1)

# can't see any good reason to access u directly
u = solution.u;
# since array notation no longer works
u(1)
# but vector notation still works
u[1].aux

# even still, odd referencing, can only access singular values at a time
days = range(extrema(t)...);
evap = solution(days).accum.evap
evap = solution(1).accum.evap
# alternatively,
days = range(solution.prob.tspan...),

# if you knew the ID, could do it this way
evap = solution(days, idxs=155).u
# or
days_j = findall(t .∈ [[days...]] != 0)
evap = solution[155, days_j]

# so explains Fabian's solution
evap  = [solution(t).accum.evap for t in days]

plot(sim_root_test.ODESolution;
    idxs = 43:3:62,
    label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "RWU (mm/day)" "XYLEM (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])


# absolutely easier to just use built-in functions

z = get_soil_(:theta, sim_root_test, depths_to_read_out_mm = [100, 400, 600], days_to_read_out_d = days)

dates_to_read_out = LWFBrook90.RelativeDaysFloat2DateTime.(
    days,sim_root_test.parametrizedSPAC.reference_date)

z.dates = dates_to_read_out