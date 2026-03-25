# compare LWFBrook90.jl with and without capacitance
# first add pull request from Github
# add https://github.com/fabern/LWFBrook90.jl#add-below-canopy-irrigation
# then switch to new model
# add https://github.com/fabern/LWFBrook90.jl#capacitance-implementation

using CSV, Dates, DataFrames
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Measures, Plots; gr()
using LWFBrook90

# functions to derive effective soil water potential
function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function get_RWU_centroid(sim)
    # borrow code from LWFBrook90 package
    solu = sim.ODESolution;

    days_to_read_out_d = unique(round.(solu.t));

    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2;

    # Compute RWU centroid
    rows_RWU_mmDay  = reduce(hcat, [solu(t).TRANI.mmday   for t in days_to_read_out_d]);

    RWU_percent = rows_RWU_mmDay ./ sum(rows_RWU_mmDay; dims = 1);
    #RWUcentroidLabel = "mean RWU depth"
    if (any(RWU_percent .< 0))
        #@warn "Some root water outfluxes detected. Centroid of RWU is  based only on uptakes."
        rows_RWU_mmDay_onlyUptake = ifelse.(rows_RWU_mmDay.>0,rows_RWU_mmDay, 0);
        RWU_percent_onlyUptake = rows_RWU_mmDay_onlyUptake ./ sum(rows_RWU_mmDay_onlyUptake; dims = 1);
        RWU_percent = RWU_percent_onlyUptake;
        
        #RWUcentroidLabel = "mean RWU depth\n(based on uptake only)"
    end

    row_RWU_centroid_mm = sum(RWU_percent .* y_center; dims=1);

    col_RWU_centroid_mm = reshape(row_RWU_centroid_mm, :);
    
    return col_RWU_centroid_mm, RWU_percent
end

function get_eff_swp(sim)
    # derive effective soil water potential based on RWU
    days, dates_out = get_dates(sim);
    swp = DataFrame(date = dates_out);
    
    swp.RWU, rwu_per = get_RWU_centroid(sim); # RWU depth and percent

    swp_all = get_soil_(:psi, sim, days_to_read_out_d=days); # swp

    swp.swp_eff .= sum(rwu_per .* Matrix(swp_all[:, Not(:time)])', dims=1)';

    return swp
end

# leaf water potential data
obs_twd_lwp = CSV.read("../../Data/Pfyn/Pfyn_TWD_LWP_2024.csv", DataFrame);
obs_lwp = CSV.read("../../Data/Pfyn/Pfyn_LWP.csv", DataFrame);
obs_lwp = obs_lwp[obs_lwp.date .>= Date(2024, 1, 1) .&& obs_lwp.date .< Date(2025, 1, 1) .&& obs_lwp.meta .== "control", :];

input_prefix = "pfynwald";
#input_path   = "LWFB_testcap/control/";
input_path   = "../../../LWFBrook90.jl/examples/PFY2024-capacitance/";

model = LWFBrook90.loadSPAC(input_path, input_prefix; simulate_isotopes = false,
    root_distribution = (beta = 0.97091, z_rootMax_m=-1.35462));

simulation = LWFBrook90.setup(model); #, requested_tspan=(8766, 9131));
LWFBrook90.simulate!(simulation)

#simulation = remakeSPAC(model; params = (STORAGEK = 1,)); #remakeSPAC returns setup
#LWFBrook90.simulate!(simulation);

# extract output
sol = simulation.ODESolution;
timesteps = unique(round.(sol.t));
dates = LWFBrook90.RelativeDaysFloat2DateTime.(timesteps, simulation.parametrizedSPAC.reference_date);

PLSTOR = [sol(t).PLSTOR.mm for t in timesteps];
PLVWC = [sol(t).PLHYD.θ for t in timesteps];
PLPSI = [sol(t).PLHYD.ψ for t in timesteps];
PLPSI_pd = [sol(t).accum.cum_pd_plpsi for t in timesteps];
PLPSI_md = [sol(t).accum.cum_md_plpsi for t in timesteps];
PLFL = [sol(t).accum.cum_d_plfl for t in timesteps];
PLRF = [sol(t).accum.cum_d_plrf for t in timesteps];
RWU = [sol(t).RWU.mmday for t in timesteps];
TRANI = reduce(hcat, [sol(t).TRANI.mmday for t in timesteps]);
PLRFI = reduce(hcat, [sol(t).PLRFI.mmday for t in timesteps]);

soil = LWFBrook90.get_soil_(:ψ, simulation);
flux = LWFBrook90.get_fluxes(simulation);
flux.date = Date.(dates);
flux[153:158, :]


Plots.plot(PLPSI_pd, label="Pre-dawn ψ")
Plots.plot(dates[153:158], PLPSI_pd[153:158], label="Pre-dawn ψ")
Plots.plot(122:303, PLPSI_pd[122:303], label="Pre-dawn ψ")
Plots.plot(dates, RWU, label="RWU")
Plots.plot(dates, PLFL, label="PLFL")
Plots.plot(dates, PLFL./(RWU.+PLFL), label="PLFL/RWU")

psi_comp = DataFrame(t=dates, PLPSI=PLPSI_pd)
psi_comp = [psi_comp; DataFrame(t=dates.+Hour(12), PLPSI=PLPSI_md)]
sort!(psi_comp, :t)
Plots.plot(psi_comp.t[305:315], psi_comp.PLPSI[305:315])
Plots.plot(psi_comp.t, psi_comp.PLPSI)

Plots.plot(PLPSI_pd, label="Pre-dawn ψ")
Plots.plot!(PLPSI_md, label="Midday ψ")

# compare plant water potential with all soil layers
soil = LWFBrook90.get_soil_(:ψ, simulation, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100, 200, 400, 600, 800, 1000, 1200]);
soil.date = Date.(dates);
soil_long = stack(soil, Not([:time, :date]), variable_name=:layer, value_name=:swp);
soil_long.layer = parse.(Int, replace.(soil_long.layer, r"[^\d]" => ""));

draw(
    data(soil_long[soil_long.date .>= Date(2024, 5, 1) .&& soil_long.date .<= Date(2024, 11, 1), :])*mapping(:date, :swp, color=:layer => nonnumeric => "Soil layer (mm)")*visual(Lines)+
    data(flux[flux.date .>= Date(2024, 5, 1) .&& flux.date .<= Date(2024, 11, 1), :])*mapping(:date, :cum_pd_plpsi)*visual(Lines, color=:red, label="Plant ψ")+
    data(obs_twd_lwp)*mapping(:date, :LWP_pd => (x -> x * 1000))*visual(Scatter, color=:black, label="TWD-derived LWP"),
    scales(X = (; label=""), Y = (; label="Water potential (kPa)"), Color= (; palette=from_continuous(:summer))),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# compare plant water potential with effective soil water potential
swp_eff = get_eff_swp(simulation);
swp_eff.swp_eff[swp_eff.swp_eff .< -2000] .= -2000;

# pre-dawn
draw(
    data(swp_eff[swp_eff.date .>= Date(2024, 5, 1) .&& swp_eff.date .<= Date(2024, 11, 1), :])*mapping(:date, :swp_eff)*visual(Lines, color=:blue, label="Effective SWP")+
    data(flux[flux.date .>= Date(2024, 5, 1) .&& flux.date .<= Date(2024, 11, 1), :])*mapping(:date, :cum_pd_plpsi)*visual(Lines, color=:red, label="Plant ψ")+
    data(obs_twd_lwp)*mapping(:date, :LWP_pd => (x -> x * 1000))*visual(Scatter, color=:black, label="TWD-derived LWP")+
    data(obs_twd_lwp)*mapping(:date, :twd_pdn => (x -> x * -500))*visual(Scatter, color=:green, label="Scaled TWD"),
    #data(obs_lwp)*mapping(:date, :predawn_mean => (x -> x * 1000))*visual(Scatter, color=:green, label="Measured LWP"),
    scales(X = (; label=""), Y = (; label="Water potential (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# midday
draw(
    data(swp_eff[swp_eff.date .>= Date(2024, 5, 1) .&& swp_eff.date .<= Date(2024, 11, 1), :])*mapping(:date, :swp_eff)*visual(Lines, color=:blue, label="Effective SWP")+
    data(flux[flux.date .>= Date(2024, 5, 1) .&& flux.date .<= Date(2024, 11, 1), :])*mapping(:date, :cum_md_plpsi)*visual(Lines, color=:red, label="Plant ψ")+
    data(obs_twd_lwp)*mapping(:date, :LWP_md => (x -> x * 1000))*visual(Scatter, color=:black, label="TWD-derived LWP")+
    data(obs_twd_lwp)*mapping(:date, :twd_md => (x -> x * -2))*visual(Scatter, color=:green, label="Scaled TWD"),
    #data(obs_lwp)*mapping(:date, :predawn_mean => (x -> x * 1000))*visual(Scatter, color=:green, label="Measured LWP"),
    scales(X = (; label=""), Y = (; label="Water potential (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# SPAC plot
p1 = Plots.plot(Date.(dates), RWU, color=:green, xlabel="", ylabel="Transpiration (mm/d)", legend=false, size=(1000, 500), left_margin=4mm);
p2 = Plots.plot(swp_eff.date, swp_eff.swp_eff, color=:black, label="Effective SWP", ylabel="Water Potential (kPa)", xlabel="", size=(1000, 500), left_margin=4mm);
p2 = Plots.plot!(p2, Date.(dates), PLPSI_pd, color=:blue, label="Plant ψ")
Plots.plot(p1, p2, layout=(2, 1))

# investigate drops in soil water potential
draw(
    data(soil_long[soil_long.date .>= Date(2024, 8, 1) .&& soil_long.date .<= Date(2024, 9, 1), :])*mapping(:date, :swp, color=:layer => nonnumeric => "Soil layer (mm)")*visual(Lines),
    scales(X = (; label=""), Y = (; label="Water potential (kPa)"), Color= (; palette=from_continuous(:summer))),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# soil water content
soilt = LWFBrook90.get_soil_(:θ, simulation, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100, 200, 400, 600, 800, 1000, 1200]);
soilt.date = Date.(dates);
soilt_long = stack(soilt, Not([:time, :date]), variable_name=:layer, value_name=:vwc);
soilt_long.layer = parse.(Int, replace.(soilt_long.layer, r"θ_m3m3_(\d+)mm" => s"\1"));

draw(
    data(soilt_long[soilt_long.date .>= Date(2024, 8, 1) .&& soilt_long.date .<= Date(2024, 9, 1), :])*mapping(:date, :vwc, color=:layer => nonnumeric => "Soil layer (mm)")*visual(Lines),
    scales(X = (; label=""), Y = (; label="VWC")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# and related conductivity
soilk = LWFBrook90.get_soil_(:K, simulation, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100, 200, 400, 600, 800, 1000, 1200]);
soilk.date = Date.(dates);
soilk_long = stack(soilk, Not([:time, :date]), variable_name=:layer, value_name=:k);
soilk_long.layer = parse.(Int, replace.(soilk_long.layer, r"[^\d]" => ""));

draw(
    data(soilk_long[soilk_long.date .>= Date(2024, 8, 1) .&& soilk_long.date .<= Date(2024, 9, 1), :])*mapping(:date, :k => (k -> log(k)), color=:layer => nonnumeric => "Soil layer (mm)")*visual(Lines),
    scales(X = (; label=""), Y = (; label="K")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# is it root water uptake 
TRANI = get_soil_(:TRANI, simulation, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100, 200, 400, 600, 800, 1000, 1200]);
TRANI.date = Date.(dates);
TRANI_long = stack(TRANI, Not([:time, :date]), variable_name=:layer, value_name=:rwu);
TRANI_long.layer = parse.(Int, replace.(TRANI_long.layer, r"[^\d]" => ""));

draw(
    data(TRANI_long[TRANI_long.date .>= Date(2024, 8, 1) .&& TRANI_long.date .<= Date(2024, 9, 1), :])*mapping(:date, :rwu, color=:layer => nonnumeric => "Soil layer (mm)")*visual(Lines),
    scales(X = (; label=""), Y = (; label="RWU")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# or plant refill
PLRFI = get_soil_(:PLRFI, simulation, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100, 200, 400, 600, 800, 1000, 1200]);
PLRFI.date = Date.(dates);
PLRFI_long = stack(PLRFI, Not([:time, :date]), variable_name=:layer, value_name=:plrf);
PLRFI_long.layer = parse.(Int, replace.(PLRFI_long.layer, r"[^\d]" => ""));

draw(
    data(PLRFI_long[PLRFI_long.date .>= Date(2024, 8, 1) .&& PLRFI_long.date .<= Date(2024, 9, 1), :])*mapping(:date, :plrf, color=:layer => nonnumeric => "Soil layer (mm)")*visual(Lines),
    scales(X = (; label=""), Y = (; label="PLRF")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# combine variables into single plot
soil_comp = get_soil_([:ψ, :θ, :K, :TRANI, :PLRFI], simulation, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100]);
rename!(soil_comp, [:time, :K, :PLRFI, :TRANI, :θ, :ψ]);
soil_comp.date = Date.(dates);
soil_plot = soil_comp[soil_comp.date .>= Date(2024, 8, 1) .&& soil_comp.date .<= Date(2024, 9, 1), :];

draw(
    data(soil_plot)*mapping(:date, :ψ => x -> x / 1000)*visual(Lines, color=:purple, label="Soil ψ")+
    data(soil_plot)*mapping(:date, :θ => x -> x * 100)*visual(Lines, color=:blue, label="Soil θ")+
    data(soil_plot)*mapping(:date, :K => k -> log(k))*visual(Lines, color=:red, label="Soil K")+
    data(soil_plot)*mapping(:date, :TRANI => x -> x * 100)*visual(Lines, color=:green, label="RWU")+
    data(soil_plot)*mapping(:date, :PLRFI => x -> x * 100)*visual(Lines, color=:orange, label="PLRF"),
    scales(X = (; label=""), Y = (; label="Var")),
    figure = (; size=(1200, 600))
)

# automated loop over parameter values for sensitivity analysis

ks_par = [0.1, 1, 2.4, 5, 10]
scen_comp = nothing;
cap_comp = nothing;

for ks ∈ ks_par

    model = remakeSPAC(model; params = (STORAGEK = ks,)); #remakeSPAC returns setup
    LWFBrook90.simulate!(model);

    # extract output
    sol = model.ODESolution;
    timesteps = unique(round.(sol.t));
    dates = LWFBrook90.RelativeDaysFloat2DateTime.(timesteps, model.parametrizedSPAC.reference_date);

    PLSTOR = [sol(t).PLSTOR.mm for t in timesteps];
    PLPSI_pd = [sol(t).accum.cum_pd_plpsi for t in timesteps];
    PLPSI_md = [sol(t).accum.cum_md_plpsi for t in timesteps];
    PLFL = [sol(t).accum.cum_d_plfl for t in timesteps];
    PLRF = [sol(t).accum.cum_d_plrf for t in timesteps];
    RWU = [sol(t).RWU.mmday for t in timesteps];
    
    soil = LWFBrook90.get_soil_(:ψ, model, days_to_read_out_d=timesteps, depths_to_read_out_mm=[100, 200, 400, 600, 800, 1000, 1200]);
    soil.date = Date.(dates);

    swp_eff = get_eff_swp(model);

    soil.ψ_kPa_eff = swp_eff.swp_eff;
    soil.RWU = RWU + PLRF; # actual RWU
    soil.AET = RWU + PLFL; # actual transpiration
    select!(soil, Not(:time))
    soil.ks .= ks;

    if isnothing(scen_comp)
        scen_comp = soil;
    else
        scen_comp = [scen_comp; soil];
    end

    plant = DataFrame(date = Date.(dates), PLPSI_pd = PLPSI_pd, PLPSI_md = PLPSI_md, PLFL = PLFL, PLRF = PLRF, RWU = RWU, PLSTOR = PLSTOR);
    plant.ks .= ks;

    if isnothing(cap_comp)
        cap_comp = plant;
    else
        cap_comp = [cap_comp; plant];
    end

end

min_psi = combine(groupby(scen_comp, :ks), names(scen_comp, r"ψ") .=> [minimum]);
mapcols(x -> minimum(x[Not(isnan.(x))]), scen_comp[scen_comp.ks .== 2, Not([:date, :RWU, :ks])])

combine(groupby(scen_comp, :ks), :RWU => sum)
combine(groupby(scen_comp, :ks), :AET => sum)

combine(groupby(cap_comp, :ks), :PLPSI_pd .=> [minimum])

# analysis

scen_comp = CSV.read("Pfyn_cap_soil.csv", DataFrame);
cap_comp = CSV.read("Pfyn_ks_plant.csv", DataFrame);

cap_comp.PLPSI_pd[cap_comp.PLPSI_pd .< -2000] .= -2000 # cap plant psi for visualization
draw(
    data(cap_comp)*mapping(:date, :PLPSI_pd, color=:cap => nonnumeric => "Capacitance")*visual(Lines),
    scales(X = (; label=""), Y = (; label="Plant ψ (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

draw(
    data(cap_comp)*mapping(:date, :PLPSI_pd, color=:ks => nonnumeric => "Conductance")*visual(Lines),
    scales(X = (; label=""), Y = (; label="Plant ψ (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

draw(
    data(cap_comp[cap_comp.cap .!= 0.1, :])*mapping(:date, :PLPSI_md, color=:cap => nonnumeric => "Capacitance")*visual(Lines),
    scales(X = (; label=""), Y = (; label="Midday plant ψ (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

draw(
    data(cap_comp)*mapping(:date, :PLPSI_md, color=:ks => nonnumeric => "Conductance")*visual(Lines),
    scales(X = (; label=""), Y = (; label="Midday plant ψ (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

# compare with model without capacitance
prior = CSV.read("Pfyn_prior.csv", DataFrame);

mapcols(x -> minimum(x[Not(isnan.(x))]), prior[:, Not([:date, :RWU, :cap])])
sum(prior.RWU)


draw(
    data(scen_comp[scen_comp.date .>= Date(2024, 8, 1) .&& scen_comp.date .< Date(2024, 9, 10), :])*mapping(:date, :AET, color=:cap => nonnumeric => "Capacitance")*visual(Lines)+
    data(prior[prior.date .>= Date(2024, 8, 1) .&& prior.date .< Date(2024, 9, 10), :])*mapping(:date, :RWU)*visual(Lines, color=:red, label="No capacitance"),
    scales(X = (; label=""), Y = (; label="Transpiration (mm/day)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)

draw(
    data(scen_comp[scen_comp.date .>= Date(2024, 8, 1) .&& scen_comp.date .< Date(2024, 9, 1), :])*mapping(:date, :ψ_kPa_eff, color=:cap => nonnumeric => "Capacitance")*visual(Lines)+
    data(prior[prior.date .>= Date(2024, 8, 1) .&& prior.date .< Date(2024, 9, 1), :])*mapping(:date, :ψ_kPa_eff)*visual(Lines, color=:red, label="No capacitance"),
    scales(X = (; label=""), Y = (; label="Effective SWP (kPa)")),
    legend = (; position = :bottom, framevisible = false), figure = (; size=(1200, 600))
)
