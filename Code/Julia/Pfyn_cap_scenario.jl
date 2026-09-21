## run best performing scenario from LWFBrook90.jl capacitance calibration
## for all Pfynwald scenarios

using CSV, DataFrames, DataFramesMeta, Dates, Statistics, RollingFunctions;
using CairoMakie, AlgebraOfGraphics, CategoricalArrays, Chain;
using Measures, Plots; gr()

include("run_LWFB90_param.jl");

function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...)[Not(end)];
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function get_swc(sim; shape = "long")
    # retrieve soil water potential data from sim
    days, dates_out = get_dates(sim);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [200, 800, 1100, 1600], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="VWC");
        z_long.depth = parse.(Int, replace.(z_long.depth, r"theta_m3m3_(\d+)mm" => s"\1")) .÷ 10; # convert depth to cm from var name
        return z_long
    else
        return z
    end

end

function get_swp(sim; shape="long")
    # retrieve soil water potential data from sim
    days, dates_out = get_dates(sim);
    
    z = get_soil_(:psi, sim, depths_to_read_out_mm = [200, 800, 1100, 1600], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="SWP");
        z_long.depth = parse.(Int, replace.(z_long.depth, r"psi_kPa_(\d+)mm" => s"\1")) .÷ 10; # convert depth to cm from var name
        
        return z_long
    else
        return z
    end

end

function get_sap(sim)
    # retrieve transpiration from sim
    
    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    z.trans = z.cum_d_tran;
    select!(z, :date, :trans);
    
    return z
end

function ann_trans(sim)
    # calculates annual transpiration
    z = get_sap(sim);

    z.year = year.(z.date);

    z_ann = combine(groupby(z, :year), :trans => sum);

    return z_ann
end

function get_RWU_centroid(sim)
    # borrow code from LWFBrook90 package
    solu = sim.ODESolution;
    saved = sim.saved_values;

    days_to_read_out_d = saved.t;

    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2;

    # Compute RWU centroid
    rows_RWU_mmDay  = reduce(hcat, [saved.saveval[t].TRANI for t in 1:(length(days_to_read_out_d)-1)]);

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

function get_PLRF_centroid(sim)
    # borrow code from LWFBrook90 package
    solu = sim.ODESolution;
    saved = sim.saved_values;

    days_to_read_out_d = saved.t;

    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2;

    # Compute PLRF centroid
    rows_PLRF_mmDay  = reduce(hcat, [saved.saveval[t].PLRFI for t in 1:(length(days_to_read_out_d)-1)]);

    PLRF_percent = rows_PLRF_mmDay ./ sum(rows_PLRF_mmDay; dims = 1);

    if (any(PLRF_percent .< 0))
        rows_PLRF_mmDay_onlyUptake = ifelse.(rows_PLRF_mmDay.>0,rows_PLRF_mmDay, 0);
        PLRF_percent_onlyUptake = rows_PLRF_mmDay_onlyUptake ./ sum(rows_PLRF_mmDay_onlyUptake; dims = 1);
        PLRF_percent = PLRF_percent_onlyUptake;
    end

    row_PLRF_centroid_mm = sum(PLRF_percent .* y_center; dims=1);

    col_PLRF_centroid_mm = reshape(row_PLRF_centroid_mm, :);
    
    return col_PLRF_centroid_mm, PLRF_percent
end

function get_REW(sim)
    # derive relate extractable soil water based on maximum root depth
    # for now, just use normalized SWAT over all layers

    days, dates_out = get_dates(sim);
    
    #swat = get_states(sim);

    #REW = (swat.SWAT_mm .- minimum(swat.SWAT_mm)) ./ (maximum(swat.SWAT_mm) - minimum(swat.SWAT_mm));

    #return REW

    # get # root zone layers from maximum rooting depth
    max_root_depth = -1.65;

    upper_lim = sim.parametrizedSPAC.soil_discretization.df.Upper_m;
    lower_lim = sim.parametrizedSPAC.soil_discretization.df.Lower_m;

    rz_layers = lower_lim[upper_lim .>= max_root_depth];
    nrz = length(rz_layers);

    # soil water content
    swc = get_soil_(:theta, sim, days_to_read_out_d=days); # soil water content
    select!(swc, Not(:time));
    swc = swc[!, 1:nrz]; # only keep soil water content for root zone layers

    # get field capacity and wilting point for available water content
    thf = sim.ODESolution.prob.p.p_soil.p_THETAF[1:nrz]; # field capacity soil water content

    # use swc to derive minimum soil water content for each layer
    # as approximation of wilting point
    # (more correct would be to derive swc from psi_crit using PTF)
    swc_min = minimum.(eachcol(swc));

    awc = thf .- swc_min; # available water capacity per layer
    
    REW_df = deepcopy(swc);
    for i in 1:nrz
        REW_df[!, i] = (swc[!, i] .- swc_min[i]) ./ awc[i];
    end

    REW = mean.(eachrow(REW_df));

    return REW
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

function get_clim(sim)

    met = get_forcing(sim);
    met.date = Date.(met.dates);

    met = select(met, :date, :vappres_kPa, :tmax_degC, :tmin_degC, :prec_mmDay);
    met.month = month.(met.date);
    met.year = year.(met.date);

    met.tmean = (met.tmax_degC .+ met.tmin_degC) ./ 2;

    # calc saturation vapor pressure
    met.Es = 0.61078 .* exp.(17.26939 .* met.tmean ./ (met.tmean .+ 237.3));
    # calc vapor pressure deficit
    met.VPD = met.Es .- met.vappres_kPa;

    return met
end

function combine_fluxes(sim)
    # get all relevant fluxes and combine into single data frame

    df_fluxes = get_fluxes(sim);
    df_fluxes.date = Date.(df_fluxes.dates);
    df_fluxes.trans = df_fluxes.cum_d_tran;
    df_fluxes.pet = df_fluxes.cum_d_ptran;
    df_fluxes.td = df_fluxes.pet .- df_fluxes.trans
    df_fluxes.PLFL = df_fluxes.cum_d_plfl;
    df_fluxes.PLRF = df_fluxes.cum_d_plrf;
    df_fluxes.plpsi_pd = df_fluxes.cum_pd_plpsi;
    df_fluxes.plpsi_md = df_fluxes.cum_md_plpsi;
    df_fluxes.PLFL_Tr_perc = df_fluxes.PLFL ./ df_fluxes.trans * 100;
    df_fluxes.RWUd, = get_RWU_centroid(sim);
    df_fluxes.PLRFd, = get_PLRF_centroid(sim);
    df_fluxes.REW = get_REW(sim);

    eff_swp = get_eff_swp(sim);
    df_fluxes.SWP = eff_swp.swp_eff;

    df_fluxes.RWUd = replace(df_fluxes.RWUd, NaN=>missing);
    df_fluxes.PLRFd = replace(df_fluxes.PLRFd, NaN=>missing);
    df_fluxes.month = month.(df_fluxes.date);
    df_fluxes.year = year.(df_fluxes.date);

    select!(df_fluxes, :date, :trans, :pet, :td, :RWU,
    :PLFL, :PLRF, :plpsi_pd, :plpsi_md, :PLFL_Tr_perc,
    :RWUd, :PLRFd, :REW, :SWP, :month, :year);

    return df_fluxes
end

function obs_fun_twd(twd_comp)

    # split into pd and md
    twd_comp_pd = twd_comp[:, [:date, :TWD_pd, :cum_pd_plpsi]];
    #lwp_comp_pd = twd_comp[:, [:date, :LWP_pd, :cum_pd_plpsi]];
    twd_comp_md = twd_comp[:, [:date, :TWD_md, :cum_md_plpsi]];
    #lwp_comp_md = twd_comp[:, [:date, :LWP_md, :cum_md_plpsi]];

    # remove missing values
    twd_comp_pd = dropmissing(twd_comp_pd);
    #lwp_comp_pd = dropmissing(lwp_comp_pd);
    twd_comp_md = dropmissing(twd_comp_md);
    #lwp_comp_md = dropmissing(lwp_comp_md);

    cc_twd_pd = cor(twd_comp_pd.cum_pd_plpsi, twd_comp_pd.TWD_pd);
    #cc_lwp_pd = cor(lwp_comp_pd.cum_pd_plpsi, lwp_comp_pd.LWP_pd);
    cc_twd_md = cor(twd_comp_md.cum_md_plpsi, twd_comp_md.TWD_md);
    #cc_lwp_md = cor(lwp_comp_md.cum_md_plpsi, lwp_comp_md.LWP_md);

    #return cc_twd_pd, cc_lwp_pd, cc_twd_md, cc_lwp_md
    return cc_twd_pd, cc_twd_md
end

function twd_combine(sim, obs_twd)
    # sim is the LWFBrook90 simulation
    # obs is the observed dendrometer data

    z_plpsi = get_plpsi(sim);

    twd_comp = sort(rightjoin(obs_twd, z_plpsi, on = :date), :date);

    return twd_comp
end

function get_plpsi(sim)
    # retrieve plant water potential from sim

    z = get_fluxes(sim);
    z.date = Date.(z.dates);
    select!(z, :date, :cum_pd_plpsi, :cum_md_plpsi);

    return z
end

# behavioral data
# dendrometer data
obs_twd = CSV.read("../../Data/Pfyn/Pfyn_twd_2011_17.csv", DataFrame);
obs_twd = obs_twd[obs_twd.month .> 4 .&& obs_twd.month .< 12, :]; # filter out winter months

obs_twd_ctr = obs_twd[obs_twd.scenario .== "control", :]; # filter for control scenario
obs_twd_irst = obs_twd[obs_twd.scenario .== "irrigation stop", :]; # filter for irrigation stop scenario

#select!(obs_twd, Not([:scenario, :year, :month, :TWD_pdn, :MDS_norm])); # drop unnecessary columns

# append root parameters to capacitance parameters
par_best_ctr = DataFrame(par_best_ctr);
par_best_irst = DataFrame(par_best_irst);
par_best_ctr.BETAROOT = [0.970584];
par_best_ctr.MAXROOTDEPTH = [-1.86983];
par_best_irst.BETAROOT = [0.965484];
par_best_irst.MAXROOTDEPTH = [-1.63381];

# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_best_ctr[1,:], Date(2010, 1, 1), Date(2024, 12, 31), "LWFB_testcap/control/", "pfynwald", "LWFB_testrun/ctr_cap/", iso=false);
sim_irst = run_LWFB90_param(par_best_irst[1,:], Date(2010, 1, 1), Date(2024, 12, 31), "LWFB_testcap/irr_stop/", "pfynwald", "LWFB_testrun/irst_cap/", irrig=true, iso=false);

# dendro
twd_comp_ctr = twd_combine(sim_ctr, obs_twd_ctr);
twd_comp_irst = twd_combine(sim_irst, obs_twd_irst);

obs_fun_twd(twd_comp_ctr)

# predawn correlations
draw(data(dropmissing(twd_comp_ctr))*mapping(:TWD_pd, :cum_pd_plpsi, color = :date => x -> dayofyear(x))*visual(Scatter, markersize=6),
    scales(X = (; label="Pre-dawn TWD"), Y= (; label="Pre-dawn plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Pre-dawn TWD vs Plant Water Potential", titlealign = :center)
)

draw(data(dropmissing(twd_comp_irst))*mapping(:TWD_pd, :cum_pd_plpsi, color = :date => x -> dayofyear(x))*visual(Scatter, markersize=6),
    scales(X = (; label="Pre-dawn TWD"), Y= (; label="Pre-dawn plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Pre-dawn TWD vs Plant Water Potential", titlealign = :center)
)


# midday correlations
draw(data(dropmissing(twd_comp_ctr))*mapping(:TWD_md, :cum_md_plpsi, color = :date => x -> dayofyear(x))*visual(Scatter, markersize=6),
    scales(X = (; label="Midday TWD"), Y= (; label="Midday plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Midday TWD vs Plant Water Potential", titlealign = :center)
)

draw(data(dropmissing(twd_comp_irst))*mapping(:TWD_md, :cum_md_plpsi, color = :date => x -> dayofyear(x))*visual(Scatter, markersize=6),
    scales(X = (; label="Midday TWD"), Y= (; label="Midday plant water potential (kPa)")),
    figure = (; size=(800, 600), title="Midday TWD vs Plant Water Potential", titlealign = :center)
)

# time series
draw(data(twd_comp_ctr)* 
    (mapping(:date, :TWD_pdn => (x -> -1 * x))*visual(Lines, color="black", label="TWD")+
    mapping(:date, :cum_pd_plpsi => (x -> x/1000))*visual(Lines, color="red", label="Model")),
    scales(X = (; label=""), Y= (; label="Pre-dawn plant water potential (MPa) / -TWDnorm")),
    legend = (; position = :bottom, framevisible = false),
    figure = (; size=(1200, 600), title="Pre-dawn TWDnorm / Plant Water Potential", titlealign = :center)
)


draw(data(twd_comp_irst)*
    (mapping(:date, :TWD_pdn => (x -> -1 * x))*visual(Lines, color="black", label="TWD")+
    mapping(:date, :cum_pd_plpsi => (x -> x/1000))*visual(Lines, color="red", label="Model")),
    scales(X = (; label=""), Y= (; label="Pre-dawn plant water potential (MPa) / -TWDnorm")),
    legend = (; position = :bottom, framevisible = false),
    figure = (; size=(1200, 600), title="Pre-dawn TWDnorm / Plant Water Potential", titlealign = :center)
)


# compare fluxes

df_fluxes_ctr = combine_fluxes(sim_ctr);
df_fluxes_irst = combine_fluxes(sim_irst);

# compare plant water potential against effective soil water potential

draw(
    data(df_fluxes_ctr)*
    (mapping(:date, :SWP => (x -> x/1000))*visual(Lines, label="SWP", linewidth=1.5)+
    mapping(:date, :plpsi_pd => (x -> x/1000))*visual(Lines, color="green", label="Plant ψ", linewidth=1.5)),
    scales(X = (; label=""), Y= (; label="Water Potentials (MPa)")),
    figure = (; size=(1200, 600), title="Comparison between Modelled Plant Water Potential and Effective Soil Water Potential (SWP)")
)

draw(
    data(df_fluxes_irst)*
    (mapping(:date, :SWP => (x -> x/1000))*visual(Lines, label="SWP", linewidth=1.5)+
    mapping(:date, :plpsi_pd => (x -> x/1000))*visual(Lines, color="green", label="Plant ψ", linewidth=1.5)),
    scales(X = (; label=""), Y= (; label="Water Potentials (MPa)")),
    figure = (; size=(1200, 600), title="Comparison between Modelled Plant Water Potential and Effective Soil Water Potential (SWP)")
)


# restrict to growing season
df_flux_grow = df_fluxes_ctr[df_fluxes_ctr.month .> 3 .&& df_fluxes_ctr.month .< 11, :];

# compare RWU depth and PLRF depth
draw(data(df_flux_grow)*
    mapping(:RWUd, :PLRFd, color=:SWP)*visual(Scatter, markersize=8),
    scales(X = (; label="Weighted RWU Depth (mm)"), Y= (; label="Weighted Plant Recharge Depth (mm)"), Color = (; label="Weighted Soil Water Potential (kPa)"))
)
