## run best performing scenario from LWFBrook90.jl calibration
## for all Pfynwald VPDrought scenarios
## and examine output against observed data

using CSV, DataFrames, Dates, Statistics, RollingFunctions, RCall;
using Plots; gr()
R"""
library(tidyverse)
"""

include("run_LWFB90_param.jl");

# function to retrieve best parameter set for control and irrigation scenarios
function par_best(met, par)
    
    # add combined NSE metrics
    met.swc_nse_com = (met.swc_nse10 + met.swc_nse40 + met.swc_nse60 + met.swc_nse80) / 4;
    met.swp_nse_com = (met.swp_nse10 + met.swp_nse80) / 2;

    met.nse_com = (met.swc_nse_com + met.swp_nse_com) / 2;

    # best overall scenario
    max_idx = argmax(met.nse_com);
    scen_best = met.scen[max_idx];

    par_best = par[scen_best, :];

    return par_best, scen_best
end

function get_dates(sim)
    days = range(sim.ODESolution.prob.tspan...);
    dates_out = LWFBrook90.RelativeDaysFloat2DateTime.(days,sim.parametrizedSPAC.reference_date);
    return days, Date.(dates_out)
end

function get_swc(sim; shape = "long")
    # retrieve soil water content data from sim
    
    days, dates_out = get_dates(sim);

    z = get_soil_(:theta, sim, depths_to_read_out_mm = [100, 800, 1200], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));

    # rename columns
    rename!(z, ["0.1", "0.8", "1.2", "date"]);

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="VWC");
        z_long.depth = parse.(Float64, z_long.depth);

        return z_long
    else
        return z
    end

end

function get_swp(sim; shape="long")
    # retrieve soil water potential data from sim
    
    days, dates_out = get_dates(sim);

    z = get_soil_(:psi, sim, depths_to_read_out_mm = [100, 800, 1200], days_to_read_out_d = days);
    z.date = dates_out;

    select!(z, Not(:time));
    
    # rename columns
    rename!(z, ["0.1", "0.8", "1.2", "date"]);

    if shape=="long"
        # reshape data
        z_long = stack(z, Not(:date), variable_name = "depth", value_name="SWP");
        z_long.depth = parse.(Float64, z_long.depth);
        
        return z_long
    else
        return z
    end

end

function get_sap(sim)
    # retrieve soil water potential data from sim
    
    days, dates_out = get_dates(sim);

    z = get_fluxes(sim);
    z.date = dates_out;
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
    
    return col_RWU_centroid_mm
end

function get_VPD(sim)

    met = get_forcing(sim);
    met.date = Date.(met.dates);

    met = met[met.date .>= Date(2024, 1, 1), :];

    met = select(met, :date, :vappres_kPa, :tmax_degC, :tmin_degC);
    met.month = month.(met.date);

    met.tmean = (met.tmax_degC .+ met.tmin_degC) ./ 2;

    # calc saturation vapor pressure
    met.Es = 0.61078 .* exp.(17.26939 .* met.tmean ./ (met.tmean .+ 237.3));
    # calc vapor pressure deficit
    met.VPD = met.Es .- met.vappres_kPa;

    return met
end

function met_comp(sim, obs)
    # sim is the LWFBrook90 simulation
    # obs is the observed soil moisture data

    z_theta = get_swc(sim);
    z_psi = get_swp(sim);

    z_sim = innerjoin(z_theta, z_psi, on = [:date, :depth]);

    filter!(:date => >(obs.date[1]), z_sim);

    z_sim.src .= "sim"; # add source column

    obs = select(obs, Not(:treatment)); # remove extra columns
    obs.src .= "obs"; # add source column

    # combine observed and simulated data
    met_comp = [obs; z_sim];

    return met_comp

end


# calibration results
met_ctr = CSV.read("LWFBcal_output/metrics_ctr_20250814.csv", DataFrame);
met_irr = CSV.read("LWFBcal_output/metrics_irr_20250814.csv", DataFrame);
par_ctr = CSV.read("LWFBcal_output/param_ctr_20250814.csv", DataFrame);
par_irr = CSV.read("LWFBcal_output/param_irr_20250814.csv", DataFrame);

par_ctr_best, scen_ctr_best = par_best(met_ctr, par_ctr);
par_irr_best, scen_irr_best = par_best(met_irr, par_irr);

met_ctr[scen_ctr_best, :]
met_irr[scen_irr_best, :]


## behavioral data
obs = CSV.read("../../Data/Pfyn/PFY_VPD_swp_swc.csv", DataFrame, missingstring="NA");
rename!(obs, :SWP_corr => :SWP);
obs.depth = -1 * obs.depth;

obs_ctr = obs[obs.treatment .== "control", :];
obs_irr = obs[obs.treatment .== "irrigation", :];
obs_irr_vpd = obs[obs.treatment .== "irrigation_vpd", :];
obs_drt = obs[obs.treatment .== "roof", :];
obs_drt_vpd = obs[obs.treatment .== "roof_vpd", :];

# include extended long-term observations
obs_ext = CSV.read("../../Data/Pfyn/PFY_swpc.csv", DataFrame);
filter!(:date => >=(Date(2024, 1, 1)), obs_ext); # filter out early dates
obs_ext_long = stack(obs_ext, Not(:date, :meta), variable_name = "depth", value_name="SWP");
obs_ext_long.depth = parse.(Int, map(x -> x[(end-3):(end-2)], obs_ext_long.depth)) / 100; # convert depth from var name
obs_ext_long.src .= "long-term obs"; # add source column

# separate scenarios
obs_ext_ctr = obs_ext_long[obs_ext_long.meta .== "control", Not(:meta)];
obs_ext_irr = obs_ext_long[obs_ext_long.meta .== "irrigation", Not(:meta)];

# include multiple roof observations
obs_roof = CSV.read("../../Data/Pfyn/soil_daily_VPDrought.csv", DataFrame,  missingstring="NA");
obs_roof = obs_roof[obs_roof.treatment .== "roof" .|| obs_roof.treatment .== "roof_vpd", :];
select!(obs_roof, :date, :treatment, :depth, :site, :descr, :SWP_corr);
rename!(obs_roof, :SWP_corr => :SWP);
obs_roof.depth = -1 * obs_roof.depth;

# hourly soil water potential
obs_swp = CSV.read("../../Data/Pfyn/soil_hourly_VPDrought.csv", DataFrame,  missingstring="NA");
obs_swp.depth = -1 * obs_swp.depth;
obs_swp.datetime = DateTime.(SubString.(obs_swp.datetime, 1, 19));
obs_swp.date = Date.(obs_swp.datetime);
obs_swp = obs_swp[obs_swp.date .>= Date("2024-04-01") .&& obs_swp.date .< Date("2025-01-01"), :];
obs_swp = obs_swp[obs_swp.treatment .!= "control", :];

obs_swp_pd = obs_swp[hour.(obs_swp.datetime) .== 6, :]; # retrieve pre-dawn swp

# leaf water potential
obs_lwp = CSV.read("../../Data/Pfyn/PFY_lwp.csv", DataFrame,  missingstring="NA");
obs_lwp.treatment = ifelse.(obs_lwp.treat2 .== "irrigated", "irrigation",
    ifelse.(obs_lwp.treat2 .== "irrigated-VPD", "irrigation_vpd",
    ifelse.(obs_lwp.treat2 .== "roof-VPD", "roof_vpd", obs_lwp.treat2)));
dropmissing!(obs_lwp);

obs_lwp_pd = obs_lwp[obs_lwp.wp .== "pd", :]; # pre-dawn lwp
obs_lwp_pd.date = Date.(DateTime.(obs_lwp_pd.MESSTIME, dateformat"y-m-d H:M:S"));
select!(obs_lwp_pd, Not([:twd_n2, :treat2, :wp, :MESSTIME]));

# run LWFBrook90.jl for all scenarios
sim_ctr = run_LWFB90_param(par_ctr_best, Date(2020, 1, 1), Date(2025, 7, 30), "LWFBinput/Pfyn_control/", "pfynwald", "LWFB_VPD/control/");
sim_irr = run_LWFB90_param(par_irr_best, Date(2020, 1, 1), Date(2025, 7, 30), "LWFBinput/Pfyn_irrigation_ambient/", "pfynwald", "LWFB_VPD/irr_amb/");
sim_irr_vpd = run_LWFB90_param(par_irr_best, Date(2020, 1, 1), Date(2025, 7, 30), "LWFBinput/Pfyn_irrigation_VPD/", "pfynwald", "LWFB_VPD/irr_vpd/");
sim_drt = run_LWFB90_param(par_ctr_best, Date(2020, 1, 1), Date(2025, 7, 30), "LWFBinput/Pfyn_drought_ambient/", "pfynwald", "LWFB_VPD/drt_amb/");
sim_drt_vpd = run_LWFB90_param(par_ctr_best, Date(2020, 1, 1), Date(2025, 7, 30), "LWFBinput/Pfyn_drought_VPD/", "pfynwald", "LWFB_VPD/drt_vpd/");

ctr_comp = met_comp(sim_ctr, obs_ctr);
irr_comp = met_comp(sim_irr, obs_irr);
irr_vpd_comp = met_comp(sim_irr_vpd, obs_irr_vpd);
drt_comp = met_comp(sim_drt, obs_drt);
drt_vpd_comp = met_comp(sim_drt_vpd, obs_drt_vpd);

# add extended observations
ctr_swp_comp = [select(ctr_comp, Not(:VWC)); obs_ext_ctr];
irr_swp_comp = [select(irr_comp, Not(:VWC)); obs_ext_irr];

R"""
rdf = $ctr_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Potential Comparison for Control")
"""

R"""
rdf = $ctr_swp_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) + theme_bw() +
    labs(title="Soil Water Potential Comparison for Control Scenario", x="", y="SMP (kPa)", color="Source") +
    theme(plot.title=element_text(hjust=.5), legend.text=element_text(size=12),legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"),axis.text=element_text(size=12), axis.title=element_text(size=14))
"""

R"""
rdf = $ctr_comp
ggplot(rdf, aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Content Comparison for Control")
"""

R"""
rdf = $irr_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation")
"""

R"""
rdf = $irr_swp_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation")
"""

R"""
rdf = $irr_comp
ggplot(rdf, aes(x=date, y=VWC, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Content Comparison for Irrigation")
"""

R"""
rdf = $irr_vpd_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) +
    labs(title="Soil Water Potential Comparison for Irrigation VPD")
"""


R"""
rdf = $drt_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) + theme_bw() + 
    labs(title="Soil Water Potential Comparison for Roof") +
    theme(legend.position=c(.1,.1), plot.title=element_text(hjust=.5))
"""

R"""
rdf = $drt_vpd_comp
ggplot(rdf, aes(x=date, y=SWP, color=src)) + geom_point(size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Roof VPD")
"""

R"""
rdf1 = $drt_comp
rdf2 = $drt_vpd_comp
ggplot(filter(rdf1, src=="sim"), aes(x=date, y=SWP, color="roof")) + geom_point(size=.5) +
    geom_point(data=filter(rdf2, src=="sim"), aes(x=date, y=SWP, color="roof_vpd"), size=.5) +
    facet_wrap(~depth, ncol=1) +
    labs(title="Soil Water Potential Comparison for Drought Scenarios")
"""

R"""
rdf1 = $drt_comp
rdf2 = $drt_vpd_comp
ggplot(rdf1, aes(x=date, y=SWP, color="roof", linetype=as.factor(src))) + geom_line(size=.5) +
    geom_line(data=rdf2, aes(x=date, y=SWP, color="roof_vpd", linetype=as.factor(src)), size=.5) +
    facet_wrap(~depth, scales="free_y", ncol=1) + theme_bw() +
    labs(title="Soil Water Potential Comparison for Drought Scenarios", x="", y="SMP (kPa)", color="Scenario", linetype="Source") + 
    theme(plot.title=element_text(hjust=.5), legend.text=element_text(size=12),legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"),axis.text=element_text(size=12), axis.title=element_text(size=14))
"""

R"""
rdf1 = $drt_comp
rdf2 = $drt_vpd_comp
rdf3 = $obs_roof
ggplot(rdf1, aes(x=date, y=SWP, color="roof", linetype=src)) + geom_line(size=1) +
    geom_line(data=rdf2, aes(x=date, y=SWP, color="roof_vpd", linetype=src), size=1) +
    stat_summary(data=rdf3, geom="ribbon", mapping=aes(x=date, y=SWP, color=treatment, fill=treatment, linetype="obs"), alpha=.3, inherit.aes=F) +
    facet_wrap(~depth, scales="free_y", ncol=1) + theme_bw() +
    labs(title="Soil Water Potential Comparison for Drought Scenarios", x="", y="SMP (kPa)", color="Scenario", fill="Scenario", linetype="Source") + 
    theme(plot.title=element_text(hjust=.5), legend.text=element_text(size=12),legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"),axis.text=element_text(size=12), axis.title=element_text(size=14))
"""

R"""
rdf1 = $drt_comp
rdf2 = $obs_roof
ggplot(filter(rdf1, date>="2025-04-25", depth==.1, src=="sim"), aes(x=date, y=SWP, color="sim")) + geom_line(size=.5) +
    geom_line(data=filter(rdf2, date>="2025-04-25", treatment=="roof", depth==.1), aes(x=date, y=SWP, group=interaction(as.factor(descr), as.factor(site)), linetype=as.factor(site), color=as.factor(descr)), size=.5) +
    labs(title="Soil Water Potential Comparison at 10 cm for Roof Scenario in 2025", x="", y="SMP (kPa)", color="Source", linetype="Plot") + theme_bw() +
    scale_color_manual(values=c("red","purple","green","blue","black")) + guides(linetype="none") +
    theme(plot.title=element_text(hjust=.5), legend.text=element_text(size=12), legend.title=element_text(size=12), axis.text=element_text(size=12), axis.title=element_text(size=14))
"""

drt_comp.treatment .= "roof";
drt_vpd_comp.treatment .= "roof_vpd";

R"""
rdf1 = $drt_comp
rdf1a = $drt_vpd_comp
rdf2 = $obs_roof
ggplot(filter(rdf1, date>="2025-04-25", depth==.1, src=="sim"), aes(x=date, y=SWP, color="sim")) + geom_line(size=.5) +
    geom_line(data=filter(rdf1a, date>="2025-04-25", depth==.1, src=="sim"), aes(x=date, y=SWP, color="sim"), size=.5) +
    geom_line(data=filter(rdf2, date>="2025-04-25", depth==.1), aes(x=date, y=SWP, group=interaction(as.factor(descr), as.factor(site)), color=as.factor(descr)), size=.5) +
    labs(title="Soil Water Potential Comparison at 10 cm for Roof Scenarios in 2025", x="", y="SMP (kPa)", color="Source", linetype="Plot") + theme_bw() +
    scale_color_manual(values=c("red","purple","orange","blue","black")) + facet_wrap(~treatment) + 
    theme(plot.title=element_text(hjust=.5), legend.text=element_text(size=12), legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"), axis.text=element_text(size=12), axis.title=element_text(size=14))
"""

# compare against pre-dawn leaf water potential

swp_irr = get_swp(sim_irr);
swp_irr.treatment .= "irrigation";

swp_irr_vpd = get_swp(sim_irr_vpd);
swp_irr_vpd.treatment .= "irrigation_vpd";

swp_drt = get_swp(sim_drt);
swp_drt.treatment .= "roof";

swp_drt_vpd = get_swp(sim_drt_vpd);
swp_drt_vpd.treatment .= "roof_vpd";

swp_vpd = [swp_irr; swp_irr_vpd; swp_drt; swp_drt_vpd];
swp_vpd = swp_vpd[swp_vpd.date .>= Date(2024, 4, 1) .&& swp_vpd.date .< Date(2025,1,1), :];

R"""
rdf1 = $swp_vpd
rdf2 = $obs_swp_pd
rdf3 = $obs_lwp_pd
ggplot(rdf2, aes(date, SWP_corr/1000, color=as.factor(scaffold), linetype=as.factor(depth), group=interaction(as.factor(scaffold), as.factor(depth))))+geom_line()+
  geom_point(data=rdf3, aes(date, wp_value/10, color=as.factor(scaffold)), inherit.aes=F)+
  geom_line(data=rdf1, aes(x=date, y=SWP/1000, linetype=as.factor(depth), color="sim"), color="black", inherit.aes=F)+
  facet_wrap(~treatment, scales="free_y")+theme_bw()+labs(x="", y="SWP, LWP (MPa)", color="Scaffold", linetype="Depth (m)")+
  ggtitle("Comparison between Observed pre-dawn Leaf Water Potential (LWP) and Soil Water Potential (SWP) with Modelled SWP")+theme(plot.title=element_text(hjust=.5))
"""

WP_comp = leftjoin(obs_lwp_pd, swp_vpd[swp_vpd.depth .== .1, :], on = [:date, :treatment]);

R"""
rdf = $WP_comp
ggplot(rdf, aes(SWP/1000, wp_value/10, color=as.factor(tree)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+theme_bw()+guides(color="none")+facet_wrap(~treatment)
"""

# compare RWU depth across scenarios

days, dates_out = get_dates(sim_ctr);

rwu_comp = DataFrame(date=dates_out);

rwu_comp.control = runmean(get_RWU_centroid(sim_ctr), 7);
rwu_comp.irrigation = runmean(get_RWU_centroid(sim_irr), 7);
rwu_comp.irrigation_vpd = runmean(get_RWU_centroid(sim_irr_vpd), 7);
rwu_comp.roof = runmean(get_RWU_centroid(sim_drt), 7);
rwu_comp.roof_vpd = runmean(get_RWU_centroid(sim_drt_vpd), 7);

rwu_comp = stack(rwu_comp, Not(:date), variable_name="treatment", value_name="RWU");
rwu_comp.RWU = replace(rwu_comp.RWU, NaN=>missing);

R"""
rdf = $rwu_comp
#rdf$treatment = factor(rdf$treatment, levels=c("roof", "roof_vpd", "control", "irrigation", "irrigation_vpd"))
rdf$treatment = factor(rdf$treatment, levels=c("control", "roof", "irrigation"))
ggplot(filter(rdf, date>="2024-01-07"), aes(x=date, y=RWU, color=treatment)) + geom_line() +
    labs(x="", title="RWU Depth Comparison across Scenarios") + theme_bw()
"""

rwu_comp = dropmissing(rwu_comp);
rwu_comp.year = year.(rwu_comp.date);
rwu_yr = combine(groupby(rwu_comp, [:year, :treatment]), :RWU .=> [mean maximum]);

R"""
rdf = $rwu_yr
ggplot(filter(rdf, year>2023), aes(x=year, y=RWU_mean, group=treatment, fill=treatment)) + geom_col(position="dodge") +
    labs(title="Mean RWU Depth Comparison Across Scenarios")
"""


# compare VPD across scenarios
met_drt = get_VPD(sim_drt);
met_drt_vpd = get_VPD(sim_drt_vpd);

R"""
rdf1 = $met_drt
rdf2 = $met_drt_vpd
ggplot(rdf1, aes(x=date, y=tmean, color="roof")) + geom_line() +
    geom_line(data=rdf2, aes(x=date, y=tmean, color="roof_vpd")) +
    labs(title="Temperature Comparison for Drought Scenarios")
"""

R"""
rdf1 = $met_drt
rdf2 = $met_drt_vpd
ggplot(rdf1, aes(x=date, y=vappres_kPa, color="roof")) + geom_line() +
    geom_line(data=rdf2, aes(x=date, y=vappres_kPa, color="roof_vpd")) +
    labs(title="Vapor Pressure Comparison for Drought Scenarios")
"""

R"""
rdf1 = $met_drt
rdf2 = $met_drt_vpd
ggplot(rdf1, aes(x=date, y=VPD, color="roof")) + geom_line() +
    geom_line(data=rdf2, aes(x=date, y=VPD, color="roof_vpd")) +
    labs(title="VPD Comparison for Drought Scenarios")
"""

# actual VPD reductions during summer months

met_drt = met_drt[met_drt.month .>= 5 .&& met_drt.month .<= 8, :];
met_drt_vpd = met_drt_vpd[met_drt_vpd.month .>= 5 .&& met_drt_vpd.month .<= 8, :];

abs_diff = met_drt.VPD - met_drt_vpd.VPD;
rel_diff = abs_diff ./ met_drt_vpd.VPD * 100;

mean(rel_diff) # 15%