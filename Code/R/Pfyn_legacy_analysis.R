## model analysis for Pfynwald legacy results

library(tidyverse)
library(lubridate)
library(patchwork)

# helper functions
recode_treatment <- function(x) {
  case_when(
    x == "control" ~ "Control",
    x == "stop"    ~ "Irrigation stop",
    x == "irrigation" ~"Irrigation",
    .default = NA)
}

# color palette
scenario_colors <- c("Control" = "#E69F00", "Irrigation" = "#56B4E9", "Irrigation stop" = "#009E73", "Obs" = "#000000")

# meteo data
met_file = "../../Code/Julia/LWFBinput/Pfyn_irrigiso_stop/pfynwald_meteoveg.csv"
met_head = read_csv(met_file, n_max=1, show_col_types=F)
met_names = colnames(met_head)
met = read_csv(met_file, skip=2, col_names=met_names, show_col_types=F)
met$tmean = (met$tmax_degC + met$tmin_degC) / 2
met$vpd = (0.61078 * exp(17.26939 * met$tmean / (met$tmean + 237.3))) - met$vappres_kPa;
met_ann = met %>% mutate(year=year(dates)) %>% filter(year>2013, year<2025) %>% group_by(year) %>% summarize(prec=sum(prec_mmDay))
ggplot(met_ann)+geom_col(aes(year, prec), fill="blue")+theme_bw()+
  scale_y_continuous(position="right")+labs(x="", y="Ann. Prec. (mm)")

# irrigation data
irr = read_csv("../../Data/Pfyn/irrigation.csv")

## behavioral data

# soil water content and soil matric potential
obs_swc <- read_csv("../../Data/Pfyn/Pfyn_swat.csv")
obs_swc = obs_swc %>% filter(date >= as.Date("2016-01-01"), date < as.Date("2025-01-01")) %>%
  mutate(scen = recode_treatment(meta))

obs_swp <- read_csv("../../Data/Pfyn/Pfyn_swp.csv")
obs_swp = obs_swp %>% filter(date >= as.Date("2016-01-01"), date < as.Date("2025-01-01")) %>%
  mutate(scen = recode_treatment(meta))

# upscaled sap flow data
obs_sap <- read_csv("../../Data/Pfyn/Pfyn_trans_2011_17.csv")

## compare parameters between scenarios

# behavioral parameter sets
par_ctr <- read_csv("../../Code/Julia/LWFBcal_output/Pfyn_ctr_param_best.csv")
par_irst <- read_csv("../../Code/Julia/LWFBcal_output/Pfyn_irst_param_best.csv")
par_irr <- read_csv("../../Code/Julia/LWFBcal_output/Pfyn_irr_param_best.csv")

# explicitly choose and name parameters for comparison
param_names = c("GLMAX", "CVPD", "R5", "PSICR","MXKPL","BETAROOT")
param_title = c("Maximum stomatal conductance (cm/s)", "VPD at 50% stomatal closure (kPa)", "Radiation at 50% stomatal closure (W/m^2)", 
                "Critical leaf water potential (MPa)", "Maximum plant conductivity (mm/d/MPa)", "Root distribution beta parameter (-)")

# create one density plot for each parameter
par_plots <- lapply(param_names, function(parameter_name, param_title) {
  
  density_data <- bind_rows(tibble(value = par_ctr[[parameter_name]], treatment = "Control"),
                            tibble(value = par_irst[[parameter_name]], treatment = "Irrigation stop"),
                            tibble(value = par_irr[[parameter_name]], treatment = "Irrigation")) %>%
    mutate(treatment = factor(treatment, levels=c("Control","Irrigation stop","Irrigation")))
  
  parameter_title = param_title[which(param_names == parameter_name)]
  
  if (parameter_name == "PSICR") {
    ggplot(density_data, aes(x=value, colour=treatment, fill=treatment))+geom_density(alpha=0.5) +
      scale_colour_manual(values=scenario_colors)+scale_fill_manual(values=scenario_colors) +
      labs(title=parameter_title, x=parameter_name, y="", color="Scenario", fill="Scenario")+theme_bw()+
      theme(legend.position="inside", legend.position.inside=c(0.75, 0.75), 
            plot.title=element_text(size=12, hjust=0.5), 
            legend.title=element_text(size=10), legend.text=element_text(size=10),
            axis.title=element_text(size = 10), axis.text=element_text(size = 10))
    
  } else {
    ggplot(density_data, aes(x=value, colour=treatment, fill=treatment))+geom_density(alpha=0.5) +
      scale_colour_manual(values=scenario_colors)+scale_fill_manual(values=scenario_colors) +
      labs(title=parameter_title, x=parameter_name, y="")+theme_bw()+
      theme(legend.position="none", plot.title=element_text(size=12, hjust=0.5), 
            axis.title=element_text(size = 10), axis.text=element_text(size = 10))
  }
  
  }, param_title)

# arrange the plots
par_grid <- wrap_plots(par_plots, nrow = 2, ncol = 3)
par_grid

# model output
flux_ctr = read_csv("../../Code/Julia/LWFBoutput/Pfyn_ctr_legacy_flux_output.csv")
soil_ctr = read_csv("../../Code/Julia/LWFBoutput/Pfyn_ctr_legacy_soil_output.csv")
flux_irst = read_csv("../../Code/Julia/LWFBoutput/Pfyn_irst_legacy_flux_output.csv")
soil_irst = read_csv("../../Code/Julia/LWFBoutput/Pfyn_irst_legacy_soil_output.csv")

flux_ctr$scen = "Control"
soil_ctr$scen = "Control"
flux_irst$scen = "Irrigation stop"
soil_irst$scen = "Irrigation stop"

flux_model = rbind(flux_ctr, flux_irst)
soil_model = rbind(soil_ctr, soil_irst)

flux_model$year = year(flux_model$date)
flux_model$month = month(flux_model$date)

# best behavioral scenarios
scen_best_ctr = 23738
scen_best_irst = 1505
scen_best_irr = 18295

flux_best = flux_model |> 
  filter((scen == "Control" & param == scen_best_ctr) | (scen == "Irrigation stop" & param == scen_best_irst))

flux_avg = flux_model |> group_by(date, scen, year, month) |> 
  summarize_all(list(mean))

# compare model output against observations

NSE = function(sim, obs) {
  NSE = 1 - (sum((obs - sim)^2) / sum((obs - mean(obs))^2))
}

obj_fun_swp = function(sim, obs) {
  # calculate NSE for soil water potential
  
  obs_10cm = na.omit(filter(obs, depth==10, date %in% sim$date))
  obs_80cm = na.omit(filter(obs, depth==80, date %in% sim$date))
  
  sim_10cm = filter(sim, depth==10, date %in% obs_10cm$date)
  sim_80cm = filter(sim, depth==80, date %in% obs_80cm$date)
  
  nse10_cal = NSE(sim_10cm$SWP[sim_10cm$date < "2022-01-01"], obs_10cm$SWP[obs_10cm$date < "2022-01-01"])
  nse80_cal = NSE(sim_80cm$SWP[sim_80cm$date < "2022-01-01"], obs_80cm$SWP[obs_80cm$date < "2022-01-01"])
  nse10_val = NSE(sim_10cm$SWP[sim_10cm$date >= "2022-01-01"], obs_10cm$SWP[obs_10cm$date >= "2022-01-01"])
  nse80_val = NSE(sim_80cm$SWP[sim_80cm$date >= "2022-01-01"], obs_80cm$SWP[obs_80cm$date >= "2022-01-01"])
  
  return (list(nse10_cal=nse10_cal, nse80_cal=nse80_cal, nse10_val=nse10_val, nse80_val=nse80_val))
}

obj_fun_swc = function(sim, obs) {
  # calculate NSE for soil water content  
  
  obs_10cm = na.omit(filter(obs, depth==10, date %in% sim$date))
  obs_80cm = na.omit(filter(obs, depth==80, date %in% sim$date))
  
  sim_10cm = filter(sim, depth==10, date %in% obs_10cm$date)
  sim_80cm = filter(sim, depth==80, date %in% obs_80cm$date)
  
  nse10_cal = NSE(sim_10cm$VWC[sim_10cm$date < "2022-01-01"], obs_10cm$VWC[obs_10cm$date < "2022-01-01"])
  nse80_cal = NSE(sim_80cm$VWC[sim_80cm$date < "2022-01-01"], obs_80cm$VWC[obs_80cm$date < "2022-01-01"])
  nse10_val = NSE(sim_10cm$VWC[sim_10cm$date >= "2022-01-01"], obs_10cm$VWC[obs_10cm$date >= "2022-01-01"])
  nse80_val = NSE(sim_80cm$VWC[sim_80cm$date >= "2022-01-01"], obs_80cm$VWC[obs_80cm$date >= "2022-01-01"])
  
  return (list(nse10_cal=nse10_cal, nse80_cal=nse80_cal, nse10_val=nse10_val, nse80_val=nse80_val))
}

# soil water potential
swp_model_10_80 = soil_model %>% select(date, scen, param, psi_100mm, psi_800mm) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  pivot_longer(cols = c(psi_100mm, psi_800mm), names_to = "depth", values_to = "SWP") %>%
  mutate(depth = if_else(depth == "psi_100mm", 10, 80))

swp_model_10_80_best = swp_model_10_80 |> 
  filter((scen == "Control" & param == scen_best_ctr) | (scen == "Irrigation stop" & param == scen_best_irst))

ctr_swp_met = obj_fun_swp(filter(swp_model_10_80_best, scen=="Control"), filter(obs_swp, scen=="Control"))
irst_swp_met = obj_fun_swp(filter(swp_model_10_80_best, scen=="Irrigation stop"), filter(obs_swp, scen=="Irrigation stop"))

swp_model_10_80_avg = swp_model_10_80 |> group_by(date, scen, depth) |> 
  summarize(SWP_mean = mean(SWP), SWP05 = quantile(SWP, 0.05), SWP95 = quantile(SWP, 0.95))

swp_model_10_80_best$depth_fact = factor(swp_model_10_80_best$depth, levels=c(10, 80), labels=c("10 cm", "80 cm"))
obs_swp$depth_fact = factor(obs_swp$depth, levels=c(10, 80), labels=c("10 cm", "80 cm"))
ggplot(swp_model_10_80_best, aes(date, SWP / 1000, color=scen, alpha="Model"))+geom_line()+
  geom_point(data=filter(obs_swp, scen != "Irrigation"), aes(date, SWP / 1000, alpha="Obs"), size=0.8, inherit.aes=F)+
  scale_alpha_manual(values=c("Model"=1, "Obs"=0.5))+
  scale_color_manual(values=scenario_colors)+
  facet_grid(depth_fact~scen, scale="free_y")+theme_bw()+
  labs(x="", y="Soil Water Potential (MPa)", color="Scenario", alpha="Source")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"), legend.position="bottom")

# soil water content
swc_model_10_80 = soil_model %>% select(date, scen, param, theta_100mm, theta_800mm) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  pivot_longer(cols = c(theta_100mm, theta_800mm), names_to = "depth", values_to = "VWC") %>%
  mutate(depth = if_else(depth == "theta_100mm", 10, 80))

swc_model_10_80_best = swc_model_10_80 |> 
  filter((scen == "Control" & param == scen_best_ctr) | (scen == "Irrigation stop" & param == scen_best_irst))

ctr_swc_met = obj_fun_swc(filter(swc_model_10_80_best, scen=="Control"), filter(obs_swc, scen=="Control"))
irst_swc_met = obj_fun_swc(filter(swc_model_10_80_best, scen=="Irrigation stop"), filter(obs_swc, scen=="Irrigation stop"))

swc_model_10_80_best$depth_fact = factor(swc_model_10_80_best$depth, levels=c(10, 80), labels=c("10 cm", "80 cm"))
obs_swc$depth_fact = factor(obs_swc$depth, levels=c(10, 80), labels=c("10 cm", "80 cm"))
ggplot(swc_model_10_80_best, aes(date, VWC, color=scen, alpha="Model"))+geom_line()+
  geom_point(data=filter(obs_swc, scen != "Irrigation"), aes(date, VWC, alpha="Obs"), size=0.8, inherit.aes=F)+
  scale_alpha_manual(values=c("Model"=1, "Obs"=0.5))+
  scale_color_manual(values=scenario_colors)+
  facet_grid(depth_fact~scen, scale="free_y")+theme_bw()+
  labs(x="", y="Soil Water Content (%)", color="Scenario", alpha="Source")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"), legend.position="bottom")

# sap flow
ggplot()+geom_rect(data=mutate(filter(irr, year %in% c(2011, 2012, 2013)), scen="Irrigation stop"), aes(x=dates, width=1, ymin=-Inf, ymax=Inf), alpha=0.5, fill="lightblue")+
  geom_point(data=filter(flux_best, year>2010), aes(date, cum_d_tran, color=scen), size=1.0, alpha=0.8, inherit.aes=F)+
  geom_line(data=filter(flux_best, year>2010), aes(date, trans_sm, color=scen), linewidth=0.8, inherit.aes=F)+
  geom_point(data=filter(obs_sap, scen!="Irrigation"), aes(date, Tr*2, color="Obs"), size=1.0, shape=16, alpha=0.5, inherit.aes=F)+
  scale_color_manual(values=scenario_colors)+
  facet_wrap(~scen, ncol=1)+theme_bw()+
  labs(x="", y="Transpiration (mm/day)", color="Scenario")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"), legend.position="bottom")


## analytical figures

# annual transpiration
flux_model %>% group_by(year, scen, param) %>% summarize(ann_trans=sum(cum_d_tran)) %>% 
  group_by(year, scen) %>% summarize_at(vars(ann_trans), list(ann_trans_mean=mean, ann_trans_sd=sd)) %>% 
  ggplot()+geom_col(aes(year, ann_trans_mean, fill=scen), position="dodge")+
  geom_errorbar(aes(x=year, ymin=ann_trans_mean-ann_trans_sd, ymax=ann_trans_mean+ann_trans_sd, group=scen), position="dodge", inherit.aes=F)

# RWU depth
ggplot()+geom_rect(data=filter(irr, year>2010, year<2014), aes(x=dates, width=1, ymin=-Inf, ymax=Inf), alpha=0.5, fill="lightblue")+
  geom_line(data=filter(flux_avg, year>2010), aes(date, rwu_depth_sm/10, color=scen), linewidth=0.8, inherit.aes=F)+
  geom_point(data=filter(flux_avg, year>2010), aes(date, rwu_depth/10, color=scen), size=1.0, inherit.aes=F)+
  theme_bw()+labs(x="", y="Weighted-Average Daily RWU Depth (cm)", color="Scenario")+
  scale_color_manual(values=scenario_colors)+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14))

# RWU depth vs trans

# combine with met data
flux_met = left_join(flux_avg, rename(met, date=dates))

# bikini filter
met_bikini = met %>% filter(tmean > 12, vpd < 1, prec_mmDay < 10)
flux_met = filter(flux_met, date %in% met_bikini$dates)

# mean rwu depth
med_rwu = flux_met |> filter(year>2013, year<2025, month>4, month<11) |> 
  group_by(scen) |> summarize(rwud_mean = mean(rwu_depth, na.rm=T))

# median root depth
beta_ctr = mean(par_ctr$BETAROOT)
beta_irst = mean(par_irst$BETAROOT)
med_rwu$root_mean = c(log(0.5)/log(beta_ctr), log(0.5)/log(beta_irst))

ggplot(filter(flux_met, year>2013, month>4, month<11), aes(cum_d_tran, rwu_depth/10, color=as.factor(month)))+geom_point(size=1.8)+
  geom_hline(data=med_rwu, aes(yintercept=rwud_mean/10, alpha="Mean RWU Depth"), color="black", linewidth=1.0, linetype="dashed")+
  geom_hline(data=med_rwu, aes(yintercept=root_mean, alpha="Mean Root Depth"), color="brown", linewidth=1.0, linetype="dashed")+
  facet_wrap(~scen)+scale_y_reverse()+theme_bw()+scale_alpha_manual(values=c(1.0, 1.0))+
  scale_color_manual(values=c("#0072B2", "#009E73", "#D55E00","#CC79A7", "#F0E442", "#56B4E9"))+
  labs(x="Daily transpiration (mm/day)", y="Weighted-Average\nDaily RWU Depth (cm)", color="Month", alpha="")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        legend.position="inside", legend.position.inside=c(0.85, 0.25),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

# color by vpd
ggplot(filter(flux_met, year>2013, month>4, month<11), aes(cum_d_tran, rwu_depth/10, color=vpd))+geom_point(size=1.8)+
  geom_hline(data=med_rwu, aes(yintercept=rwud_mean/10, alpha="Mean RWU Depth"), color="black", linewidth=1.0, linetype="dashed")+
  geom_hline(data=med_rwu, aes(yintercept=root_mean, alpha="Mean Root Depth"), color="brown", linewidth=1.0, linetype="dashed")+
  facet_wrap(~scen)+scale_y_reverse()+theme_bw()+scale_alpha_manual(values=c(1.0, 1.0))+scale_color_viridis_c()+
  labs(x="Daily transpiration (mm/day)", y="Weighted-Average\nDaily RWU Depth (cm)", color="VPD (kPa)", alpha="")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        legend.position="inside", legend.position.inside=c(0.85, 0.25),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))


