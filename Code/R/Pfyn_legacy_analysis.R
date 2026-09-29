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
met = mutate(met, year = year(dates), month=month(dates), day=day(dates))
met$tmean = (met$tmax_degC + met$tmin_degC) / 2
met$vpd = (0.61078 * exp(17.26939 * met$tmean / (met$tmean + 237.3))) - met$vappres_kPa;

# precip plots
# daily
met = met %>% mutate(date_fill=as.Date(paste("2000", month, day, sep="-")))
ggplot(filter(met, year>2013, year<2025), aes(date_fill, prec_mmDay))+geom_col()+
  theme_bw()+scale_x_date(date_labels="%b")+
  facet_wrap(~year)+labs(x="",y="Daily Precipitation (mm/day)")

# annual total
met_ann = met %>% filter(year>2013, year<2025) %>% group_by(year) %>% summarize(prec=sum(prec_mmDay))
ggplot(met_ann)+geom_col(aes(year, prec), fill="blue")+theme_bw()+
  scale_y_continuous(position="right")+labs(x="", y="Ann. Prec. (mm)")

# annual cum sum
met_list = split(met, met$year)
met_cs = do.call(rbind, lapply(met_list, function(df) {
  df$precip_cs = cumsum(df$prec_mmDay)
  return(df)
}))
ggplot(filter(met_cs, year>2013, year<2025), aes(date_fill, precip_cs, group=as.factor(year), color=as.factor(year)))+
  geom_line()+theme_bw()+scale_x_date(date_labels="%b %d")+labs(x="", y="Precipitation (mm)", color="Year")
ggplot(filter(met_cs, year>2010), aes(date_fill, precip_cs))+geom_line()+
  facet_wrap(~year)+theme_bw()+scale_x_date(date_labels="%b")+labs(x="",y="Precipitation (mm)")

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
param_names = c("GLMAX", "CVPD", "PSICR","MXKPL","BETAROOT", "MAXROOTDEPTH")
param_title = c("a) Maximum stomatal conductance (m/s)", "b) VPD at 50% stomatal closure (kPa)", "c) Critical leaf water potential (MPa)",
                "d) Maximum plant conductivity (mm/d/MPa)", "e) Root distribution beta parameter (-)", "f) Maximum rooting depth (m)")

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
  summarize_all(list(mean)) |> ungroup()

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
ggplot()+geom_col(data=mutate(filter(irr, year %in% c(2011, 2012, 2013)), scen="Irrigation stop"), aes(x=dates, y=Inf), alpha=0.5, fill="lightblue")+
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
ggplot()+geom_col(data=filter(irr, year>2010, year<2014), aes(x=dates, y=Inf), alpha=0.5, fill="lightblue")+
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

summer_cols <- c("#0072B2", "#56B4E9", "#CC79A7", "#D55E00", "#E69F00", "#F0E442")
summer_mons = c("May", "June", "July", "August", "September", "October")

ggplot(filter(flux_met, year>2013, month>4, month<11), aes(cum_d_tran, rwu_depth/10, color=as.factor(month)))+geom_point(size=1.8)+
  geom_hline(data=med_rwu, aes(yintercept=rwud_mean/10, alpha="Mean RWU Depth"), color="black", linewidth=1.0, linetype="dashed")+
  geom_hline(data=med_rwu, aes(yintercept=root_mean, alpha="Mean Root Depth"), color="brown", linewidth=1.0, linetype="dashed")+
  facet_wrap(~scen)+scale_y_reverse()+theme_bw()+scale_alpha_manual(values=c(1.0, 1.0))+
  scale_color_manual(values=summer_cols, labels=summer_mons)+
  labs(x="Daily Transpiration (mm/day)", y="Weighted-Average\nDaily RWU Depth (cm)", color="Month", alpha="")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        legend.position="inside", legend.position.inside=c(0.85, 0.25),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

# color by vpd
ggplot(filter(flux_met, year>2013, month>4, month<11), aes(cum_d_tran, rwu_depth/10, color=vpd))+geom_point(size=1.8)+
  geom_hline(data=med_rwu, aes(yintercept=rwud_mean/10, alpha="Mean RWU Depth"), color="black", linewidth=1.0, linetype="dashed")+
  geom_hline(data=med_rwu, aes(yintercept=root_mean, alpha="Mean Root Depth"), color="brown", linewidth=1.0, linetype="dashed")+
  facet_wrap(~scen)+scale_y_reverse()+theme_bw()+scale_alpha_manual(values=c(1.0, 1.0))+
  scale_color_distiller(palette="YlOrRd", direction=1)+
  labs(x="Daily transpiration (mm/day)", y="Weighted-Average\nDaily RWU Depth (cm)", color="VPD (kPa)", alpha="")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        legend.position="inside", legend.position.inside=c(0.85, 0.25),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))


# RWU distribution

# calculate RWU percent from trans fluxes
get_RWU_percent <- function(df) {
  
  # Identify layer-specific root water uptake columns
  trans_cols <- names(df)[str_detect(names(df), "^TRANS_\\d+mm$")]
  
  # Extract depths from column names
  depths <- as.integer(str_match(trans_cols, "^TRANS_(\\d+)mm$")[, 2])
  
  # Extract RWU matrix
  rwu <- as.matrix(df[, trans_cols])

  # Remove negative uptake
  rwu[rwu < 0] <- 0
  
  # Calculate total uptake at each timestep
  rwu_total <- rowSums(rwu, na.rm = FALSE)
  
  # Avoid division by zero
  rwu_total[rwu_total == 0] <- NA
  
  # Convert each soil-layer uptake to fraction of total uptake
  rwu_percent <- rwu / rwu_total
  
  # Convert back to data frame
  rwu_percent <- as.data.frame(rwu_percent)
  
  # Preserve depth-based column names
  names(rwu_percent) <- paste0("TRANS_",depths,"mm")
  
  # Add date
  rwu_percent <- rwu_percent %>% mutate(date = df$date, .before = 1)
  
  return(rwu_percent)
}

# Group depth-specific variables into custom depth bins
mod_depth_bins <- function(df, depth_bins) {
  
  # Identify depth columns
  depth_cols <- setdiff(names(df),"date")
  
  # Extract depth from column names
  depths <- str_match(depth_cols,"_(\\d+)mm$")[, 2]
  
  # Convert to long format
  df_long <- df %>% rename_with(~ depths, all_of(depth_cols)) %>%
    pivot_longer(cols = -date, names_to = "depth", values_to = "value") %>%
    mutate(depth = as.integer(depth),
      depth_bin = cut(depth - 1, breaks = depth_bins, include.lowest=T, right=T, dig.lab=5))
  
  # Aggregate within depth bins
  df_long_d <- df_long %>% group_by(date, depth_bin) %>%
      summarise(RWU_per = sum(value, na.rm = TRUE), .groups = "drop")

  return(df_long_d)
}

# Calculate RWU fractions
ctr_rwu_per = get_RWU_percent(filter(flux_best, year>2013, scen=="Control"))
irst_rwu_per = get_RWU_percent(filter(flux_best, year>2013, scen=="Irrigation stop"))

# Define custom depth bins
depth_bins = c(0,200,400,600,800,1000,1200,1600,2000)
depth_bin_levels = cut(depth_bins[-1] - 1, depth_bins, include.lowest=T, right=T, dig.lab=5)

# Sum RWU fractions within the new depth bins
ctr_rwu_dist <- mod_depth_bins(ctr_rwu_per, depth_bins)
irst_rwu_dist <- mod_depth_bins(irst_rwu_per, depth_bins)

# Monthly mean RWU distribution
ctr_rwu_dist_mon <- ctr_rwu_dist %>% mutate(month = month(date)) %>% group_by(month, depth_bin) %>%
  summarise(RWU_per = mean(RWU_per, na.rm = TRUE), .groups = "drop") %>%
  mutate(scen = "Control")
irst_rwu_dist_mon <- irst_rwu_dist %>% mutate(month = month(date)) %>% group_by(month, depth_bin) %>%
  summarise(RWU_per = mean(RWU_per, na.rm = TRUE), .groups = "drop") %>%
  mutate(scen = "Irrigation stop")

rwu_dist_mon = rbind(ctr_rwu_dist_mon, irst_rwu_dist_mon)

# root fractions
root_frac_ctr = beta_ctr^(depth_bins[-length(depth_bins)] / 10) - beta_ctr^(depth_bins[-1] / 10)
root_frac_irst = beta_irst^(depth_bins[-length(depth_bins)] / 10) - beta_ctr^(depth_bins[-1] / 10)

df_root = rbind(tibble(depth_bin = depth_bin_levels, root_frac = root_frac_ctr * 100, scen = "Control"),
                tibble(depth_bin = depth_bin_levels, root_frac = root_frac_irst * 100, scen = "Irrigation stop"))

ggplot(filter(rwu_dist_mon, month>4, month<11))+
  geom_path(aes(x=RWU_per*100, y=fct_rev(depth_bin), color=as.factor(month), group=as.factor(month)), linewidth=0.8)+
  geom_path(data=df_root, aes(x=root_frac, y=fct_rev(depth_bin), group=scen, alpha="Beta Root Distribution"), color="black", linewidth=1.2, inherit.aes=F) +
  facet_wrap(~scen)+theme_bw()+scale_alpha_manual(values=c(0.8))+
  scale_color_manual(values=summer_cols, labels=summer_mons)+
  labs(x="% Contribution to Transpiration & Cumulative Root Percentage", y="RWU Depth Bin (mm)", colour="Month", alpha="")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        legend.position="inside", legend.position.inside=c(0.85, 0.25),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))+
  guides(color=guide_legend(order = 1), alpha=guide_legend(order = 2))


# water balance
flux_wb = flux_avg |> mutate(int_loss = cum_d_irvp + cum_d_isvp) |> 
  select(month, year, scen, cum_d_prec, int_loss, cum_d_slvp, cum_d_tran, td, flow) |>
  group_by(month, year, scen) |> summarize_all(list(sum)) |> filter(year > 2013) |> 
  group_by(month, scen) |> summarize_all(list(mean)) |> ungroup() |> select(-year)

flux_wb_long = flux_wb |> pivot_longer(-c(month, scen), names_to="Flux", values_to="Val")

flux_wb_long$Flux = factor(flux_wb_long$Flux, levels=c("cum_d_prec", "td", "cum_d_tran", "int_loss", "cum_d_slvp", "flow"),
                          labels=c("Precipitation", "Transpiration deficit", "Actual transpiration", "Interception loss", "Soil evaporation", "Drainage"))

wb_colors = c("Transpiration deficit" = "red2", "Actual transpiration" = "lightgreen",
              "Interception loss" = "forestgreen", "Soil evaporation" = "khaki3", 
              "Drainage" = "darkblue", "Precipitation" = "black")

ggplot(filter(flux_wb_long, Flux!="Precipitation"), aes(month, Val, fill=Flux, color=Flux))+geom_col()+
  geom_line(data=filter(flux_wb_long, Flux=="Precipitation"), aes(month, Val, color=Flux))+
  facet_wrap(~scen)+theme_bw()+labs(x="", y="Water flux (mm/month)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend(reverse=T))+
  scale_color_manual(values=wb_colors, guide=guide_legend(reverse=T))+
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(1,4,7,10)])+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

# annual
flux_wb_yr = flux_avg |> mutate(int_loss = cum_d_irvp + cum_d_isvp) |> 
  select(year, scen, cum_d_prec, int_loss, cum_d_slvp, cum_d_tran, td, flow) |>
  group_by(year, scen) |> summarize_all(list(sum)) |> filter(year>2013)

flux_wb_yr_long = flux_wb_yr |> pivot_longer(-c(year, scen), names_to="Flux", values_to="Val")

flux_wb_yr_long$Flux = factor(flux_wb_yr_long$Flux, levels=c("cum_d_prec", "td", "cum_d_tran", "int_loss", "cum_d_slvp", "flow"),
                           labels=c("Precipitation", "Transpiration deficit", "Actual transpiration", "Interception loss", "Soil evaporation", "Drainage"))

ggplot(filter(flux_wb_yr_long, Flux!="Precipitation"), aes(year, Val, fill=Flux, color=Flux))+geom_col()+
  geom_line(data=filter(flux_wb_yr_long, Flux=="Precipitation"), aes(year, Val, color=Flux))+
  facet_wrap(~scen)+theme_bw()+labs(x="", y="Annual water flux (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend(reverse=T))+
  scale_color_manual(values=wb_colors, guide=guide_legend(reverse=T))+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

# yearly scenario diff
flux_wb_yr_comp = flux_wb_yr_long |> pivot_wider(names_from=scen, values_from=Val) |> 
  mutate(scen_diff = `Irrigation stop` - Control) |> filter(Flux != "Precipitation")

ggplot(flux_wb_yr_comp, aes(year, scen_diff, fill=Flux, color=Flux))+geom_col()+
  theme_bw()+labs(x="", y="Difference in annual water flux\n Irrigation stop - Control (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend())+
  scale_color_manual(values=wb_colors, guide=guide_legend())+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

# monthly scenario diff
flux_wb_mo_comp = flux_wb_long |> pivot_wider(names_from=scen, values_from=Val) |> 
  mutate(scen_diff = `Irrigation stop` - Control) |> filter(Flux != "Precipitation")

ggplot(flux_wb_mo_comp, aes(month, scen_diff, fill=Flux))+geom_col()+
  theme_bw()+labs(x="", y="Difference in monthly water flux\n Irrigation stop - Control (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend())+
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(1,4,7,10)])+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))


# water balance sensitivity to LAI of irrigation stop
flux_wb_irst = read_csv("../Julia/LWFB_LAIsens/flux_wb_irst_laisens.csv")
flux_wb_irst_yr = flux_wb_irst |> mutate(year=year(date)) |> select(-date) |> 
  group_by(year, sens) |> summarize_all(list(sum)) |> filter(year>2013)

flux_wb_irst_sens = flux_wb_irst_yr |> pivot_longer(-c(year, sens), names_to="Flux", values_to="Val") |> 
  pivot_wider(names_from=sens, values_from=Val) |> filter(Flux != "cum_d_prec")

flux_wb_irst_sens$Flux = factor(flux_wb_irst_sens$Flux, levels=c("cum_d_prec", "td", "cum_d_tran", "int_loss", "cum_d_slvp", "flow"),
                              labels=c("Precipitation", "Transpiration deficit", "Actual transpiration", "Interception loss", "Soil evaporation", "Drainage"))

flux_wb_irst_sens = left_join(flux_wb_irst_sens, select(flux_wb_yr_comp, year, Flux, Control))
flux_wb_irst_sens = mutate(flux_wb_irst_sens, pos_diff = positive - Control,
                           neg_diff = negative - Control, def_diff = default - Control)

flux_wb_irst_long = flux_wb_irst_sens |> select(year, Flux, contains("diff")) |> 
  pivot_longer(-c(year, Flux), names_to="Sens", names_pattern="(.*)_diff")

flux_wb_irst_long$Sens = with(flux_wb_irst_long, 
                              case_when(Sens == "def" ~ "Default", 
                                        Sens == "pos" ~ "Pos. legacy effect",
                                        Sens == "neg" ~ "Neg. legacy effect"))

# show 3 trajectories explicitly
ggplot(flux_wb_irst_long, aes(year, value, fill=Flux))+geom_col()+facet_wrap(~Sens)+
  theme_bw()+labs(x="", y="Difference in annual water flux\nIrrigation stop - Control (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend())+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))


# water balance comparison against climate change scenarios
flux_ctr_cc = read_csv("../../Code/Julia/LWFBoutput/Pfyn_ctr_legacy_cc_flux_output.csv")
flux_irst_cc = read_csv("../../Code/Julia/LWFBoutput/Pfyn_irst_legacy_cc_flux_output.csv")
flux_irr_cc = read_csv("../../Code/Julia/LWFBoutput/Pfyn_irr_legacy_cc_flux_output.csv")
flux_irr = read_csv("../../Code/Julia/LWFBoutput/Pfyn_irr_legacy_flux_output.csv")

# control
flux_mo_ctr_cc = flux_ctr_cc |> mutate(year=year(date), month=month(date), int_loss = cum_d_irvp + cum_d_isvp) |> 
  group_by(param, month, year) |> summarize_at(vars(cum_d_tran, td, int_loss, cum_d_slvp, flow), list(sum))
flux_mo_ctr = flux_ctr |> mutate(year=year(date), month=month(date), int_loss = cum_d_irvp + cum_d_isvp) |> 
  group_by(param, month, year) |> summarize_at(vars(cum_d_tran, td, int_loss, cum_d_slvp, flow), list(sum))

# monthly flux differences  
flux_mo_ctr_comp_cc = flux_mo_ctr_cc
flux_mo_ctr_comp_cc[-(1:3)] = flux_mo_ctr_cc[-(1:3)] - flux_mo_ctr[-(1:3)]
flux_mo_ctr_comp_cc = flux_mo_ctr_comp_cc |> filter(year>2013) |> 
  group_by(month) |> summarize_all(list(mean)) |> ungroup() |> 
  select(-c(param, year)) |> mutate(scen="Control")

# irrigation stop
flux_mo_irst_cc = flux_irst_cc |> mutate(year=year(date), month=month(date), int_loss = cum_d_irvp + cum_d_isvp) |> 
  group_by(param, month, year) |> summarize_at(vars(cum_d_tran, td, int_loss, cum_d_slvp, flow), list(sum))
flux_mo_irst = flux_irst |> mutate(year=year(date), month=month(date), int_loss = cum_d_irvp + cum_d_isvp) |> 
  group_by(param, month, year) |> summarize_at(vars(cum_d_tran, td, int_loss, cum_d_slvp, flow), list(sum))

# monthly flux differences
flux_mo_irst_comp_cc = flux_mo_irst_cc
flux_mo_irst_comp_cc[-(1:3)] = flux_mo_irst_cc[-(1:3)] - flux_mo_irst[-(1:3)]
flux_mo_irst_comp_cc = flux_mo_irst_comp_cc |> filter(year>2013) |> 
  group_by(month) |> summarize_all(list(mean)) |> ungroup() |> 
  select(-c(param, year)) |> mutate(scen="Irrigation stop")

# irrigation
flux_mo_irr_cc = flux_irr_cc |> mutate(year=year(date), month=month(date), int_loss = cum_d_irvp + cum_d_isvp) |> 
  group_by(param, month, year) |> summarize_at(vars(cum_d_tran, td, int_loss, cum_d_slvp, flow), list(sum))
flux_mo_irr = flux_irr |> mutate(year=year(date), month=month(date), int_loss = cum_d_irvp + cum_d_isvp) |> 
  group_by(param, month, year) |> summarize_at(vars(cum_d_tran, td, int_loss, cum_d_slvp, flow), list(sum))

# monthly flux differences
flux_mo_irr_comp_cc = flux_mo_irr_cc
flux_mo_irr_comp_cc[-(1:3)] = flux_mo_irr_cc[-(1:3)] - flux_mo_irr[-(1:3)]
flux_mo_irr_comp_cc = flux_mo_irr_comp_cc |> filter(year>2013) |> 
  group_by(month) |> summarize_all(list(mean)) |> ungroup() |> 
  select(-c(param, year)) |> mutate(scen="Irrigation")


# monthly water balance differences for all treatments
flux_mo_comp_cc = rbind(flux_mo_ctr_comp_cc, flux_mo_irst_comp_cc)
flux_mo_comp_cc = rbind(flux_mo_comp_cc, flux_mo_irr_comp_cc)
flux_mo_comp_cc = flux_mo_comp_cc |> pivot_longer(-c(month, scen), names_to="Flux", values_to="Val")

flux_mo_comp_cc$Flux = factor(flux_mo_comp_cc$Flux, levels=c("td", "cum_d_tran", "int_loss", "cum_d_slvp", "flow"),
                              labels=c("Transpiration deficit", "Actual transpiration", "Interception loss", "Soil evaporation", "Drainage"))
flux_mo_comp_cc$scen = factor(flux_mo_comp_cc$scen, levels=c("Control", "Irrigation stop", "Irrigation"))

ggplot(flux_mo_comp_cc, aes(month, Val, fill=Flux, color=Flux))+geom_col()+
  facet_wrap(~scen)+theme_bw()+labs(x="", y="Difference in monthly water flux\n Future - Historical (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend())+
  scale_color_manual(values=wb_colors, guide=guide_legend())+
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(1,4,7,10)])+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))


# monthly differences between irrigation stop and control
flux_mo_ctr_cc2 = flux_mo_ctr_cc |> filter(year>2020) |> 
  group_by(month) |> summarize_all(list(mean)) |> ungroup() |> 
  select(-c(param, year))

flux_mo_irst_cc2 = flux_mo_irst_cc |> filter(year>2020) |> 
  group_by(month) |> summarize_all(list(mean)) |> ungroup() |> 
  select(-c(param, year))

flux_mo_comp_cc2 = flux_mo_ctr_cc2
flux_mo_comp_cc2[-1] = flux_mo_irst_cc2[-1] - flux_mo_ctr_cc2[-1]
flux_mo_comp_cc2 = flux_mo_comp_cc2 |> pivot_longer(-month, names_to="Flux", values_to="Val")

flux_mo_comp_cc2$Flux = factor(flux_mo_comp_cc2$Flux, levels=c("td", "cum_d_tran", "int_loss", "cum_d_slvp", "flow"),
                              labels=c("Transpiration deficit", "Actual transpiration", "Interception loss", "Soil evaporation", "Drainage"))

ggplot(flux_mo_comp_cc2, aes(month, Val, fill=Flux, color=Flux))+geom_col()+
  theme_bw()+labs(x="", y="Difference in monthly water flux\n Irrigation stop - Control (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend())+
  scale_color_manual(values=wb_colors, guide=guide_legend())+
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(1,4,7,10)])+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

# yearly differences between irrigation stop and control
flux_yr_ctr_cc2 = flux_mo_ctr_cc |> filter(year>2013) |> 
  group_by(year, param) |> summarize_all(list(sum)) |> ungroup() |> 
  select(-c(param, month)) |> group_by(year) |> summarize_all(list(mean))

flux_yr_irst_cc2 = flux_mo_irst_cc |> filter(year>2013) |> 
  group_by(year, param) |> summarize_all(list(sum)) |> ungroup() |> 
  select(-c(param, month)) |> group_by(year) |> summarize_all(list(mean))

flux_yr_comp_cc2 = flux_yr_ctr_cc2
flux_yr_comp_cc2[-1] = flux_yr_irst_cc2[-1] - flux_yr_ctr_cc2[-1]
flux_yr_comp_cc2 = flux_yr_comp_cc2 |> pivot_longer(-year, names_to="Flux", values_to="Val")

flux_yr_comp_cc2$Flux = factor(flux_yr_comp_cc2$Flux, levels=c("td", "cum_d_tran", "int_loss", "cum_d_slvp", "flow"),
                               labels=c("Transpiration deficit", "Actual transpiration", "Interception loss", "Soil evaporation", "Drainage"))

ggplot(flux_yr_comp_cc2, aes(year, Val, fill=Flux, color=Flux))+geom_col()+
  theme_bw()+labs(x="", y="Difference in annual water flux\n Irrigation stop - Control (mm)")+
  scale_fill_manual(values=wb_colors, guide=guide_legend(reverse=T))+
  scale_color_manual(values=wb_colors, guide=guide_legend(reverse=T))+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12, face="bold"))

