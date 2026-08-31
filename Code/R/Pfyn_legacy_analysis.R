## model analysis for Pfynwald legacy results

library(tidyverse)
library(lubridate)

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

# behavioral data

# soil water content and soil matric potential
obs_swc <- read_csv("../../Data/Pfyn/Pfyn_swat.csv")
obs_swc = obs_swc %>% filter(date >= as.Date("2016-01-01"), date < as.Date("2025-01-01")) %>%
  mutate(scen = recode_treatment(meta))

obs_swp <- read_csv("../../Data/Pfyn/Pfyn_swp.csv")
obs_swp = obs_swp %>% filter(date >= as.Date("2016-01-01"), date < as.Date("2025-01-01")) %>%
  mutate(scen = recode_treatment(meta))

# Up-scaled sap flow data
obs_sap <- read_csv("../../Data/Pfyn/Pfyn_trans_2011_17.csv")


# behavioral parameter sets

par_ctr <- read_csv("../../Code/Julia/LWFBcal_output/Pfyn_ctr_param_best.csv")
par_irst <- read_csv("../../Code/Julia/LWFBcal_output/Pfyn_irst_param_best.csv")
par_irr <- read_csv("../../Code/Julia/LWFBcal_output/Pfyn_irr_param_best.csv")

# compare parameters between scenarios

# Explicitly exclude the scenario column
parameter_names <- setdiff(names(par_ctr), "scen")

scenario_levels <- c(
  "Control",
  "Irrigation stop",
  "Irrigation"
)


# Create one density plot for each parameter
par_plots <- lapply(parameter_names, function(parameter_name) {
  
  density_data <- bind_rows(
    tibble(
      value = par_ctr[[parameter_name]],
      treatment = "Control"
    ),
    tibble(
      value = par_irst[[parameter_name]],
      treatment = "Irrigation stop"
    ),
    tibble(
      value = par_irr[[parameter_name]],
      treatment = "Irrigation"
    )
  ) %>%
    mutate(
      treatment = factor(
        treatment,
        levels = scenario_levels
      )
    )
  
  ggplot(
    density_data,
    aes(
      x = value,
      colour = treatment,
      fill = treatment
    )
  ) +
    geom_density(
      alpha = 0.5,
      na.rm = TRUE
    ) +
    scale_colour_manual(values = scenario_colours) +
    scale_fill_manual(values = scenario_colours) +
    labs(
      title = parameter_name,
      x = NULL,
      y = "Density"
    ) +
    theme_bw(base_size = 6) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 8),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6)
    )
})


# Arrange the plots in four rows and seven columns
par_grid <- wrap_plots(
  par_plots,
  nrow = 4,
  ncol = 7
)

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

# soil water potential
swp_model_10_80 = soil_model %>% select(date, scen, param, psi_100mm, psi_800mm) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  pivot_longer(cols = c(psi_100mm, psi_800mm), names_to = "depth", values_to = "SWP") %>%
  mutate(depth = if_else(depth == "psi_100mm", 10, 80))

swp_model_10_80_best = swp_model_10_80 |> 
  filter((scen == "Control" & param == scen_best_ctr) | (scen == "Irrigation stop" & param == scen_best_irst))

swp_model_10_80_avg = swp_model_10_80 |> group_by(date, scen, depth) |> 
  summarize(SWP_mean = mean(SWP), SWP05 = quantile(SWP, 0.05), SWP95 = quantile(SWP, 0.95))

ggplot(swp_model_10_80_best, aes(date, SWP, color=as.factor(depth), alpha="Model"))+geom_line()+
  geom_point(data=filter(obs_swp, scen != "Irrigation"), aes(date, SWP, alpha="Obs"), size=0.8, inherit.aes=F)+
  scale_alpha_manual(values=c("Model"=1, "Obs"=0.5))+
  facet_grid(depth~scen, scale="free_y")+theme_bw()+
  labs(x="", y="Soil Water Potential (kPa)", color="Depth", alpha="Source")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12))


# soil water content
swc_model_10_80 = soil_model %>% select(date, scen, param, theta_100mm, theta_800mm) %>%
  filter(date >= as.Date("2014-01-01")) %>%
  pivot_longer(cols = c(theta_100mm, theta_800mm), names_to = "depth", values_to = "VWC") %>%
  mutate(depth = if_else(depth == "theta_100mm", 10, 80))

swc_model_10_80_best = swc_model_10_80 |> 
  filter((scen == "Control" & param == scen_best_ctr) | (scen == "Irrigation stop" & param == scen_best_irst))

ggplot(swc_model_10_80_best, aes(date, VWC, color=as.factor(depth), alpha="Model"))+geom_line()+
  geom_point(data=filter(obs_swc, scen != "Irrigation"), aes(date, VWC, alpha="Obs"), size=0.8, inherit.aes=F)+
  scale_alpha_manual(values=c("Model"=1, "Obs"=0.5))+
  facet_grid(depth~scen, scale="free_y")+theme_bw()+
  labs(x="", y="Soil Water Content (%)", color="Depth", alpha="Source")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12))


# sap flow
ggplot(filter(flux_best, year>2010), aes(date, cum_d_tran, color=scen))+geom_point(size=0.5)+
  #geom_line(aes(date, trans_sm, color=scen, alpha="Model"), linewidth=0.5)+
  geom_point(data=filter(obs_sap, scen!="Irrigation"), aes(date, Tr*2, color="Obs"), size=0.5, alpha=0.5, inherit.aes=F)+
  scale_color_manual(values=scenario_colors)+
  facet_wrap(~scen, ncol=1)+theme_bw()+
  labs(x="", y="Transpiration (mm/day)", color="Scenario")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12))


# RWU depth vs trans
med_rwu = flux_avg |> filter(year>2013, year<2025, month>3, month<11) |> 
  group_by(scen) |> summarize(rwud_mean = mean(rwu_depth, na.rm=T))

ggplot(filter(flux_avg, year>2013, month>3, month<11), aes(cum_d_tran, rwu_depth, color=as.factor(month)))+
  facet_wrap(~scen)
