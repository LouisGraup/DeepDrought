# calculate PET using multiple methods
# from Pfynwald met data

library(tidyverse)
library(lubridate)
library(humidity)
library(SPEI)

met = read_csv("../../Data/Pfyn/meteo/meteo_irr_Control.csv")

# thornthwaite method uses monthly temperatures

met_mth_yr = met %>% mutate(month=month(dates), year=year(dates)) %>% 
  group_by(month, year) %>% summarize(tmean=mean(tmean))

met_mth_yr$pet_th = thornthwaite(met_mth_yr$tmean, 46.301624) # monthly PET

met_mth_yr$day_th = met_mth_yr$pet_th / 30 # average daily value

met_mth = met_mth_yr %>% select(-year) %>%
  group_by(month) %>% summarize_all(list(mean)) # monthly summary

met_yr = met_mth_yr %>% group_by(year) %>% summarize(pet_yr=sum(pet_th)) %>% 
  filter(year < 2025)


# penman-monteith method

# take met from Richard's data
sap_vpd = readRDS("../../Data/Pfyn/data for Pfynwald/2024/VPDrought_sfden_2025-06-18.RDS")
sap_vpd$Date.Time = as.POSIXct(sap_vpd$Date.Time, tz="CET", format="%Y-%m-%d %H:%M:%S")
sap_vpd$date = as.Date(sap_vpd$Date.Time, tz="CET")

met_hourly = sap_vpd %>% filter(Total_Sap_Flow_kg_h > 0) %>% 
  mutate(datehour=floor_date(Date.Time, "1 hour")) %>% 
  group_by(datehour, Treatment) %>% 
  summarize_at(vars(sr_wm.2, ws_ms.1, temp_degreeC, vpd_kPa), list(mean)) %>% 
  filter(Treatment == "control")

# also need sap flow for gs calculation
sap_vpd_hourly = sap_vpd %>% filter(Total_Sap_Flow_kg_h > 0) %>% 
  mutate(datehour=floor_date(Date.Time, "1 hour")) %>% 
  group_by(datehour, Tree_id, Treatment) %>% summarize(Total_Sap_Flow=mean(Total_Sap_Flow_kg_h, na.rm=T))
sap_vpd_hourly$date = as.Date(sap_vpd_hourly$datehour, tz="CET")

### go get tree specs from sap_flow.R

sap_vpd_ctr = filter(sap_vpd_hourly, Tree_id %in% tree_specs_ctr$tree)
sap_vpd_ctr = left_join(sap_vpd_ctr, select(rename(tree_specs_ctr, Tree_id=tree), Tree_id, meta, sapwood_area, sapwood_area_bin))

# derive sap flux density per unit sapwood area
sap_vpd_ctr$Sap_Flux_Density_kg_cm2_h = sap_vpd_ctr$Total_Sap_Flow / sap_vpd_ctr$sapwood_area
sap_vpd_ctr$Sap_Flux_Density_g_m2_h = sap_vpd_ctr$Sap_Flux_Density_kg_cm2_h * 1e7

sap_tree = filter(sap_vpd_ctr, Tree_id=="1037")
met_hourly = filter(met_hourly, datehour %in% sap_tree$datehour)

# variables for below
e = sap_tree$Sap_Flux_Density_g_m2_h # transpiration (g/m²/h)
t = met_hourly$temp_degreeC # air temperature (°C)
vpd = met_hourly$vpd_kPa # vapor pressure deficit (kPa)
rad = met_hourly$sr_wm.2 * .8 # net radiation assuming albedo = 0.2 as model default (W/m²)
wind = met_hourly$ws_ms.1 # wind speed (m/s)

elv = 615 # elevation (m)

# constants
R      <- 8.314        # J/mol/K
cp_g   <- 1.01         # J/g/K
cp_kg  <- 1013         # J/kg/K
r_air  <- 287.05       # J/kg/K
r_wat  <- 461.5        # J/kg/K
g      <- 9.81         # m/s²
p0     <- 101325       # Pa (sea-level)

# temperature
tK <- t + 273.15

# ==== PRESSURE ADJUSTED FOR ELEVATION ====
tp <- tK + 0.005 * elv / 2
P_Pa <- p0 * exp(-g * elv / (r_air * tp))

# ==== LATENT HEAT OF VAPORIZATION (J/g) ====
lambda_g <- 2500.8 - 2.36 * t + 0.0016 * t^2 - 0.00006 * t^3  # J/g
lambda_kg <- lambda_g * 1000  # J/kg

# ==== PSYCHROMETRIC CONSTANT (Pa/K) ====
gamma_pa <- (P_Pa / 100) * cp_g / lambda_g / (r_air / r_wat) * 100

# ==== AIR DENSITY ====
rho_air_kg <- P_Pa / (r_air * tK)
rho_air_g <- rho_air_kg * 1000

# ==== VPD ====
VPD_Pa <- vpd * 1000

# ==== TRANSPIRATION RATE ====
e_g_m2_s <- e / 3600

# ==== STOMATAL CONDUCTANCE (m/s) ====
gs_m_s <- (lambda_g * e_g_m2_s * gamma_pa) / (rho_air_g * cp_g * VPD_Pa)

# Adjust wind to 2 m height if needed
wind_height = 10
if (wind_height != 2) {
  wind <- wind * log(67.8 * 2 - 5.42) / log(67.8 * wind_height - 5.42)
}

# Slope of saturation vapor pressure curve (Delta, kPa/°C)
Delta <- 4098 * (0.6108 * exp((17.27 * t) / (t + 237.3))) / ((t + 237.3)^2)

# Aerodynamic resistance
ra <- 208 / wind # s/m
ca = 1 / ra

# PET (assuming m/h)
PET = (Delta * rad + rho_air_kg * cp_kg * VPD_Pa / ra) /
  (lambda_kg * (Delta + gamma_pa * (1 + ca / gs_m_s)))
PET_mm_h = PET * 1000

PET_df = data.frame(date=as.Date(sap_tree$datehour), PET_mm_h=PET_mm_h)
PET_daily = PET_df %>% group_by(date) %>% summarize(PET=sum(PET_mm_h))

ggplot(PET_daily, aes(date, PET))+geom_line()+theme_bw()+
  labs(x="", y="PET (mm/day)")
