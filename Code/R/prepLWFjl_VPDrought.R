
## prepare LWFBrook90.jl inputs for Pfynwald VPDrought scenarios

library(tidyverse)
library(lubridate)

library(LWFBrook90R)
source("C:/Users/grauplou/Documents/LWFBrook90.jl/src/generate_LWFBrook90jl_Input_mod.R", echo=F)


## meteo inputs for control and vpd manipulation experiment

meteo_Con = read_csv("../../Data/Pfyn/meteo/meteo_irr_Control.csv")
meteo_VPD = read_csv("../../Data/Pfyn/meteo/meteo_irr_VPD.csv")


# separate treatments

# control
meteo_cont = meteo_Con %>% rename(prec = precip_ctrl) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# irrigation stop
meteo_irrstp = meteo_Con %>% rename(prec = precip_irrstp) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# irrigation under ambient conditions
meteo_irr_con = meteo_Con %>% rename(prec = precip_irr) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# drought under ambient conditions
meteo_roof_con = meteo_Con %>% rename(prec = precip_roof) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# irrigation under manipulated VPD
meteo_irr_vpd = meteo_VPD %>% rename(prec = precip_irr) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# drought under manipulated VPD
meteo_roof_vpd = meteo_VPD %>% rename(prec = precip_roof) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)


# site data for tree heights

site_df = read_csv("../../Data/Pfyn/siteproperties.csv")
ht_cont = mean(filter(site_df, BEZKM == "control")$height_m)
ht_irr = mean(filter(site_df, BEZKM == "irrigation")$height_m)
ht_stp = mean(filter(site_df, BEZKM == "stop")$height_m)

# modeled LAI data from prepLWFjl.R

LAI_df = read_csv("../../Data/Pfyn/LAI_ext.csv")

# soil data

soil_df = read_csv("../../Data/Pfyn/soil_hydraulic.csv")


## use LWFBrook90R to prepare input files for julia version

# options common to all scenarios
opt = set_optionsLWFB90(startdate=as.Date("2000-01-01"), enddate=as.Date("2025-06-30"), 
                        root_method="soilvar", budburst_method="Menzel", 
                        leaffall_method="vonWilpert", lai_method="Coupmodel")

# parameters for LAI and climate vary by scenario

# control scenario
par_c = set_paramLWFB90(maxlai=LAI_df$LAI_ctrl, winlaifrac=.6, height=ht_cont, height_ini=ht_cont, 
                        coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, bypar=1, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_control","pfynwald",".", options_b90=opt, param_b90=par_c, climate=meteo_cont, soil=soil_df)


# drought scenarios use same parameters as control

# drought scenario under ambient climate
generate_LWFBrook90jl_Input("Pfyn_drought_ambient","pfynwald",".", options_b90=opt, param_b90=par_c, climate=meteo_roof_con, soil=soil_df)

# drought scenario under manipulated VPD
generate_LWFBrook90jl_Input("Pfyn_drought_VPD","pfynwald",".", options_b90=opt, param_b90=par_c, climate=meteo_roof_vpd, soil=soil_df)


# irrigation stop scenario
par_irst = set_paramLWFB90(maxlai=LAI_df$LAI_irrstp, winlaifrac=.6, height=ht_stp, height_ini=ht_stp, 
                           coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, bypar=1, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_irr_stop","pfynwald",".", options_b90=opt, param_b90=par_irst, climate=meteo_irrstp, soil=soil_df)


# irrigation scenarios use same parameters
par_ir = set_paramLWFB90(maxlai=LAI_df$LAI_irr, winlaifrac=.6, height=ht_irr, height_ini=ht_irr, 
                         coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, bypar=1, budburst_species="Pinus sylvestris")

# irrigation under ambient climate
generate_LWFBrook90jl_Input("Pfyn_irrigation_ambient","pfynwald",".", options_b90=opt, param_b90=par_ir, climate=meteo_irr_con, soil=soil_df)

# irrigation under manipulated VPD
generate_LWFBrook90jl_Input("Pfyn_irrigation_VPD","pfynwald",".", options_b90=opt, param_b90=par_ir, climate=meteo_irr_vpd, soil=soil_df)

