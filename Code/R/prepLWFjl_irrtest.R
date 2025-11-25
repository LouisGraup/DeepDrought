
## prepare LWFBrook90.jl inputs for Pfynwald VPDrought scenarios

library(tidyverse)
library(lubridate)

library(LWFBrook90R)
source("C:/Users/grauplou/Documents/LWFBrook90.jl/src/generate_LWFBrook90jl_Input_mod.R", echo=F)

# function to extend time series with leaf out timing for incomplete year
extend_meteoveg = function(meteo_scen, scen) {
  
  # append scenario name to file path
  filename = paste0(scen,"/pfynwald_meteoveg.csv")
  
  # retrieve file header and column names
  met_head = read_csv(filename, n_max=1, show_col_types=F)
  met_names = colnames(met_head)
  
  # read in file with column names and skip header
  met = read_csv(filename, skip=2, col_names=met_names, show_col_types=F)
  
  # extract substitution year for veg data
  veg_fill = met %>% mutate(year=year(dates), month=month(dates)) %>% filter(year==2014, month < 8) %>% select(densef_percent:sai_percent)
  
  # retrieve most recent data from meteo measurements
  meteo_fill = meteo_scen %>% filter(dates >= "2025-01-01") %>% select(-tmean)
  meteo_fill[,2:6] = round(meteo_fill[,2:6], 5)
  
  # combine meteo with veg and add column names
  meteo_fill = cbind(meteo_fill, veg_fill)
  colnames(meteo_fill) = met_names
  
  # append recent data to LWFBrook90 output
  meteo_rep = rbind(met, meteo_fill)
  meteo_rep$dates = as.character(meteo_rep$dates)
  
  # create output data frame with file header and replace in folder
  meteo_out = rbind(met_head, meteo_rep)
  write_csv(meteo_out, filename)
  
}


## meteo inputs for control and vpd manipulation experiment

meteo_Con = read_csv("../../Data/Pfyn/meteo/meteo_irr_Control.csv")
meteo_VPD = read_csv("../../Data/Pfyn/meteo/meteo_irr_VPD.csv")

# use existing control treatment and add irrigation column

met_head = read_csv("../Julia/LWFBinput/Pfyn_irrigation_ambient/pfynwald_meteoveg.csv", n_max=1, show_col_types=F)
met_names = colnames(met_head)

# read in file with column names and skip header
#metveg = read_csv("../Julia/LWFBinput/Pfyn_irrigation_ambient/pfynwald_meteoveg.csv", skip=2, col_names=met_names, show_col_types=F)
metveg = read_csv("../Julia/LWFBinput/Pfyn_irr_stop/pfynwald_meteoveg.csv", skip=2, col_names=met_names, show_col_types=F)

met = metveg[, 1:7]
veg = metveg[, 8:11]

met$prec_mmDay = meteo_Con$precip_ctrl
met$irrig_mmDay = meteo_Con$irrig_mm
#met$irrig_mmDay = if_else(met$dates <= "2014-01-01", meteo_Con$irrig_mm, 0) # irr stop

metirr = cbind(met, veg)

write_csv(metirr, "../Julia/LWFBinput/Pfyn_irrigation_test/pfynwald_meteoveg.csv")


# isotopes

iso = read_csv("../../Data/Pfyn/Rhone_iso2000_23.csv")
iso$date = as.Date(iso$DateTime)

iso_wide = iso %>% select(-DateTime) %>% pivot_wider(names_from=Symbol, values_from=Amount)

iso_wide = rbind(iso_wide, data.frame(date="2024-01-01", O18=median(iso_wide$O18), H2=median(iso_wide$H2)))

write_csv(iso_wide, "../Julia/LWFBinput/Pfyn_irrigation_test/pfynwald_irrigiso.csv")

# separate treatments

# control
meteo_cont = meteo_Con %>% rename(prec = precip_ctrl) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# irrigation stop
meteo_irrstp = meteo_Con %>% rename(prec = precip_irrstp) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# irrigation under ambient conditions
meteo_irr_con = meteo_Con %>% rename(prec = precip_ctrl, irrig=irrig_mm) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec, irrig)

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
LAI_df = LAI_df[1:25,]

# soil data

soil_df = read_csv("../../Data/Pfyn/soil_hydraulic.csv")


## use LWFBrook90R to prepare input files for julia version

# options common to all scenarios
opt = set_optionsLWFB90(startdate=as.Date("2000-01-01"), enddate=as.Date("2024-12-31"), 
                        root_method="soilvar", budburst_method="Menzel", 
                        leaffall_method="vonWilpert", lai_method="Coupmodel")

# parameters for LAI and climate vary by scenario

# control scenario
par_c = set_paramLWFB90(maxlai=LAI_df$LAI_ctrl, winlaifrac=.6, height=ht_cont, height_ini=ht_cont, 
                        coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, bypar=1, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_control","pfynwald",".", options_b90=opt, param_b90=par_c, climate=filter(meteo_cont, dates<"2025-01-01"), soil=soil_df)
#extend_meteoveg(meteo_cont, "Pfyn_control")

# drought scenarios use same parameters as control

# drought scenario under ambient climate
generate_LWFBrook90jl_Input("Pfyn_drought_ambient","pfynwald",".", options_b90=opt, param_b90=par_c, climate=filter(meteo_roof_con, dates<"2025-01-01"), soil=soil_df)
#extend_meteoveg(meteo_roof_con, "Pfyn_drought_ambient")

# drought scenario under manipulated VPD
generate_LWFBrook90jl_Input("Pfyn_drought_VPD","pfynwald",".", options_b90=opt, param_b90=par_c, climate=filter(meteo_roof_vpd, dates<"2025-01-01"), soil=soil_df)
#extend_meteoveg(meteo_roof_vpd, "Pfyn_drought_VPD")

# irrigation stop scenario
par_irst = set_paramLWFB90(maxlai=LAI_df$LAI_irrstp, winlaifrac=.6, height=ht_stp, height_ini=ht_stp, 
                           coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, bypar=1, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_irr_stop","pfynwald",".", options_b90=opt, param_b90=par_irst, climate=filter(meteo_irrstp, dates<"2025-01-01"), soil=soil_df)
#extend_meteoveg(meteo_irrstp, "Pfyn_irr_stop")

# irrigation scenarios use same parameters
par_ir = set_paramLWFB90(maxlai=LAI_df$LAI_irr, winlaifrac=.6, height=ht_irr, height_ini=ht_irr, 
                         coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, bypar=1, budburst_species="Pinus sylvestris")

# irrigation under ambient climate
generate_LWFBrook90jl_Input("Pfyn_irrigation_test","pfynwald",".", options_b90=opt, param_b90=par_ir, climate=filter(meteo_irr_con, dates<"2025-01-01"), soil=soil_df)
#extend_meteoveg(meteo_irr_con, "Pfyn_irrigation_test")

# irrigation under manipulated VPD
generate_LWFBrook90jl_Input("Pfyn_irrigation_VPD","pfynwald",".", options_b90=opt, param_b90=par_ir, climate=filter(meteo_irr_vpd, dates<"2025-01-01"), soil=soil_df)
#extend_meteoveg(meteo_irr_vpd, "Pfyn_irrigation_VPD")
