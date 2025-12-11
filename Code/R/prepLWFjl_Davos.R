
## use helper functions to create input files for LWFBrook90.jl

library(tidyverse)
library(lubridate)

## read in meteo data

meteo = read.csv("../../Data/Davos/meteo/Davos_meteo.csv")
meteo$dates = as.Date(meteo$date)

meteo = meteo %>% select(-ET, -date)

soil_df = read.csv("../../Data/Pfyn/soil_hydraulic.csv") # dummy filler

## use LWFBrook90R to prepare input files for julia version

library(LWFBrook90R)
source("C:/Users/grauplou/Documents/LWFBrook90.jl/src/generate_LWFBrook90jl_Input_mod.R", echo=F)

opt = set_optionsLWFB90(startdate=as.Date("2007-01-01"), enddate=as.Date("2022-12-31"), 
                        root_method="soilvar", budburst_method="Menzel", 
                        leaffall_method="vonWilpert", lai_method="Coupmodel")

# control scenario
par = set_paramLWFB90(maxlai=3, winlaifrac=.6, height=25, height_ini=25, 
                      coords_y=47, budburst_species="Picea abies (frueh)")
generate_LWFBrook90jl_Input("Davos","davos",".", options_b90=opt, param_b90=par, climate=meteo, soil=soil_df)

