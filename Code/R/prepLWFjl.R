
## use helper functions to create input files for LWFBrook90.jl

library(tidyverse)
library(lubridate)
library(plantecophys)

# read in meteo data

meteoCH = read.table("./Data/MeteoSchweiz/order_129912_data.txt", sep=";", header=T, skip=2)
# data frame for LWFBrook90R
meteoLWF<-meteoCH%>%
  mutate( dates = as.Date(as.character(time), format="%Y%m%d"),
          tmin = tre200dn,
          tmax = tre200jx,
          tmean = tre200d0,
          prec = rka150d0,
          relhum = ure200d0,
          globrad = ((24*60*60)/1000000)*gre000d0, # convert from W/m2 to MJ/day/m2
          windspeed = fkl010d0,
          vappres = RHtoVPD(relhum, tmean, Pa = 101)) %>%
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)
meteoLWF = filter(meteoLWF, dates<as.Date("2025-01-01"))

# read in site data

site_df = read.csv("./Data/Pfyn/siteproperties.csv")
ht_cont = mean(filter(site_df, BEZKM == "control")$height_m)

# read in lai data

lai_df = read.csv("./Data/Pfyn/PFY_lai.csv")
lai_cont = mean(filter(lai_df, year==2004, trt=="control")$LAI)

# read in soil and root data

soil_df = read.csv("./Data/Pfyn/soil_hydraulic.csv")


library(LWFBrook90R)
source("~/Documents/WSL/DeepDrought/LWFBrook90.jl/src/generate_LWFBrook90jl_Input.R", echo=F)

opt = set_optionsLWFB90(startdate=as.Date("2000-01-01"), enddate=as.Date("2000-12-31"), 
                        root_method="soilvar", budburst_method="Menzel", 
                        leaffall_method="vonWilpert", lai_method="Coupmodel")
par = set_paramLWFB90(maxlai=lai_cont, winlaifrac=.6, height=ht_cont, height_ini=ht_cont, 
                      coords_x=7.330278, coords_y=46.218611, budburst_species="Pinus sylvestris")

generate_LWFBrook90jl_Input("LWFBinput","pfynwald",".", options_b90=opt, param_b90=par, climate=meteoLWF, soil=soil_df)

# test run Pfynwald in R
par = set_paramLWFB90(maxlai=lai_cont, winlaifrac=.6, height=ht_cont, height_ini=ht_cont, 
                      coords_x=7.330278, coords_y=46.218611, budburst_species="Pinus sylvestris",
                      maxrootdepth=-1, frintlai=0, frintsai=0, fsintlai=0, fsintsai=0)

res <- run_LWFB90(options_b90=opt, param_b90=par, climate=meteoLWF, soil=soil_df)
out = res$output

summary(out)

ggplot(out, aes(doy, rfal/100), fill="black")+geom_col()+geom_line(aes(doy, cumsum(balerr), color="WB"))

ggplot(out, aes(doy, irvp, color="IRVP"))+geom_line()+geom_line(aes(doy, intr, color="IntR"))+
  geom_line(aes(doy, rint, color="RInt"))

wb_model = select(out, doy, rfal, sfal, evap, seep, flow, swat, snow, gwat, intr, ints)
wb_model = mutate(wb_model, inflow=rfal+sfal, outflow=evap+seep+flow, 
                  swat_diff=c(0,diff(swat)), snodiff=c(0,diff(snow)),
                  intr_diff=c(0,diff(intr)), ints_diff=c(0,diff(ints)),
                  gwat_diff=c(0,diff(gwat)))
wb_model$wb_error = with(wb_model, inflow - outflow - swat_diff - snodiff - intr_diff - ints_diff - gwat_diff)
wb_model$wb_error[1] = 0
wb_model$wb_cs = cumsum(wb_model$wb_error)

ggplot(wb_model, aes(doy, rfal), fill="black")+geom_col()+geom_line(aes(doy, evap, color="ET"))+
  geom_line(aes(doy, seep, color="Seep"))+geom_line(aes(doy, flow, color="Flow"))+
  geom_line(aes(doy, wb_cs, color="WB"))

ggplot(wb_model, aes(doy, swat/100, color="SWAT"))+geom_line()+geom_line(aes(doy, snow, color="Snow"))+
  geom_line(aes(doy, gwat, color="GWAT"))+geom_line(aes(doy, intr, color="IntR"))+geom_line(aes(doy, ints, color="IntS"))+
  geom_line(aes(doy, wb_cs, color="WB"))

