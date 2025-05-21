
## use helper functions to create input files for LWFBrook90.jl

library(tidyverse)
library(lubridate)
library(plantecophys)

# read in meteo data
# data from MeteoSwiss stations Sion and Sierre (precip post-2013)
# combined with irrigation data in Pfynwald

meteo = read.csv("../../Data/Pfyn/meteo_irr.csv")
meteo$dates = as.Date(meteo$dates)

meteo_cont = meteo %>% rename(prec = precip_ctrl) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)
meteo_irr = meteo %>% rename(prec = precip_irr) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)
meteo_irrstp = meteo %>% rename(prec = precip_irrstp) %>% 
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec)

# yearly sum of precip for lai comparison
meteoy = mutate(meteo_cont, year=year(dates))
meteoy = meteoy %>% group_by(year) %>% summarize_at(vars(prec), list(sum))

# look back at previous years
meteoy = meteoy %>% mutate(prev_prec = 0, prev2_prec = 0)
meteoy$prev_prec[2:25] = meteoy$prec[1:24]
meteoy$prev2_prec[3:25] = meteoy$prec[1:23] + meteoy$prec[2:24]

# read in site data

site_df = read.csv("../../Data/Pfyn/siteproperties.csv")
ht_cont = mean(filter(site_df, BEZKM == "control")$height_m)
ht_irr = mean(filter(site_df, BEZKM == "irrigation")$height_m)
ht_stp = mean(filter(site_df, BEZKM == "stop")$height_m)

# read in lai data

lai_df = read.csv("../../Data/Pfyn/PFY_lai.csv")
lai_cont = filter(lai_df, trt=="control") %>% group_by(year) %>% 
  summarize_at(vars(LAI), list(mean))
lai_irr = filter(lai_df, trt=="irrig.") %>% group_by(year) %>% 
  summarize_at(vars(LAI), list(mean))

# control plots

lai_lm_c = lm(LAI ~ poly(year, 2, raw=TRUE), lai_cont)
lai_fit_c = predict(lai_lm_c)
lai_cont$pred = lai_fit_c

ggplot(lai_cont, aes(year, LAI))+geom_point()+
  geom_line(aes(year, pred), color="red")

lai_ext <- data.frame(year=c(2000:2020))
lai_ext$lai_pred<-predict(lai_lm_c, newdata=lai_ext,type="response")
lai_ext$lai_pred[1:5] = lai_cont$LAI[1]
#lai_ext = left_join(lai_ext, lai_cont)
#lai_ext$LAI[is.na(lai_ext$LAI)] = lai_ext$lai_pred

#ggplot(lai_ext, aes(year, LAI))+geom_point()

# quick analysis against annual precipitation
lai_cont = left_join(lai_cont, meteoy)
plot(lai_cont$prec, lai_cont$LAI)
plot(lai_cont$prev_prec, lai_cont$LAI)
plot(lai_cont$prev2_prec, lai_cont$LAI)

# none are very informative
summary(lm(LAI ~ prec, lai_cont))
summary(lm(LAI ~ prev_prec, lai_cont))
summary(lm(LAI ~ prev2_prec, lai_cont))


# irrigated plots

lai_lm_i = lm(LAI ~ poly(year, 2, raw=TRUE), lai_irr)
lai_fit_i = predict(lai_lm_i)
lai_irr$pred = lai_fit_i

ggplot(lai_irr, aes(year, LAI))+geom_point()+
  geom_line(aes(year, pred), color="red")

lai_ext_irr <- data.frame(year=c(2000:2020))
lai_ext_irr$lai_pred<-predict(lai_lm_i, newdata=lai_ext_irr,type="response")
lai_ext_irr$lai_pred[1:4] = lai_cont$LAI[1]
lai_ext_irr = left_join(lai_ext_irr, lai_irr)
lai_ext_irr$LAI[is.na(lai_ext_irr$LAI)] = lai_ext_irr$lai_pred[is.na(lai_ext_irr$LAI)]

ggplot(lai_ext_irr, aes(year, LAI))+geom_point()


# irrigation stop plots
# use irrigated LAI until 2014, then average of irrigated and control

lai_ext_irrstp = select(lai_ext_irr, year, LAI) %>% rename(LAI_irr=LAI)
lai_ext_irrstp$LAI_irr_c_avg = (lai_ext_irrstp$LAI_irr + lai_ext$lai_pred) / 2
lai_ext_irrstp$LAI_irrstp = with(lai_ext_irrstp, ifelse(year<2014, LAI_irr, LAI_irr_c_avg))


# read in soil and root data

soil_df = read.csv("../../Data/Pfyn/soil_hydraulic.csv")


library(LWFBrook90R)
source("C:/Users/grauplou/Documents/LWFBrook90.jl/src/generate_LWFBrook90jl_Input_mod.R", echo=F)

opt = set_optionsLWFB90(startdate=as.Date("2000-01-01"), enddate=as.Date("2024-12-31"), 
                        root_method="soilvar", budburst_method="Menzel", 
                        leaffall_method="vonWilpert", lai_method="Coupmodel")

# control scenario
par_c = set_paramLWFB90(maxlai=lai_ext$lai_pred, winlaifrac=.6, height=ht_cont, height_ini=ht_cont, 
                      coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_control","pfynwald",".", options_b90=opt, param_b90=par_c, climate=meteo_cont, soil=soil_df)

# irrigation scenario
par_ir = set_paramLWFB90(maxlai=lai_ext_irr$LAI, winlaifrac=.6, height=ht_irr, height_ini=ht_irr, 
                        coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_irrigation","pfynwald",".", options_b90=opt, param_b90=par_ir, climate=meteo_irr, soil=soil_df)

# irrigation stop scenario
par_irst = set_paramLWFB90(maxlai=lai_ext_irrstp$LAI_irrstp, winlaifrac=.6, height=ht_stp, height_ini=ht_stp, 
                         coords_x=7.611329, coords_y=46.301624, eslope=7.6, aspect=299, budburst_species="Pinus sylvestris")
generate_LWFBrook90jl_Input("Pfyn_irr_stop","pfynwald",".", options_b90=opt, param_b90=par_irst, climate=meteo_irrstp, soil=soil_df)



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

