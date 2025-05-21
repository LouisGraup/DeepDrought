library(rLWFpg)
library(tidyverse)
library(lubridate)
library(plantecophys)

## retrieve in-situ data from Pfynwald

con <- db_connect(username = "grauplou", password = rstudioapi::askForPassword("enter password"), db_host = "pgdblwf", db_name = "lwf")

# Information about measurements
messvar.df <- db_tbl(con, schema = "ada", table = "v_messvar", retrieve = TRUE) %>%
  dplyr::filter( project_name %in% 'Pfynwald irrigation', varfreq_name != "1 Minutes")

messvar.tbl <- db_tbl(con, schema = "ada", table = "v_messvar") %>%
  dplyr::filter( project_name %in% 'Pfynwald irrigation')

pfyndata.tbl = db_tbl(con, table="pfyn_messdat")

#Messvar_ID Information
#134 -> Atmospheric precipitation weighing bucket, 10 min
#160 -> Air relative humidity 2 m, control, 10 min
#158 -> Air temperature 2 m, control, 10 min
#162 -> Atmospheric vapor pressure deficit 2 m, 10 min
#133 -> Photosynthetic photon flux density, 10 min
#131 -> Radiation shortwave incoming over canopy, 10 min
#125 -> Wind direction over canopy, 10 min
#124 -> Wind speed 2 m, 10 min

#62 -> Soil water matric potential, 10 min (T21_SMP)
#63 -> Soil temperature, 10 min (T21_SMP)
#64 -> Stem sapflow, 10 min
#66 -> Stem temperature, 10 min
#76 -> Soil water matric potential, 20 min (TM_SMP)
#77 -> Soil temperature, 20 min (TM_SMP)
#79 -> Standard water temperature, 10 min


# meteo data
meteo.df = messvar.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  filter(messvar_id %in% c(124,131,134,158,162)) %>% 
  select(messvar_id, messtime, messval, varname_name, messvar_name, varunit_symbol, varvpos, treatment) %>% 
  collect()

meteo.wide = meteo.df %>% pivot_wider(id_cols = messtime, names_from=messvar_name, values_from=messval)
meteo.wide$WS_2m = ifelse(meteo.wide$WS_2m>2, 1, meteo.wide$WS_2m)

# aggregate to daily
meteo.wide$dates = as.Date(meteo.wide$messtime)
meteo.day = meteo.wide %>% group_by(dates) %>% 
  summarize(globrad_MJDayM2=mean(GlobalRad_oc_Avg)*.0864, # convert W/m2 to MJ/day/m2
            tmax_degC=max(TA_2m_Avg, na.rm=T), tmin_degC=min(TA_2m_Avg, na.rm=T),
            vappres_kPa=mean(VPD_2m_Avg, na.rm=T), windspeed_ms=mean(WS_2m, na.rm=T),
            prec_mmDay=sum(RaineH3_amount_dif_Tot, na.rm=T))

dateseq = seq.Date(as.Date("2021-06-07"), Sys.Date()-1, by=1)

meteo = data.frame(dates=dateseq) %>% left_join(meteo.day)



## combine with MeteoSwiss data

# convert MeteoSwiss to LWFBrook90 input..............................................
#gre000d0  W/m2                                 Global radiation; daily mean
#tre200dn  ?C                                   Air temperature 2 m above ground; daily minimum
#tre200d0  ?C                                   Air temperature 2 m above ground; daily mean
#tre200jx  ?C                                   Air temperature 2 m above ground; daily maximum (6 UTC up to 18 UTC)
#rka150d0  mm                                   Precipitation; daily total 0 UTC - 0 UTC
#ure200d0  %                                    Relative air humidity; 2 m above ground; daily mean
#pva200d0  hP                                   Vapor pressure 2 m above ground; daily mean
#fkl010d0  m/s                                  Wind speed scalar; daily mean

# Sion station (SIO) lat. 46° 13' 07'' lon. 07° 19' 49'' elev. 482 m
meteoCH = read.table("../../Data/MeteoSwiss/Sion/order_129912_data.txt", sep=";", header=T, skip=2)

meteoCH2<-meteoCH%>%
  mutate( dates = as.Date(as.character(time), format="%Y%m%d"),
          tmin_degC = tre200dn,
          tmax_degC = tre200jx,
          tmean = tre200d0,
          prec_mmDay = rka150d0,
          relhum = ure200d0,
          globrad_MJDayM2 = ((24*60*60)/1000000)*gre000d0, # convert from W/m2 to MJ/day/m2
          windspeed_ms = fkl010d0,
          vappres_kPa = RHtoVPD(relhum, tmean, Pa = 101)) %>%
  select(dates, globrad_MJDayM2, tmax_degC, tmin_degC, vappres_kPa, windspeed_ms, prec_mmDay)

meteo_comp = inner_join(meteoCH2, meteo, by="dates")

ggplot(meteo_comp, aes(tmax_degC.x, tmax_degC.y))+geom_point()
ggplot(meteo_comp, aes(tmin_degC.x, tmin_degC.y))+geom_point()
ggplot(meteo_comp, aes(vappres_kPa.x, vappres_kPa.y))+geom_point()
ggplot(meteo_comp, aes(prec_mmDay.x, prec_mmDay.y))+geom_point()
ggplot(meteo_comp, aes(globrad_MJDayM2.x, globrad_MJDayM2.y))+geom_point()
ggplot(meteo_comp, aes(windspeed_ms.x, windspeed_ms.y))+geom_point() # consider differences

