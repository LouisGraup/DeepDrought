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
  filter(messvar_id %in% c(124,131,134,158,160,162)) %>% 
  select(messvar_id, messtime, messval, varname_name, messvar_name, varunit_symbol, varvpos, treatment) %>% 
  collect()

meteo.wide = meteo.df %>% pivot_wider(id_cols = messtime, names_from=messvar_name, values_from=messval)
meteo.wide$WS_2m = ifelse(meteo.wide$WS_2m>2, 1, meteo.wide$WS_2m)

# aggregate to daily
meteo.wide$dates = as.Date(meteo.wide$messtime)
meteo.day = meteo.wide %>% group_by(dates) %>% 
  summarize(globrad=mean(GlobalRad_oc_Avg)*.0864, # convert W/m2 to MJ/day/m2
            tmax=max(TA_2m_Avg, na.rm=T), tmin=min(TA_2m_Avg, na.rm=T),
            VPD=mean(VPD_2m_Avg, na.rm=T), RH=mean(RH_2m_Avg, na.rm=T),
            windspeed=mean(WS_2m, na.rm=T), prec=sum(RaineH3_amount_dif_Tot, na.rm=T))

dateseq = seq.Date(as.Date("2021-06-07"), Sys.Date()-1, by=1)

meteo_Pfyn = data.frame(dates=dateseq) %>% left_join(meteo.day)


## compare with MeteoSwiss data

# MeteoSwiss metadata
#gre000d0  W/m2                                 Global radiation; daily mean
#tre200dn  ?C                                   Air temperature 2 m above ground; daily minimum
#tre200d0  ?C                                   Air temperature 2 m above ground; daily mean
#tre200jx  ?C                                   Air temperature 2 m above ground; daily maximum (6 UTC up to 18 UTC)
#rka150d0  mm                                   Precipitation; daily total 0 UTC - 0 UTC
#ure200d0  %                                    Relative air humidity; 2 m above ground; daily mean
#pva200d0  hP                                   Vapor pressure 2 m above ground; daily mean
#fkl010d0  m/s                                  Wind speed scalar; daily mean

# Sion station (SIO) lat. 46° 13' 07'' lon. 07° 19' 49'' elev. 482 m
meteoSion = read.table("../../Data/MeteoSwiss/Sion/order_129912_data.txt", sep=";", header=T, skip=2)

meteoSion<-meteoSion %>%
  mutate( dates = as.Date(as.character(time), format="%Y%m%d"),
          tmin = tre200dn,
          tmax = tre200jx,
          tmean = tre200d0,
          prec = rka150d0,
          RH = ure200d0,
          globrad = ((24*60*60)/1000000)*gre000d0, # convert from W/m2 to MJ/day/m2
          windspeed = fkl010d0,
          VPD = RHtoVPD(RH, tmean, Pa = 101)) %>%
  select(dates, globrad, tmax, tmin, RH, VPD, windspeed, prec, stn)

meteo_comp = inner_join(meteoSion, meteo_Pfyn, by="dates")

ggplot(meteo_comp, aes(tmax.x, tmax.y))+geom_point()
ggplot(meteo_comp, aes(tmin.x, tmin.y))+geom_point()
ggplot(meteo_comp, aes(VPD.x, VPD.y))+geom_point()
ggplot(meteo_comp, aes(RH.x, RH.y))+geom_point()
ggplot(meteo_comp, aes(prec.x, prec.y))+geom_point()
ggplot(meteo_comp, aes(globrad.x, globrad.y))+geom_point()
ggplot(meteo_comp, aes(windspeed.x, windspeed.y))+geom_point() # consider differences



## compare Sierre and Sion precipitation from MeteoSwiss
## as well as Sierre automatic and historical manual measurements

Sierre = read.table("../../Data/MeteoSwiss/Sierre/order_130615_data.txt", sep=";", header=T, skip=2)
Sierre_hist = read.table("../../Data/MeteoSwiss/Sierre/ogd-nime_sre_d_historical.csv", sep=";", header=T)

Sierre$dates = as.Date(as.character(Sierre$time), format="%Y%m%d")
Sierre_hist$dates = as.Date(Sierre_hist$reference_timestamp, format="%d.%m.%Y %H:%M")

# remove NA values
Sierre = Sierre %>% rename(precip=rka150d0) %>% filter(precip != "-") %>% 
  select(dates, precip, stn)
Sierre$precip = as.numeric(Sierre$precip)

Sierre_hist = Sierre_hist %>% rename(precip=rre150d0, stn=station_abbr) %>% 
  select(dates, precip, stn)

# first compare Sion and Sierre historical
meteo_comp_ext = inner_join(Sierre_hist, select(meteoSion, dates, prec))
ggplot(meteo_comp_ext, aes(prec, precip))+geom_point()+labs(x="Sion",y="Sierre")

meteo_comp_ext$year = year(meteo_comp_ext$dates)
meteo_comp_ext_yr = meteo_comp_ext %>% group_by(year) %>% summarize_at(vars(precip, prec), list(sum))
t.test(meteo_comp_ext_yr$precip, meteo_comp_ext_yr$prec)

meteo_comp_ext_year = meteo_comp_ext_yr %>% pivot_longer(-year)
ggplot(meteo_comp_ext_year, aes(year, value, fill=name))+geom_col(position="dodge")

# match time period to Sierre station
Sion_comp = meteoSion %>% filter(dates %in% Sierre$dates) %>% rename(precip=prec) %>% 
  select(dates, precip, stn)

Sierre_hist_comp = Sierre_hist %>% filter(dates %in% Sierre$dates)

meteo_comp = rbind(Sion_comp, Sierre)
meteo_comp = rbind(meteo_comp, Sierre_hist_comp)
meteo_comp$year = year(meteo_comp$dates)

meteo_comp = filter(meteo_comp, year>2013)

# compare yearly sums
meteo_yr = meteo_comp %>% group_by(year, stn) %>% summarize_at(vars(precip), list(sum))
ggplot(meteo_yr, aes(year, precip, fill=stn))+geom_col(position="dodge")

meteo_yr2 = meteo_yr %>% pivot_wider(names_from=stn, values_from=precip)
t.test(meteo_yr2$SIO, meteo_yr2$VSSIE) # yearly means are similar
t.test(meteo_yr2$SRE, meteo_yr2$VSSIE) # yearly means are similar between Sierre stations

# compare daily values
meteo_wide = meteo_comp %>% pivot_wider(names_from=stn, values_from=precip)
ggplot(meteo_wide, aes(VSSIE, SIO))+geom_point()
ggplot(meteo_wide, aes(VSSIE, SRE))+geom_point()
ggplot(meteo_wide, aes(SIO, SRE))+geom_point()
t.test(meteo_wide$SIO, meteo_wide$VSSIE) # daily means are similar
t.test(meteo_wide$SRE, meteo_wide$VSSIE) # daily means are similar between Sierre stations


## compare to in-situ Pfynwald data

Sion_Pcomp = Sion_comp %>% select(dates, precip) %>% rename(precip_Sion=precip)
Sierre_Pcomp = Sierre %>% select(dates, precip) %>% rename(precip_Sierre=precip)
Sierre_hist_Pcomp = Sierre_hist_comp %>% select(dates, precip) %>% rename(precip_Sierre_hist=precip)

P_comp = left_join(meteo_Pfyn, Sion_Pcomp)
P_comp = left_join(P_comp, Sierre_Pcomp)
P_comp = left_join(P_comp, Sierre_hist_Pcomp)

# compare daily values
ggplot(P_comp, aes(prec, precip_Sion))+geom_point()
ggplot(P_comp, aes(prec, precip_Sierre))+geom_point()
ggplot(P_comp, aes(prec, precip_Sierre_hist))+geom_point()

# compare yearly values
P_comp$year = year(P_comp$dates)
P_yr = P_comp %>% group_by(year) %>% summarize_at(vars(prec, precip_Sion, precip_Sierre, precip_Sierre_hist), list(sum))
P_yr_long = P_yr %>% pivot_longer(-year) %>% filter(year<2025)

ggplot(P_yr_long, aes(year, value, fill=name))+geom_col(position="dodge")
t.test(P_yr$prec, P_yr$precip_Sion)
t.test(P_yr$prec, P_yr$precip_Sierre)
t.test(P_yr$prec, P_yr$precip_Sierre_hist)