# combine meteo data with irrigation from database

library(rLWFpg)
library(tidyverse)
library(lubridate)
library(humidity)

con <- db_connect(username = "grauplou", password = rstudioapi::askForPassword("enter password"), db_host = "pgdblwf", db_name = "lwf")

messvar.tbl <- db_tbl(con, schema = "ada", table = "v_messvar") %>%
  dplyr::filter( project_name %in% 'Pfynwald irrigation')

# View(collect(messvar.tbl))

pfyndata.tbl = db_tbl(con, table="pfyn_messdat")

## irrigation amounts from database
irr.df = messvar.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  filter(variable_id %in% c(263), sensor_id != 1) %>% 
  select(installation_name, messtime, messval, messvar_name, messvar_id) %>% 
  collect()
irr.df$dates = as.Date(irr.df$messtime)

# further filter out manual amounts and remove extra columns
irr.df = irr.df %>% filter(messvar_id %in% 8099:8102) %>% 
  select(-c("messtime","messvar_id", "messvar_name")) %>% rename(irrig_mm = messval)

# sum up multiple amounts on single days for each plot
irr.df = irr.df %>% group_by(dates, installation_name) %>% summarize_at(vars(irrig_mm), list(sum))

# average over treatment plots
irr = irr.df %>% group_by(dates) %>% summarize_at(vars(irrig_mm), list(mean))

dbDisconnect(con)

## Sierre and Sion precipitation from MeteoSwiss
## Sierre has automatic and historical manual measurements

Sion = read.table("../../Data/MeteoSwiss/Sion/order_129912_data.txt", sep=";", header=T, skip=2)
Sierre = read.table("../../Data/MeteoSwiss/Sierre/order_130615_data.txt", sep=";", header=T, skip=2)
Sierre_hist = read.table("../../Data/MeteoSwiss/Sierre/ogd-nime_sre_d_historical.csv", sep=";", header=T)

Sion$dates = as.Date(as.character(Sion$time), format="%Y%m%d")
Sierre$dates = as.Date(as.character(Sierre$time), format="%Y%m%d")
Sierre_hist$dates = as.Date(Sierre_hist$reference_timestamp, format="%d.%m.%Y %H:%M")

# remove NA values
Sierre = Sierre %>% rename(precip=rka150d0) %>% filter(precip != "-") %>% select(dates, precip)
Sierre$precip = as.numeric(Sierre$precip)

Sierre_hist = Sierre_hist %>% rename(precip=rre150d0, stn=station_abbr) %>% 
  select(dates, precip)

## first, Sion data for long-term meteo inputs

meteoLWF = Sion %>%
  mutate( tmin = tre200dn,
          tmax = tre200jx,
          tmean = tre200d0,
          prec_Sion = rka150d0,
          relhum = ure200d0,
          globrad = ((24*60*60)/1000000)*gre000d0, # convert from W/m2 to MJ/day/m2
          windspeed = fkl010d0) %>%
  select(dates, globrad, tmax, tmin, tmean, windspeed, relhum, prec_Sion) %>% 
  filter(dates<as.Date("2024-01-01"))

# derive actual vapor pressure
meteoLWF$Es = SVP.ClaCla(meteoLWF$tmean+273.15) # saturation vapor pressure hPa
meteoLWF$Esmean = ( SVP.ClaCla(meteoLWF$tmax+273.15) + SVP.ClaCla(meteoLWF$tmin+273.15) ) / 2
meteoLWF$vappres = WVP2(meteoLWF$relhum, meteoLWF$Esmean) / 1000 # actual vapor pressure kPa

meteoLWF = select(meteoLWF, -c(Es, Esmean, relhum))


## then, Sierre precipitation (manual until 2014, then automatic) until 2024

meteoLWF = left_join(meteoLWF, rename(Sierre_hist, precip_Sierre_hist=precip))
meteoLWF = left_join(meteoLWF, rename(Sierre, precip_Sierre=precip))

# fill in missing values in Sierre column with Sion data
meteoLWF$precip_ctrl = ifelse(meteoLWF$dates<"2014-01-01", meteoLWF$precip_Sierre_hist, meteoLWF$precip_Sierre)
meteoLWF$precip_ctrl[is.na(meteoLWF$precip_ctrl)] = meteoLWF$prec_Sion[is.na(meteoLWF$precip_ctrl)]

meteoLWF = select(meteoLWF, -c(prec_Sion, precip_Sierre, precip_Sierre_hist))


## last, in-situ data from 2024

meteo_ins_VPD = read_csv("../../Data/Pfyn/meteo/Pfyn_meteo0124_0725_insitu_VPD.csv")
meteo_ins_Con = read_csv("../../Data/Pfyn/meteo/Pfyn_meteo0124_0725_insitu_Control.csv")

# VPD is reduced VPD from manipulation experiment, also includes changes in temp
meteoLWF_Con = rbind(meteoLWF, rename(meteo_ins_Con, precip_ctrl=precip))
meteoLWF_VPD = rbind(meteoLWF, rename(meteo_ins_VPD, precip_ctrl=precip))

# assume 50% reduction in P from roof for additional soil treatments
meteoLWF_Con$precip_roof = with(meteoLWF_Con, if_else(dates >= "2024-01-01", precip_ctrl * .5, precip_ctrl))
meteoLWF_VPD$precip_roof = with(meteoLWF_VPD, if_else(dates >= "2024-01-01", precip_ctrl * .5, precip_ctrl))


## now using precip column to combine with irrigation
meteoLWF_Con = left_join(meteoLWF_Con, irr)
meteoLWF_VPD = left_join(meteoLWF_VPD, irr)

meteoLWF_Con$irrig_mm[is.na(meteoLWF_Con$irrig_mm)] = 0
meteoLWF_VPD$irrig_mm[is.na(meteoLWF_VPD$irrig_mm)] = 0

meteoLWF_Con$precip_irr = meteoLWF_Con$precip_ctrl + meteoLWF_Con$irrig_mm
meteoLWF_VPD$precip_irr = meteoLWF_VPD$precip_ctrl + meteoLWF_VPD$irrig_mm

# irrigation stop experiment only applies to atmospheric control treatment
meteoLWF_Con$precip_irrstp = with(meteoLWF_Con, if_else(dates < "2014-01-01", precip_irr, precip_ctrl))

#write_csv(meteoLWF_Con, "../../Data/Pfyn/meteo/meteo_irr_Control.csv")
#write_csv(meteoLWF_VPD, "../../Data/Pfyn/meteo/meteo_irr_VPD.csv")

## plot irrigation amounts
meteo_plot = meteoLWF_Con %>% filter(dates >= "2003-01-01") %>% select(-c(2:7,9)) %>% pivot_longer(-dates, names_to="treatment", values_to="mm")

# daily
ggplot(meteo_plot, aes(dates, mm, color=treatment)) + geom_line() + facet_wrap(~treatment, ncol=1)

meteo_plot$year = year(meteo_plot$dates)
meteo_plot$month = month(meteo_plot$dates)
meteo_plot2 = meteo_plot %>% group_by(year, month, treatment) %>% summarize_at(vars(mm), list(sum))
meteo_plot2$treatment = factor(meteo_plot2$treatment, levels=c("precip_irr", "irrig_mm", "precip_irrstp", "precip_ctrl"))
meteo_plot2$date = as.Date(paste(meteo_plot2$year, meteo_plot2$month, 1, sep="-"))

# monthly
ggplot(filter(meteo_plot2, treatment %in% c("precip_ctrl", "irrig_mm")), aes(x=date, y=mm, fill=treatment))+geom_bar(position="stack", stat="identity")

# yearly comparison
meteo_yr = meteo_plot %>% group_by(year, treatment) %>% summarize_at(vars(mm), list(sum)) 
meteo_yr_pre14 = meteo_yr %>% filter(year<2014) %>% group_by(treatment) %>% summarize_at(vars(mm), list(mean))
meteo_yr_post14 = meteo_yr %>% filter(year>=2014) %>% group_by(treatment) %>% summarize_at(vars(mm), list(mean))
meteo_yr = meteo_yr %>% group_by(treatment) %>% summarize_at(vars(mm), list(mean))

meteo_yr

# get dates
meteoLWF_Con$irrig_dates = ifelse(meteoLWF_Con$irrig_mm > 0, meteoLWF_Con$dates, NA)
meteoLWF_Con$year = year(meteoLWF_Con$dates)
meteo_irrig = select(meteoLWF_Con, irrig_dates, year)
meteo_irrig = na.omit(meteo_irrig)
met_irrig = meteo_irrig %>% group_by(year) %>% 
  summarize(on=min(irrig_dates), off=max(irrig_dates))
met_irrig$on = as.Date(met_irrig$on)
met_irrig$off = as.Date(met_irrig$off)

on = c("2003-06-19", "2004-05-15", "2005-04-23", "2006-05-06", "2007-05-04", "2008-05-15", 
       "2009-05-14", "2010-06-17", "2011-05-14", "2012-05-11", "2013-05-17", "2014-05-19",
       "2015-05-12", "2016-05-30", "2017-05-16", "2018-05-08", "2019-05-17", "2020-05-26",
       "2021-05-18", "2022-05-11", "2023-05-15", "2024-05-26", "2025-04-30")

off = c("2003-10-21", "2004-10-26", "2005-10-04", "2006-10-25", "2007-10-02", "2008-10-14",
        "2009-10-12", "2010-10-01", "2011-10-16", "2012-10-02", "2013-09-23", "2014-10-01",
        "2015-10-05", "2016-09-26", "2017-10-09", "2018-09-27", "2019-10-15", "2020-10-19",
        "2021-10-11", "2022-10-21", "2023-10-18", "2024-10-24", "2025-07-31")

irr = data.frame(on=as.Date(on), off=as.Date(off))
irr$year = year(irr$on)