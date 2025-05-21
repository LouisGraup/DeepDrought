# combine meteo data with irrigation from database

library(rLWFpg)
library(tidyverse)
library(lubridate)
library(plantecophys)

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


## compare Sierre and Sion precipitation from MeteoSwiss

Sion = read.table("../../Data/MeteoSwiss/Sion/order_129912_data.txt", sep=";", header=T, skip=2)
Sierre = read.table("../../Data/MeteoSwiss/Sierre/order_130615_data.txt", sep=";", header=T, skip=2)

Sion$dates = as.Date(as.character(Sion$time), format="%Y%m%d")
Sierre$dates = as.Date(as.character(Sierre$time), format="%Y%m%d")

# remove NA values
Sierre = Sierre %>% rename(precip=rka150d0) %>% filter(precip != "-") %>% 
  select(dates, precip, stn)
Sierre$precip = as.numeric(Sierre$precip)

# match time period
Sion_comp = Sion %>%  filter(dates %in% Sierre$dates) %>% rename(precip=rka150d0) %>% 
  select(dates, precip, stn)

meteo_comp = rbind(Sion_comp, Sierre)
meteo_comp$year = year(meteo_comp$dates)

# compare yearly sums
meteo_yr = meteo_comp %>% group_by(year, stn) %>% summarize_at(vars(precip), list(sum))
ggplot(meteo_yr, aes(year, precip, fill=stn))+geom_col(position="dodge")

meteo_yr2 = meteo_yr %>% pivot_wider(names_from=stn, values_from=precip)
t.test(meteo_yr2$SIO, meteo_yr2$VSSIE) # yearly means are similar

# compare daily values
meteo_wide = meteo_comp %>% pivot_wider(names_from=stn, values_from=precip)
ggplot(meteo_wide, aes(VSSIE, SIO))+geom_point()
t.test(meteo_wide$SIO, meteo_wide$VSSIE) # daily means are similar


## use Sion precip until 2014, then Sierre

meteoLWF = Sion %>%
  mutate( tmin = tre200dn,
          tmax = tre200jx,
          tmean = tre200d0,
          prec_Sion = rka150d0,
          relhum = ure200d0,
          globrad = ((24*60*60)/1000000)*gre000d0, # convert from W/m2 to MJ/day/m2
          windspeed = fkl010d0,
          vappres = RHtoVPD(relhum, tmean, Pa = 101)) %>%
  select(dates, globrad, tmax, tmin, tmean, vappres, windspeed, prec_Sion) %>% 
  filter(dates<as.Date("2025-01-01"))

meteoLWF = left_join(meteoLWF, select(rename(Sierre, precip_Sierre=precip), -stn))

# fill in missing values in Sierre column with Sion data
meteoLWF$precip_ctrl = meteoLWF$precip_Sierre
meteoLWF$precip_ctrl[is.na(meteoLWF$precip_ctrl)] = meteoLWF$prec_Sion[is.na(meteoLWF$precip_ctrl)]

# now using precip column to combine with irrigation
meteoLWF = left_join(meteoLWF, irr)
meteoLWF$irrig_mm[is.na(meteoLWF$irrig_mm)] = 0
meteoLWF$precip_irr = meteoLWF$precip_ctrl + meteoLWF$irrig_mm

# irrigation stop
meteoLWF$precip_irrstp = if_else(meteoLWF$dates < "2014-01-01", meteoLWF$precip_irr, meteoLWF$precip_ctrl)

#write_csv(meteoLWF, "../../Data/Pfyn/meteo_irr.csv")

## plot irrigation amounts
meteo_plot = meteoLWF %>% filter(dates >= "2003-01-01") %>% select(-c(2:9,11)) %>% pivot_longer(-dates, names_to="treatment", values_to="mm")

# daily
ggplot(meteo_plot, aes(dates, mm, color=treatment)) + geom_line() + facet_wrap(~treatment, ncol=1)

meteo_plot$year = year(meteo_plot$dates)
meteo_plot$month = month(meteo_plot$dates)
meteo_plot2 = meteo_plot %>% group_by(year, month, treatment) %>% summarize_at(vars(mm), list(sum))
meteo_plot2$treatment = factor(meteo_plot2$treatment, levels=c("precip_irr", "precip_irrstp", "precip_ctrl"))
meteo_plot2$date = as.Date(paste(meteo_plot2$year, meteo_plot2$month, 1, sep="-"))

# monthly
ggplot(filter(meteo_plot2, treatment!="precip_irrstp"), aes(x=date, y=mm, fill=treatment))+geom_bar(position="stack", stat="identity")

# yearly comparison
meteo_yr = meteo_plot %>% group_by(year, treatment) %>% summarize_at(vars(mm), list(sum)) %>% 
  group_by(treatment) %>% summarize_at(vars(mm), list(mean))
meteo_yr
