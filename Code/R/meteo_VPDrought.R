## retrieve and organize Pfynwald meteo data for VPDrought experiment
## code borrowed from Stefan Hunziker

library(lubridate)
library(tidyverse)
library(rLWFpg)
library(humidity)


start_date <- as.POSIXct("2024-01-01 00:00:00", tz="GMT")
end_date <- as.POSIXct("2025-01-01 00:00:00", tz="GMT")

conn <- db_connect('grauplou', rstudioapi::askForPassword("Database password"), db_host = 'pgdblwf', db_name = 'lwf')

metadata.tbl <- db_tbl(conn, 'ada', 'v_messvar', retrieve = F) %>%
  filter(project_name %in% "Pfynwald irrigation",
         (grepl("PFY_SC", installation_name)),
         varname_name %in% c("Air temperature", "Atmospheric vapour pressure deficit"),
         varfreq_name %in% '1 Minutes')
metadata.df = collect(metadata.tbl)

# only consider center in-canopy measurements
meta_filter.tbl = metadata.tbl %>% filter(grepl("ic_ce", messvar_name))

pfyndata.tbl = db_tbl(conn, table="pfyn_messdat", retrieve = F)

# combine data with metadata columns and filter dates
# and average over replicates for each treatment
meteo.tbl = meta_filter.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  select(messtime, messval, messvar_name, treatment) %>% 
  filter(messtime >= start_date, messtime < end_date) %>% 
  group_by(messtime, messvar_name, treatment) %>% summarize(messval=mean(messval))

# reshape into vpd and air temp columns
meteo.wide = meteo.tbl %>% pivot_wider(values_from=messval, names_from=messvar_name) %>% 
  mutate(dates = as.Date(messtime))

# summarize daily values for VPD and temperature
meteo.daily = meteo.wide %>% group_by(dates, treatment) %>% summarize(VPD=mean(VPD_ic_ce), tmax=max(TA_ic_ce), tmin=min(TA_ic_ce), tmean=mean(TA_ic_ce))

meteo.df = collect(meteo.daily)

#write_csv(meteo.df, "Pfyn_micro_clim_2024.csv")

# drop roof scenario and aggregate into vpd and control scenarios
meteo.df = filter(meteo.df, treatment != "roof")
meteo.df$scen = ifelse(grepl("vpd", meteo.df$treatment),"vpd","con")

meteo.scen = meteo.df %>% select(-treatment) %>% group_by(dates, scen) %>% summarize_all(list(mean))

# actual relative reduction of VPD between control and vpd scenarios
vpd_diff = (meteo.scen$VPD[meteo.scen$scen=="vpd"] - meteo.scen$VPD[meteo.scen$scen=="con"]) / meteo.scen$VPD[meteo.scen$scen=="con"]
mean(vpd_diff)


## need to convert VPD to actual vapor pressure

# first calculate SVP using temperature and convert to kPa
meteo.scen$Es = SVP.ClaCla(meteo.scen$tmean+273.15) / 10

# then actual vapor pressure is just SVP - VPD
meteo.scen$Ea = meteo.scen$Es - meteo.scen$VPD

ggplot(meteo.scen, aes(dates, Ea, color=scen))+geom_line()+geom_line(aes(dates, Es, color=scen), linetype="dashed")

# separate scenarios
meteo.ctr = meteo.scen %>% filter(scen == "con") %>% select(-c(scen, VPD, Es)) %>% 
  rename(vappres = Ea)

meteo.vpd = meteo.scen %>% filter(scen == "vpd") %>% select(-c(scen, VPD, Es)) %>% 
  rename(vappres = Ea)


## retrieve remaining in-situ meteo data
# precipitation, radiation, wind speed

messvar.tbl <- db_tbl(conn, schema = "ada", table = "v_messvar", retrieve = F) %>%
  dplyr::filter(project_name %in% 'Pfynwald irrigation', messvar_id %in% c(124, 131, 134))

meteovar.df = messvar.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  select(messtime, messval, messvar_name) %>% 
  filter(messtime >= start_date, messtime < end_date) %>% 
  pivot_wider(names_from=messvar_name, values_from=messval) %>% 
  mutate(dates = as.Date(messtime)) %>% collect()

meteovar.df$WS_2m = ifelse(meteovar.df$WS_2m>6, 1, meteovar.df$WS_2m)

meteovar.daily = meteovar.df %>% group_by(dates) %>% 
  summarize(precip = sum(RaineH3_amount_dif_Tot, na.rm=T),
            globrad = mean(GlobalRad_oc_Avg, na.rm=T)*.0864,
            windspeed=mean(WS_2m, na.rm=T))


# combine meteo data
meteo_Pfyn_VPD = left_join(meteovar.daily, meteo.vpd)
meteo_Pfyn_Con = left_join(meteovar.daily, meteo.ctr)

#write_csv(meteo_Pfyn_VPD, "Pfyn_meteo2024_insitu_VPD.csv")
#write_csv(meteo_Pfyn_Con, "Pfyn_meteo2024_insitu_Control.csv")

dbDisconnect(conn)