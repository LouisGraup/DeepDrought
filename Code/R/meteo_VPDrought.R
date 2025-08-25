## retrieve and organize Pfynwald meteo data for VPDrought experiment
## code borrowed from Stefan Hunziker

library(lubridate)
library(tidyverse)
library(rLWFpg)
library(humidity)


start_date <- as.POSIXct("2024-01-01 00:00:00", tz="GMT")
end_date <- as.POSIXct("2025-08-01 00:00:00", tz="GMT")

conn <- db_connect('grauplou', rstudioapi::askForPassword("Database password"), db_host = 'pgdblwf', db_name = 'lwf')

metadata.tbl <- db_tbl(conn, 'ada', 'v_messvar', retrieve = F) %>%
  filter(project_name %in% "Pfynwald irrigation",
         (grepl("PFY_SC", installation_name)),
         varname_name %in% c("Air temperature", "Atmospheric vapour pressure deficit"),
         varfreq_name %in% '10 Minutes')
metadata.df = collect(metadata.tbl)

# only consider center in-canopy measurements for control scenario
meta_filter_ctr.tbl = metadata.tbl %>% 
  filter(grepl("ic_ce", messvar_name), treatment == "control")

# now consider above-canopy front left measurements for strongest VPD reduction
meta_filter_vpd.tbl = metadata.tbl %>% filter(grepl("oc_fl", messvar_name))

pfyndata.tbl = db_tbl(conn, table="pfyn_messdat", retrieve = F)

# combine data with metadata columns and filter dates
# and average over replicates for each treatment
meteo_ctr.df = meta_filter_ctr.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  select(messtime, messval, messvar_name, treatment) %>% 
  filter(messtime >= start_date, messtime < end_date) %>% 
  group_by(messtime, messvar_name, treatment) %>% 
  summarize(messval=mean(messval), .groups = "drop") %>% collect()

meteo_ctr.df = meteo_ctr.df[order(meteo_ctr.df$messtime),] # sort by messtime

# VPD measurements only have one instance
meteo_vpd.df = meta_filter_vpd.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  select(messtime, messval, messvar_name, treatment) %>% 
  filter(messtime >= start_date, messtime < end_date) %>% 
  collect()

# combine into single data frame with consistent messvar_names
meteo_ctr.df$messvar_name = if_else(grepl("TA",meteo_ctr.df$messvar_name),"TA","VPD")
meteo_vpd.df$messvar_name = if_else(grepl("TA",meteo_vpd.df$messvar_name),"TA","VPD")

meteo.df = rbind(meteo_ctr.df, meteo_vpd.df)

# reshape into vpd and air temp columns
meteo.wide = meteo.df %>% pivot_wider(values_from=messval, names_from=messvar_name) %>% 
  mutate(dates = as.Date(messtime))


## need to convert VPD to actual vapor pressure

# first calculate SVP using temperature and convert to kPa
meteo.wide$Es = SVP.ClaCla(meteo.wide$TA+273.15) / 10

# then actual vapor pressure is just SVP - VPD
meteo.wide$Ea = meteo.wide$Es - meteo.wide$VPD

# summarize daily values for VPD and temperature

meteo.daily = meteo.wide %>% group_by(dates, treatment) %>% 
  summarize(VPD=mean(VPD), tmax=max(TA), tmin=min(TA), tmean=mean(TA),
            Es=mean(Es), Ea=mean(Ea))

#write_csv(meteo.daily, "../../Data/Pfyn/meteo/Pfyn_micro_clim_0124_0725.csv")

# drop roof scenario and aggregate into vpd and control scenarios
#meteo.daily = filter(meteo.daily, treatment != "roof")
meteo.daily$scen = ifelse(grepl("vpd", meteo.daily$treatment),"vpd","con")

#meteo.scen = meteo.daily %>% select(-treatment) %>% group_by(dates, scen) %>% summarize_all(list(mean))
meteo.scen = meteo.daily %>% select(-treatment)

# actual relative reduction of VPD between control and vpd scenarios
vpd_diff = (meteo.scen$VPD[meteo.scen$scen=="vpd"] - meteo.scen$VPD[meteo.scen$scen=="con"]) / meteo.scen$VPD[meteo.scen$scen=="con"]
mean(vpd_diff)


ggplot(meteo.scen, aes(dates, Ea, color=scen))+geom_line()+geom_line(aes(dates, Es, color=scen), linetype="dashed")

# separate scenarios
meteo.ctr = meteo.scen %>% filter(scen == "con") %>% select(-c(scen, VPD, Es)) %>% 
  rename(vappres = Ea)

meteo.vpd = meteo.scen %>% filter(scen == "vpd") %>% select(-c(scen, VPD, Es)) %>% 
  rename(vappres = Ea)

# fill in missing dates
date_seq = seq.Date(from=as.Date("2024-01-01"), to=as.Date("2025-07-31"), by=1)

meteo.vpd = left_join(data.frame(dates=date_seq), meteo.vpd)
meteo.vpd[549:550, 2:5] = meteo.ctr[549:550, 2:5]

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

#write_csv(meteo_Pfyn_VPD, "../../Data/Pfyn/meteo/Pfyn_meteo0124_0725_insitu_VPD.csv")
#write_csv(meteo_Pfyn_Con, "../../Data/Pfyn/meteo/Pfyn_meteo0124_0725_insitu_Control.csv")

dbDisconnect(conn)