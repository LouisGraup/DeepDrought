## retrieve and organize Pfynwald soil water data for VPDrought experiment

library(lubridate)
library(tidyverse)
library(rLWFpg)
library(humidity)


start_date <- as.POSIXct("2023-08-01 00:00:00", tz="GMT")
end_date <- as.POSIXct("2025-07-01 00:00:00", tz="GMT")

conn <- db_connect('grauplou', rstudioapi::askForPassword("Database password"), db_host = 'pgdblwf', db_name = 'lwf')


messvar.df = db_tbl(conn, 'ada', 'v_messvar', retrieve = T) %>%
  filter(project_name %in% "Pfynwald irrigation")

soilvar.tbl = db_tbl(conn, 'ada', 'v_messvar', retrieve = F) %>%
  filter(project_name == "Pfynwald irrigation", vartable_name == "pfyn_messdat",
         varname_name %in% c("Soil water content", "Soil water matric potential", "Soil temperature"),
         varfreq_name == "10 Minutes", grepl("PFY_SC", installation_name))
soilvar.df = collect(soilvar.tbl)

pfyndata.tbl = db_tbl(conn, table="pfyn_messdat", retrieve = F)

# summarize soil water measurements to daily time step
pfyn_daily.tbl = pfyndata.tbl %>% 
  filter(messvar_id %in% soilvar.df$messvar_id,
         messtime >= start_date, messtime < end_date) %>% 
  mutate(dates = as.Date(messtime)) %>% group_by(dates, messvar_id) %>% 
  summarize(messval = mean(messval))

pfyn_daily.df = soilvar.tbl %>% inner_join(pfyn_daily.tbl, by="messvar_id") %>% 
  select(dates, messval, messvar_id, messvar_name, varname_name, varvpos, type, scaffold, treatment, descr) %>% collect()

# triple replicates for most measurements
SWC_df = pfyn_daily.df %>% filter(varname_name == "Soil water content")
SWC_df$varvpos = factor(SWC_df$varvpos, levels = c(-0.1, -0.8, -1.2))

SMP_df = pfyn_daily.df %>% filter(varname_name == "Soil water matric potential")
SMP_df$varvpos = factor(SMP_df$varvpos, levels = c(-0.01, -0.10, -0.80, -1.20))

# plot mean and standard deviation per treatment and depth
ggplot(SWC_df, aes(dates, messval, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~varvpos, ncol=1)+theme_bw()+
  labs(x="", y="SWC")

# plot mean and standard deviation per treatment and depth
ggplot(filter(SMP_df, type!="MPS 6", varvpos != -0.01), aes(dates, messval/1000, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~varvpos, scales="free_y", ncol=1)+theme_bw()+
  labs(x="", y="SMP (MPa)")

# plot mean and standard deviation per depth for roof treatments
ggplot(filter(SMP_df, treatment %in% c("roof_vpd", "roof"), type!="MPS 6", varvpos != -0.01), aes(dates, messval/1000, color=varvpos, fill=varvpos))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~treatment)+theme_bw()+labs(x="", y="SMP (MPa)")

# plot roof treatments 
ggplot(filter(SMP_df, type!="MPS 6", !is.na(descr)), aes(dates, messval/1000, group=as.factor(messvar_id)))+geom_line(aes(linetype=as.factor(varvpos), color=as.factor(scaffold)))+
  facet_grid(descr~treatment)+theme_bw()+labs(x="", y="SMP (MPa)")

# plot roof treatments from recent sensors only
ggplot(filter(SMP_df, messvar_id>10000), aes(dates, messval/1000, group=as.factor(messvar_id)))+geom_line(aes(color=as.factor(scaffold)))+
  facet_grid(descr~treatment)+theme_bw()+labs(x="", y="SMP (MPa)")


# examine sub-daily SMP measurements for temp. correction

SMP_test = pfyndata.tbl %>% filter(messvar_id == 5672, messtime >= start_date, messtime < end_date) %>% collect()

ggplot(filter(SMP_test, messtime >= as.POSIXct("2024-05-01 00:00:00", tz="GMT"), messtime < as.POSIXct("2024-05-04 00:00:00", tz="GMT")), aes(messtime, messval))+geom_line()