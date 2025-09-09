# analyze sap flow data in Pfynwald

library(tidyverse)
library(lubridate)

# sap flow data from Zweifel et al. (2020)

sf = read_tsv("../../Data/Pfyn/pfynwald_sapflow.tab", skip=33)

# lengthen data frame and separate column names into multiple variables
sf_long = sf %>% pivot_longer(cols=-1, names_to=c(".value", "tree_id","scenario"), 
                              names_pattern="^(.*?) \\[.*?\\] \\(at tree (\\d+),\\s*([^\\)]+)\\)")

# remove plot from scenario names
sf_long$scenario <- sf_long$scenario %>%
  str_replace(" plot$", "") %>%
  str_trim()

# remove spaces from column names
sf_long = sf_long %>% rename(deltaT='delta T', Tree_SR='Tree SR', datetime='Date/Time')

ggplot(sf_long, aes(x=datetime, y=Tree_SR, color=as.factor(tree_id)))+geom_line()+
  facet_wrap(~tree_id)

# summarize daily sap flow
sf_long$date = as.Date(sf_long$datetime)

sf_daily = sf_long %>% select(-datetime) %>% group_by(date, tree_id, scenario) %>% 
  summarize(sfd=mean(Tree_SR, na.rm=T)) %>% mutate(year=year(date))

ggplot(sf_daily, aes(date, sfd, color=as.factor(tree_id)))+geom_line(alpha=.8)+guides(color="none")+
  facet_wrap(~scenario, ncol=1)+theme_bw()+labs(x="",y="Mean Daily Sap Flux [kg/h]", color="Tree ID")

ggplot(filter(sf_daily, scenario!="irrigation"), aes(date, sfd, color=as.factor(scenario)))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()

# group by scenario
sf_meta = sf_daily %>% group_by(date, scenario) %>% summarize(sfd=mean(sfd, na.rm=T))
sf_meta$year = year(sf_meta$date)

ggplot(filter(sf_meta, scenario!="irrigation"), aes(date, sfd, color=as.factor(scenario)))+geom_line()+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()+
  theme(legend.position="inside",legend.position.inside=c(.1,.95))+
  labs(x="", y="Mean Daily Sap Flux [kg/h]", color="Scenario")


# sap flow data from Richard Peters
sap = readRDS("../../Data/Pfyn/PFY_sfd_cleaned.Rds")
sap$tree = factor(sap$tree)

# remove outliers
sap$.annotation = ifelse(is.na(sap$.annotation), "OK", sap$.annotation)
sap = filter(sap, .annotation!="Outlier")

sap = select(sap, -c(.dcrkey, .annotation, selection_count, k))

ggplot(sap, aes(timestamp, sfd, color=tree))+geom_line()+facet_wrap(~tree)

# group to daily
sap$date = as.Date(sap$timestamp)
sap_daily = sap %>% select(-timestamp) %>% group_by(date, tree) %>% summarize_all(list(mean))

# group by scenario
sap_daily$meta = case_when(grepl("Control", sap_daily$tree) ~ "Control",
                           grepl("Irrigation", sap_daily$tree) ~ "Irrigation",
                           grepl("Stop", sap_daily$tree) ~ "Irrigation Stop")

ggplot(sap_daily, aes(date, sfd, color=tree))+geom_line()+
  geom_ribbon(aes(ymin=q025, ymax=q975, fill=tree), alpha=0.2)+
  facet_wrap(~tree)+theme_bw()+labs(x="", y="Mean Daily Sap Flux Density")+
  theme(legend.position="none", strip.text=element_text(face="bold", size=12))

ggplot(sap_daily, aes(date, sfd, color=tree))+geom_line()+
  facet_wrap(~meta)+theme_bw()+theme(legend.position="none")+
  labs(x="", y="Mean Daily Sap Flux Density")
  
# filter out outliers
sap_daily_test = sap_daily %>% filter(!(tree=="Pfynwald_Irr_Stop_02_16_ch6" & date > "2021-05-15" & date < "2021-06-01"))

sap_daily_meta = sap_daily_test %>% group_by(date, meta) %>% summarize(sfd=mean(sfd, na.rm=T))
#write_csv(sap_daily_meta, "../../Data/Pfyn/PFY_sap.csv")

ggplot(sap_daily_meta, aes(date, sfd, color=meta))+geom_point()

sap_daily_meta$year = year(sap_daily_meta$date)

ggplot(sap_daily_meta, aes(date, sfd, color=meta))+geom_line()+
  facet_wrap(~year, scales="free_x", ncol=1)+theme_bw()+
  labs(x="", y="Mean Daily Sap Flux Density", color="Scenario")+
  theme(legend.position="inside", legend.position.inside=c(.1, .1))


# TreeNet data ... useless
TN_ctr = read_csv("../../Data/Pfyn/TreeNet/tn_timeseries_Pfyn_sap_control.csv")
TN_irr = read_csv("../../Data/Pfyn/TreeNet/tn_timeseries_Pfyn_sap_irrigation.csv")

TN_ctr_meta = read_csv("../../Data/Pfyn/TreeNet/tn_metadata_Pfyn_sap_control.csv")
TN_irr_meta = read_csv("../../Data/Pfyn/TreeNet/tn_metadata_Pfyn_sap_irrigation.csv")

TN_ctr$date = as.Date(TN_ctr$ts)
TN_ctr_daily = TN_ctr %>% group_by(date, series_id) %>% summarize(value=mean(value))

ggplot(TN_ctr_daily, aes(date, value, color=as.factor(series_id)))+geom_line()+
  facet_wrap(~series_id, scales="free")+theme_bw()+theme(legend.position="none")+
  labs(x="", y="")


TN_irr$date = as.Date(TN_irr$ts)
TN_irr_daily = TN_irr %>% group_by(date, series_id) %>% summarize(value=mean(value))

ggplot(TN_irr_daily, aes(date, value, color=as.factor(series_id)))+geom_line()+
  facet_wrap(~series_id, scales="free")+theme_bw()+theme(legend.position="none")+
  labs(x="", y="")

# first attempt at processing with TREX
library(TREX)

# one time series at a time
raw = TN_ctr %>% filter(series_id==932) %>% select(-series_id) %>% rename(timestamp=ts)

#input = is.trex(raw, tz="Etc/GMT-1", time.format="%Y-%m-%d %H:%M:%S", solar.time=F)

# manual check is.trex 
# (circumvent problem with as.character(timestamp) removing midnight time)

time.format="%Y-%m-%d %H:%M:%S"
tz="Etc/GMT-1"
data = raw

# convert time to UTC
timestamp <-
  base::as.POSIXct(
    as.character(data$timestamp),
    format = time.format,
    tz = tz
  )
timestamp<-format(timestamp,tz="UTC",usetz=F)

# fix NA values from midnight (now 23h)
timestamp[is.na(timestamp)] = as.character(as.POSIXct(timestamp[which(is.na(timestamp))-1]) + 600)

# return function output
input = zoo::zoo(data$value, order.by = base::as.POSIXct(paste0(left(timestamp,16),":00"),tz="UTC"))



