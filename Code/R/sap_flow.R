# analyze sap flow data in Pfynwald

library(tidyverse)
library(lubridate)
library(zoo)

# sap flow data from Zweifel et al. (2020)

sf = read_tsv("../../Data/Pfyn/pfynwald_sapflow.tab", skip=33)

# lengthen data frame and separate column names into multiple variables
sf_long = sf %>% pivot_longer(cols=-1, names_to=c(".value", "tree_id","scenario"), 
                              names_pattern="^(.*?) \\[.*?\\] \\(at tree (\\d+),\\s*([^\\)]+)\\)")

# remove plot from scenario names
sf_long$scenario <- sf_long$scenario %>% str_replace(" plot$", "") %>% str_trim()

# remove spaces from column names
sf_long = sf_long %>% rename(deltaT='delta T', Tree_SR='Tree SR', datetime='Date/Time')

ggplot(sf_long, aes(x=datetime, y=Tree_SR, color=as.factor(tree_id)))+geom_line()+
  facet_wrap(~tree_id)

# summarize daily sap flow
sf_long$date = as.Date(sf_long$datetime)

sf_daily = sf_long %>% select(-datetime) %>% group_by(date, tree_id, scenario) %>% 
  summarize(sfd=mean(Tree_SR, na.rm=T)) %>% mutate(year=year(date))
sf_daily = na.omit(sf_daily)

ggplot(sf_daily, aes(date, sfd, color=as.factor(tree_id)))+geom_line(alpha=.8)+guides(color="none")+
  facet_wrap(~scenario, ncol=1)+theme_bw()+labs(x="",y="Mean Daily Sap Flux [kg/h]", color="Tree ID")

ggplot(filter(sf_daily, scenario!="irrigation"), aes(date, sfd, color=as.factor(scenario)))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()

# group by scenario
sf_meta = sf_daily %>% filter(tree_id != "125") %>% group_by(date, scenario) %>% summarize(sfd=mean(sfd, na.rm=T))
sf_meta$year = year(sf_meta$date)

ggplot(filter(sf_meta, scenario!="irrigation"), aes(date, sfd, color=as.factor(scenario)))+geom_line()+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()+
  theme(legend.position="inside",legend.position.inside=c(.1,.95))+
  labs(x="", y="Mean Daily Sap Flux [kg/h]", color="Scenario")

# write_csv(sf_meta, "Pfyn_sap_2011_17.csv")


# sap flow data from Richard Peters

setwd("../../Data/Pfyn/data for Pfynwald/")

# VPDrought 2024 data

sap_vpd = readRDS("2024/VPDrought_SF_L3_2025-02-02.RDS")

sap_vpd$Date.Time = as.POSIXct(sap_vpd$Date.Time, tz="CET", format="%Y-%m-%d %H:%M:%S")
sap_vpd$date = as.Date(sap_vpd$Date.Time, tz="CET")

ggplot(sap_vpd, aes(Date.Time, Total_Sap_Flow, color=as.factor(Tree_id)))+geom_line()+
  facet_wrap(~Treatment)+theme_bw()+guides(color="none")

sap_vpd_daily = sap_vpd %>% filter(Total_Sap_Flow>0) %>% group_by(date, Tree_id, Treatment) %>% 
  summarize(Total_Sap_Flow=mean(Total_Sap_Flow, na.rm=T))

ggplot(sap_vpd_daily, aes(date, Total_Sap_Flow, color=as.factor(Tree_id)))+geom_line()+
  facet_wrap(~Treatment)+theme_bw()+guides(color="none")

sap_vpd_meta = sap_vpd_daily %>% group_by(date, Treatment) %>% 
  summarize(Total_Sap_Flow=mean(Total_Sap_Flow, na.rm=T))

ggplot(sap_vpd_meta, aes(date, Total_Sap_Flow, color=as.factor(Treatment)))+geom_line()+
  theme_bw()+labs(x="",color="Treatment")

# write_csv(sap_vpd_meta, "../Pfyn_sap_vpd.csv")


# 2021 - 2022 data

sap = readRDS("2021-2022/PFY_sfd_cleaned.Rds")
sap$tree = factor(sap$tree)

# remove outliers
sap$.annotation = ifelse(is.na(sap$.annotation), "OK", sap$.annotation)
sap = filter(sap, .annotation!="Outlier")

sap = select(sap, -c(.dcrkey, .annotation, selection_count, k))

ggplot(sap, aes(timestamp, sfd, color=tree))+geom_line()+facet_wrap(~tree)+
  guides(color="none")+theme_bw()

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
sap_daily_test = na.omit(sap_daily_test)

sap_daily_meta = sap_daily_test %>% group_by(date, meta) %>% summarize(sfd=mean(sfd, na.rm=T))
#write_csv(sap_daily_meta, "../../Data/Pfyn/PFY_sap.csv")

ggplot(sap_daily_meta, aes(date, sfd, color=meta))+geom_point()

sap_daily_meta$year = year(sap_daily_meta$date)

ggplot(sap_daily_meta, aes(date, sfd, color=meta))+geom_line()+
  facet_wrap(~year, scales="free_x", ncol=1)+theme_bw()+
  labs(x="", y="Mean Daily Sap Flux Density", color="Scenario")+
  theme(legend.position="inside", legend.position.inside=c(.1, .1))


# older data

file_list = list.files("older/sap flow/", "*.RDS")

# import individual files into list
sap_list = sapply(file_list, function(x) readRDS(paste0("older/sap flow/", x)))

# get tree IDs from file names
tree_ids = gsub(".*_([0-9]+)\\.RDS$", "\\1", file_list)

# convert from zoo objects to data frames in list
sap_list_df = lapply(sap_list, function(x) data.frame(timestamp=index(x), coredata(x)))

# add tree ID as column
sap_list_df = map2(sap_list_df, tree_ids, ~cbind(.x, TreeID = .y))

# combine into single data frame
sap_df = do.call(rbind, sap_list_df)

sap_df = na.omit(sap_df)

# aggregate to daily
sap_df$date = as.Date(sap_df$timestamp)
sap_daily = sap_df %>% group_by(date, TreeID) %>% summarize_at(vars(sfd, q025, q975), list(mean))

ggplot(sap_daily, aes(date, sfd, color=as.factor(TreeID)))+geom_line()+
  facet_wrap(~TreeID)+theme_bw()+guides(color="none")

# write_csv(sap_daily, "older/sap flow/Pfyn_sap_old.csv")


# TreeNet data

setwd("../TreeNet")

TN_ctr = read_csv("tn_timeseries_Pfyn_sap_control.csv")
TN_irr = read_csv("tn_timeseries_Pfyn_sap_irrigation.csv")

TN_ctr_meta = read_csv("tn_metadata_Pfyn_sap_control.csv")
TN_irr_meta = read_csv("tn_metadata_Pfyn_sap_irrigation.csv")

# process raw data with TREX
library(TREX)
library(shiny)
library(plotly)

# one time series at a time
raw = TN_ctr %>% filter(series_id==932) %>% select(-series_id) %>% rename(timestamp=ts)

#input = is.trex(raw, tz="Etc/GMT-1", time.format="%Y-%m-%d %H:%M:%S", solar.time=F)

# manual check is.trex 
# (circumvent problem with as.character(timestamp) removing midnight time)

time.format="%Y-%m-%d %H:%M:%S"
tz="Etc/GMT-1"
data = raw

data = filter(data, value > 0, value < 20)

# return function output
input = zoo::zoo(data$value, order.by = data$timestamp)

# clip input and aggregate to hourly
input_clip <- dt.steps(input=input, 
                  start="2025-01-01 00:00:00",
                  end="2025-07-13 23:50:00",
                  time.int=60,
                  max.gap=120,
                  df=TRUE)

# save(input_clip, file="input24.RData")

outlier()

input_list = readRDS("input24_Cleaned.Rds")
# input_clean = input_list$series_input
out = rbind(input_list$selected_data_auto, input_list$selected_data_manual)

input_clean = input_clip

# manually remove outliers
input_clean$value[input_clean$timestamp %in% as.character(out$timestamp)] = NA

# fix midnight values
ts = as.POSIXct(input_clean$timestamp, tz="UTC", format=time.format)
ts[is.na(ts)] = ts[which(is.na(ts))+1] - 3600
input_clean$timestamp = ts # format(ts, format=time.format)

input_clean = na.omit(input_clean)
input_fix = zoo(input_clean$value, order.by=input_clean$timestamp)

input_long = rbind(input_long, input_clean)
input_final = rbind(input_final, input_fix)

# calculate maximum difference

output.max <- tdm_dt.max(input_final,
                    methods = "pd", # c("pd", "mw", "dr"),
                    max.days = 10)

plot(output.max$input, ylab = "dC")
lines(output.max$max.pd, col = "green")
# lines(output.max$max.dr, col = "orange")

# calculate sap flux density
output.data<- tdm_cal.sfd(output.max, wood="Coniferous")
sfd_data <- output.data$sfd.pd$sfd

plot(sfd_data)

sfd_df = data.frame(timeseries=index(output.data$sfd.pd), coredata(output.data$sfd.pd))
# write_csv(sfd_df, "Pfyn_sap_tree109.csv")