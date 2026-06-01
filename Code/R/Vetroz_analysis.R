# combine and analyze Vetroz data
# data processing below

library(tidyverse)
library(lubridate)
library(gridExtra)

setwd("../../Data/Daten_Vetroz_Lorenz")

## leaf water potential
SWP_daily = read_csv("SWP_daily.csv")

LWP = read_csv("LWP_gs.csv")
LWP$Date = as.Date(LWP$Date, format="%m/%d/%y")
LWP$year = year(LWP$Date)
LWP$pd_md = if_else(LWP$Time < hms("06:00:00"), "pd", "md")
LWP$LWP_mean = LWP$LWP_mean / -10 # convert to MPa

# compare bagged and unbagged
ggplot(filter(LWP, Date < "2021-10-01", pd_md=="md"), aes(Date, LWP_mean, fill=Method, group=interaction(Date, Method)))+geom_boxplot()+
  facet_wrap(~TreeNr)+theme_bw()

LWP_BWP = LWP %>% filter(pd_md=="md") %>% select(Date, TreeNr, LWP_mean, Method) %>% 
  group_by(Date, TreeNr, Method) %>% summarize(LWP=mean(LWP_mean)) %>% 
  pivot_wider(names_from=Method, values_from=LWP) %>% 
  rename(LWP = unbagged, BWP = bagged) %>% mutate(WP_diff = BWP - LWP)

ggplot(LWP_BWP, aes(LWP, WP_diff))+geom_point()+
  stat_smooth()+theme_bw()+
  labs(x="Unbagged Leaf Water Potential (MPa)", y="Difference between bagged and unbagged LWP (MPa)")


# compare against soil water potential
ggplot(filter(LWP, TreeNr %in% c(3, 8)), aes(SMP_mean / 1000, LWP_mean, color=pd_md))+geom_point()+
  geom_abline(slope=1, intercept=0)+facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Soil Water Potential (MPa)", y = "Leaf Water Potential (MPa)")

# sap flow and leaf water potential
ggplot(filter(LWP, TreeNr %in% c(3, 8)), aes(LWP_mean, sapflow))+geom_point()+
  facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Leaf Water Potential (MPa)", y = "Sapflow")

# TWD and LWP
ggplot(LWP, aes(TWD, LWP_mean, color=pd_md))+geom_point()+
  facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Tree Water Deficit (MPa)", y = "Leaf Water Potential (MPa)")

# gs and LWP
ggplot(LWP, aes(gs, LWP_mean))+geom_point()+
  facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Stomatal Conductance (mol m^-2 s^-1)", y = "Leaf Water Potential (MPa)")


## process soil water potential data
SWP = read_csv("SWP_20_80_110_160cm_hourly_2014_2022.csv")
SWP1cm = read_csv("SWP_1cm_Vetroz.csv")

# need to add time to midnight
SWP$MEAS_DATE = ifelse(str_length(SWP$MEAS_DATE) == 10,
                   paste(SWP$MEAS_DATE, "00:00:01"), SWP$MEAS_DATE)
# convert to date-time and extract depths
SWP$datetime = as.POSIXct(SWP$MEAS_DATE, format="%d.%m.%Y %H:%M:%S", tz="UTC")
SWP$date = as.Date(SWP$datetime)
SWP$depth = parse_number(str_extract(SWP$NAME, "\\d+\\s*cm"))
SWP = rename(SWP, value = `VALUE kPa`)
SWP$value[SWP$value == -9999999] = NA

# group to daily
SWP_daily = SWP %>% group_by(date, depth) %>% summarize(SWP=mean(value))

ggplot(SWP_daily, aes(date, SWP, color=as.factor(depth)))+geom_line()+
  labs(x="", y="SWP (kPa)", color="Depth (cm)")+theme_bw()

# add 1cm data
SWP1cm$datetime = as.POSIXct(SWP1cm$Date, format="%d.%m.%Y %H:%M", tz="UTC")
SWP1cm_daily = SWP1cm %>% mutate(date=as.Date(datetime)) %>% 
  rename(value = SWP_dry_corr_kPa) %>% group_by(date) %>% 
  summarize(SWP=mean(value)) %>% mutate(depth=1)

ggplot(SWP1cm_daily, aes(date, log10(SWP*-10)))+geom_line()+labs(x="", y="pF")+theme_bw()

SWP1cm_daily$SWP[SWP1cm_daily$SWP < -4000] = NA # cap at -4000 kPa

SWP_daily = rbind(SWP_daily, SWP1cm_daily)

SWP_daily$year = year(SWP_daily$date)
SWP_daily$month = month(SWP_daily$date)

# remove NAs for output
SWP_daily = na.omit(SWP_daily)
#write_csv(SWP_daily, "SWP_daily.csv")

## sap flow processing
sap_files = list.files(path = ".", pattern="Sapflow_node")
sap_ids = str_extract(sap_files, "\\d+\\s*") # extract node ids from sap_files
sap = lapply(sap_files, read_csv)

# filter only relevant columns
sap_list = lapply(sap, function(df) {
  df %>% select(Date, Time, `Sap Flow Out (kg/hr)`) %>% rename(sapflow = `Sap Flow Out (kg/hr)`) })

# apply node ids to list elements
names(sap_list) = sap_ids

# combine sap flow list into single data frame with individual IDs
sap_df = bind_rows(sap_list, .id="node_id")

# combine date and time columns
sap_df$datetime = as.POSIXct(paste(sap_df$Date, sap_df$Time), format="%m/%d/%Y %H:%M:%S")

# define root and stem sensors
sap_df = sap_df %>% mutate(sensor_loc = ifelse(node_id %in% c("3", "8"), "stem", "root"))

# group to daily
sap_daily = sap_df %>% mutate(datehour=floor_date(datetime, "1 hour")) %>%
  group_by(datehour, node_id, sensor_loc) %>% summarize(sapflow=mean(sapflow, na.rm=T), .groups="drop") %>%
  mutate(date=as.Date(datehour)) %>% select(-datehour) %>%
  group_by(date, node_id, sensor_loc) %>% summarize_at(vars(sapflow), list(sum))

ggplot(sap_daily, aes(date, sapflow, color=sensor_loc))+geom_line()+facet_wrap(~node_id)

# fix zero offset for some sensors
sap_daily$sapflow = ifelse(sap_daily$node_id %in% c("3", "7", "9"), sap_daily$sapflow - 1,
                           ifelse(sap_daily$node_id == "6", sap_daily$sapflow - 0.5, sap_daily$sapflow))

ggplot(sap_daily, aes(date, sapflow, color=sensor_loc))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", aes(fill=sensor_loc), alpha=0.3)+
  labs(x="", y="Daily Sap Flow (kg/day)", color="Location", fill="Location")+theme_bw()

# mean sap flow
sap_mean = sap_daily %>% group_by(date, sensor_loc) %>% summarize(sapflow=mean(sapflow, na.rm=T)) %>%
  filter(date >= as.Date("2021-05-01"), date < as.Date("2023-12-19"))

ggplot(sap_mean, aes(date, sapflow, color=sensor_loc))+geom_line()+
  labs(x="", y="Daily Sap Flow (kg/day)", color="Location")+theme_bw()

sap_mean$month = month(sap_mean$date)
sap_mean$year = year(sap_mean$date)
#write_csv(sap_mean, "sapflow_daily.csv")

sap_daily = read_csv("sapflow_daily.csv")
ggplot(sap_daily, aes(date, sapflow, color=sensor_loc))+geom_line()+
  labs(x="", y="Daily Sap Flow (kg/day)", color="Location")+theme_bw()

## process dendrometer data

# function to normalize pre-dawn and midday TWD with MDS
NORM_DENDRO = function(den_long) {
  # expect den_long with datetime, tree_name, tree_id, sensor_loc, and TWD columns
  
  # formatting
  den_long = na.omit(den_long)
  den_long$date = as.Date(den_long$datetime)
  den_long$hour = hour(den_long$datetime)
  
  # aggregate to daily with pre-dawn
  den_pd = den_long %>% filter(hour <= 6) %>% 
    group_by(date, tree_name, tree_id, sensor_loc) %>% summarize(TWD_pd=min(TWD, na.rm=T))
  
  # aggregate to daily with mid-day
  den_md = den_long %>% filter(hour >= 12, hour <= 17) %>% 
    group_by(date, tree_name, tree_id, sensor_loc) %>% summarize(TWD_md=max(TWD, na.rm=T))
  
  den_daily = full_join(den_pd, den_md)
  
  # calculate maximum daily shrinkage
  den_daily$MDS = den_daily$TWD_md - den_daily$TWD_pd
  
  # calculate 99th percentile MDS over entire time series
  den_mds = den_daily %>% mutate(month=month(date)) %>% 
    filter(month>4, month<11, MDS>0) %>%  # only accept positive shrinkage during summer months
    group_by(tree_name) %>% summarize(mds_max = quantile(MDS, .99))
  
  # join max MDS to daily data frame
  den_daily = left_join(den_daily, den_mds)
  
  # normalize pre-dawn TWD and MDS with max MDS
  den_daily$TWD_pdn = den_daily$TWD_pd / den_daily$mds_max
  den_daily$MDS_norm = den_daily$MDS / den_daily$mds_max
  
  return(den_daily)
}


# TreeNet data
TN = read_csv("TreeNet_dendro/tn_timeseries.csv")
TN_meta = read_csv("TreeNet_dendro/tn_metadata.csv")

TN_meta$sensor_loc = if_else(grepl("Root", TN_meta$tree_name), "root", "stem")
TN_meta$tree_id = str_extract(TN_meta$tree_name, "\\d+\\s*")

TN = TN %>% filter(frost == FALSE) %>%
  select(-c(frost, flags, gro_start, gro_end)) %>% 
  left_join(select(TN_meta, series_id, tree_name, tree_id, sensor_loc)) %>% 
  rename(datetime=ts, TWD=twd) %>% select(-series_id)

ggplot(TN, aes(datetime, value, color=tree_id, group=tree_name))+geom_line()+
  facet_wrap(~sensor_loc, scales="free")+guides(color="none")+theme_bw()+
  labs(x="", y="Radius")

TWD_norm = NORM_DENDRO(TN)

# plots
ggplot(TWD_norm, aes(date, TWD_pd, color=tree_id, group=tree_name))+geom_line()+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(date, TWD_md, color=tree_id, group=tree_name))+geom_line()+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(TWD_pd, MDS, color=tree_name))+geom_point(size=.5, alpha=.5)+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(date, TWD_pdn, color=tree_id, group=tree_name))+geom_line()+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(date, TWD_pdn, group=tree_name))+geom_line(color="green")+
  geom_line(aes(date, MDS_norm, group=tree_name), color="black")+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()

ggplot(TWD_norm)+stat_summary(geom="line", fun=mean, aes(date, TWD_pdn, color="TWD_pdn"))+
  stat_summary(geom="ribbon", aes(date, TWD_pdn, color="TWD_pdn", fill="TWD_pdn"), alpha=0.2)+
  stat_summary(geom="line", fun=mean, aes(date, MDS_norm, color="MDS_norm"))+
  stat_summary(geom="ribbon", aes(date, MDS_norm, color="MDS_norm", fill="MDS_norm"), alpha=0.2)+
  scale_colour_manual(name="Var", values=c("MDS_norm"="black","TWD_pdn"="green"))+
  scale_fill_manual(name="Var", values=c("MDS_norm"="grey","TWD_pdn"="green"))+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()


# mean daily TWD

TWD_daily = TWD_norm %>% select(-c(tree_name, tree_id, mds_max)) %>% 
  group_by(date, sensor_loc) %>% summarize_all(list(mean))

ggplot(TWD_daily, aes(date, TWD_pdn))+geom_line(color="green")+
  geom_line(aes(date, MDS_norm), color="black")+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()

# exclude winter months
TWD_daily$year = year(TWD_daily$date)
TWD_daily$month = month(TWD_daily$date)
TWD_daily[TWD_daily$month < 5 | TWD_daily$month > 11, c("TWD_pdn", "MDS_norm")] = NA


# stacked plots
p_sap = sap_daily %>% mutate(ddate=as.Date(paste(2000, month, day(date), sep="-"))) %>% 
  filter(date < "2022-11-01", month > 4, month < 11) %>% 
  ggplot()+geom_line(aes(ddate, sapflow, color=sensor_loc))+
  labs(x="", y="Daily Sap Flow (kg/day)", color="Location")+theme_bw()+
  scale_x_date(date_breaks="1 month", date_labels="%b")+facet_wrap(~year)+
  theme(legend.position="inside", legend.position.inside=c(0.45,0.75))
p_twd = TWD_daily %>% mutate(ddate=as.Date(paste(2000, month, day(date), sep="-"))) %>%
  filter(date >= "2021-05-07", date < "2022-11-01", month>4, month<11) %>% 
  ggplot()+geom_line(aes(ddate, TWD_pdn, color=sensor_loc, linetype="TWD_pdn"))+
  geom_line(aes(ddate, MDS_norm, color=sensor_loc, linetype="MDS_norm"))+
  geom_hline(yintercept=0, color="darkgrey")+
  theme_bw()+labs(x="", y="Norm. TWD / MDS")+guides(color="none")+
  scale_x_date(date_breaks="1 month", date_labels="%b")+facet_wrap(~year)+
  scale_linetype_manual(name="Var", values=c("MDS_norm"=3,"TWD_pdn"=1))+
  theme(legend.position="inside", legend.position.inside=c(0.1,0.75))
p_wp = SWP_daily %>% mutate(ddate=as.Date(paste(2000, month, day(date), sep="-"))) %>%
  filter(date >= "2021-05-07", date < "2022-11-01", month>4, month<11) %>% 
  ggplot()+geom_line(aes(ddate, SWP/1000, color=as.factor(depth)))+
  geom_point(data=filter(mutate(LWP, ddate=as.Date(paste(2000, month(Date), day(Date), sep="-"))), pd_md=="pd"), aes(ddate, LWP_mean), color="red", inherit.aes=F)+
  labs(x="", y="SWP / pre-dawn LWP (MPa)", color="Depth (cm)")+theme_bw()+
  scale_x_date(date_breaks="1 month", date_labels="%b")+facet_wrap(~year)+
  theme(legend.position="inside", legend.position.inside=c(0.05,0.4))

grid.arrange(p_sap, p_twd, p_wp)

