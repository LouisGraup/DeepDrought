# combine and analyze Vetroz data
# data processing below

library(tidyverse)
library(lubridate)

setwd("../../Data/Daten_Vetroz_Lorenz")

## leaf water potential
SWP_daily = read_csv("SWP_daily.csv")

LWP = read_csv("LWP_gs.csv")
LWP$Date = as.Date(LWP$Date, format="%m/%d/%y")
LWP$pd_md = if_else(LWP$Time < hms("06:00:00"), "pd", "md")
LWP$LWP_mean = LWP$LWP_mean / -10 # convert to MPa

# compare bagged and unbagged
ggplot(filter(LWP, Date < "2021-10-01", pd_md=="md"), aes(Date, LWP_mean, color=Method, group=interaction(Date, Method)))+geom_boxplot()+
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

# remove NAs for output
SWP_daily = na.omit(SWP_daily)
#write_csv(SWP_daily, "SWP_daily.csv")

## sap flow processing
# sap_files = list.files(path = ".", pattern="Sapflow_node")
# sap = lapply(sap_files, read_csv)
# 
# # filter only relevant columns
# sap_list = lapply(sap, function(df) {
#   df %>% select(Date, Time, `Sap Flow Out (kg/hr)`) %>% rename(sapflow = `Sap Flow Out (kg/hr)`) })
# 
# # combine sap flow list into single data frame with individual IDs
# sap_df = bind_rows(sap_list, .id="node_id")
# 
# # combine date and time columns
# sap_df$datetime = as.POSIXct(paste(sap_df$Date, sap_df$Time), format="%m/%d/%Y %H:%M:%S")
# 
# # group to daily
# sap_daily = sap_df %>% mutate(datehour=floor_date(datetime, "1 hour")) %>%
#   group_by(datehour, node_id) %>% summarize(sapflow=mean(sapflow, na.rm=T), .groups="drop") %>%
#   mutate(date=as.Date(datehour)) %>% select(-datehour) %>%
#   group_by(date, node_id) %>% summarize_at(vars(sapflow), list(sum))
# 
# ggplot(sap_daily, aes(date, sapflow, color=node_id))+geom_line()+facet_wrap(~node_id)
# 
# # fix zero offset for some sensors
# sap_daily$sapflow = ifelse(sap_daily$node_id %in% c("4", "8", "10"), sap_daily$sapflow - 1, 
#                            ifelse(sap_daily$node_id == "7", sap_daily$sapflow - 0.5, sap_daily$sapflow))
# 
# # mean sap flow
# sap_mean = sap_daily %>% group_by(date) %>% summarize(sapflow=mean(sapflow, na.rm=T)) %>%
#   filter(date >= as.Date("2021-05-01"), date < as.Date("2023-12-19"))
# 
# ggplot(sap_mean, aes(date, sapflow))+geom_line()+labs(x="", y="Daily Sap Flow (kg/day)")+theme_bw()
# 
# sap_mean$month = month(sap_mean$date)
#write_csv(sap_mean, "sapflow_daily.csv")

# sap_daily = read_csv("sapflow_daily.csv")
# ggplot(sap_daily, aes(date, sapflow))+geom_line()+labs(x="", y="Daily Sap Flow (kg/day)")+theme_bw()

# TreeNet data
TN = read_csv("TreeNet_dendro/tn_timeseries.csv")
meta = read_csv("TreeNet_dendro/tn_metadata.csv")

# long-term dendro data for 2 stems
TWD12 = read_csv("TWD_hourly_2014_2022_Nodes_1_2.csv")
Stem_rad12 = read_csv("Stem_radius_hourly_2014_2022_Nodes_1_2.csv")
