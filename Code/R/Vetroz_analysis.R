# combine and analyze Vetroz data

library(tidyverse)
library(lubridate)

setwd("../../Data/Daten_Vetroz_Lorenz")

## soil water potential data
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
