# combine and analyze Vetroz data
# data processing below

library(tidyverse)
library(lubridate)
library(gridExtra)

setwd("../../Data/Daten_Vetroz_Lorenz")

## meteo data
VPD = read_csv("Vetroz_Eiche_VPD.csv")
VPD$datetime = as.POSIXct(VPD$MEAS_DATE, format="%d.%m.%Y %H:%M:%S", tz="UTC")
VPD$date = as.Date(VPD$datetime)

# daily mean VPD
VPD_daily = VPD %>% group_by(date) %>% summarize(VPD = mean(VALUE))
VPD_daily$year = year(VPD_daily$date)
VPD_daily$month = month(VPD_daily$date)

## leaf water potential
LWP = read_csv("LWP_gs.csv")
LWP$year = year(LWP$Date)
LWP$LWP_mean = LWP$LWP_mean / -10 # convert to MPa

# compare bagged and unbagged
ggplot(filter(LWP, Date < "2021-10-01", pd_md=="md"), aes(Date, LWP_mean, fill=Method, group=interaction(Date, Method)))+
  geom_point(shape=21, size=2, position=position_dodge(width=5))+
  geom_point(data=filter(LWP, Date < "2021-10-01", pd_md=="pd"), aes(Date, LWP_mean, fill="pre-dawn"), size=2, shape=21)+
  scale_fill_manual(values=c("bagged"="green", "unbagged"="red", "pre-dawn"="black"))+
  facet_wrap(~TreeNr)+theme_bw()+labs(x="2021", y="Leaf Water Potential (MPa)")

LWP_BWP = LWP %>% filter(pd_md=="md") %>% select(Date, TreeNr, LWP_mean, Method) %>% 
  group_by(Date, TreeNr, Method) %>% summarize(LWP=mean(LWP_mean)) %>% 
  pivot_wider(names_from=Method, values_from=LWP) %>% 
  rename(LWP = unbagged, BWP = bagged) %>% mutate(WP_diff = BWP - LWP)

ggplot(LWP_BWP, aes(BWP, LWP))+geom_point()+geom_abline(slope=1, intercept=0)+
  stat_smooth(method="lm")+theme_bw()+
  labs(x="Bagged Leaf Water Potential (MPa)", y="Unbagged Leaf Water Potential (MPa)")

ggplot(LWP_BWP, aes(LWP, WP_diff))+geom_point()+
  stat_smooth()+theme_bw()+
  labs(x="Unbagged Leaf Water Potential (MPa)", y="Difference between bagged and unbagged LWP (MPa)")

# compared against pre-dawn lwp
LWPpd = LWP %>% filter(pd_md=="pd") %>% select(Date, TreeNr, LWP_mean) %>% rename(LWP_pd=LWP_mean)
LWP_BWP = left_join(LWP_BWP, LWPpd)

ggplot(LWP_BWP, aes(LWP_pd, WP_diff))+geom_point()+
  stat_smooth()+theme_bw()+
  labs(x="Pre-dawn Leaf Water Potential (MPa)", y="Difference between bagged and unbagged LWP (MPa)")

ggplot(LWP_BWP, aes(LWP_pd, LWP))+geom_point()+
  theme_bw()+geom_abline(slope=1, intercept=0)+
  labs(x="Pre-dawn Leaf Water Potential (MPa)", y="Midday Leaf Water Potential (MPa)")

# compare against soil water potential
ggplot(filter(LWP, TreeNr %in% c(3, 8)), aes(SMP_mean / 1000, LWP_mean, color=pd_md))+geom_point(size=2)+
  geom_abline(slope=1, intercept=0)+facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Soil Water Potential (MPa)", y = "Leaf Water Potential (MPa)")

# sap flow and leaf water potential
ggplot(filter(LWP, TreeNr %in% c(3, 8)), aes(LWP_mean, sapflow, shape=Method, color=VPD_kPa))+geom_point(size=2)+
  facet_wrap(~TreeNr)+theme_bw()+scale_color_viridis_c(name="VPD (kPa)")+
  theme(legend.position="inside", legend.position.inside=c(.6,.7))+
  labs(x="Leaf Water Potential (MPa)", y = "Sapflow", color="VPD (kPa)")

# TWD and LWP
ggplot(filter(LWP, Method=="unbagged"), aes(TWD, LWP_mean, color=pd_md))+geom_point()+
  geom_smooth(method="lm")+facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Tree Water Deficit", y = "Leaf Water Potential (MPa)")

# gs and LWP
ggplot(LWP, aes(LWP_mean, gs, shape=Method, color=VPD_kPa))+geom_point(size=2)+facet_wrap(~TreeNr)+
  theme_bw()+scale_color_viridis_c(name="VPD (kPa)")+
  theme(legend.position="inside", legend.position.inside=c(.05,.25))+
  labs(x = "Leaf Water Potential (MPa)", y="Stomatal Conductance (mol m^-2 s^-1)")


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

SWP1cm_daily$SWP[SWP1cm_daily$SWP < -15000] = NA # cap at -15000 kPa

SWP_daily = rbind(SWP_daily, SWP1cm_daily)

SWP_daily$year = year(SWP_daily$date)
SWP_daily$month = month(SWP_daily$date)

# derive weighted-average SWP
depth_weights = data.frame(depth=c(1, 20, 80, 110, 160), weight=c(5/160, 45/160, 40/160, 40/160, 30/160))
SWP_mean = SWP_daily %>% left_join(depth_weights) %>% group_by(date, year, month) %>% 
  summarize(SWP = sum(SWP * weight, na.rm=T)) %>% filter(date < as.Date("2023-01-01"))
ggplot(SWP_mean, aes(date, SWP))+geom_line()+
  labs(x="", y="Weighted-average SWP (kPa)")+theme_bw()

# remove NAs for output
SWP_daily = na.omit(SWP_daily)
#write_csv(SWP_daily, "SWP_daily.csv")

#SWP_daily = read_csv("SWP_daily.csv")

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

# normalize sap flow by treatment-specific max sap flow
sap_max = sap_mean %>% group_by(sensor_loc) %>% summarize(max_sap=max(sapflow))
sap_mean = sap_mean %>% left_join(sap_max) %>% mutate(sapflow_norm = sapflow / max_sap) %>% select(-max_sap)

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
    group_by(date, tree_name, tree_id, sensor_loc) %>% summarize(TWD_pd=min(TWD, na.rm=T), .groups="drop")
  
  # aggregate to daily with mid-day
  den_md = den_long %>% filter(hour >= 12, hour <= 17) %>% 
    group_by(date, tree_name, tree_id, sensor_loc) %>% summarize(TWD_md=max(TWD, na.rm=T), .groups="drop")
  
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

# correlate predawn LWP against TWD per tree
LWP_twd_pd = LWP %>% filter(pd_md == "pd") %>%
  select(date = Date, TreeNr, LWP_mean, pd_md) %>%
  inner_join(TWD_norm %>% filter(sensor_loc == "stem") %>%
               mutate(TreeNr = as.numeric(tree_id)) %>%
               select(date, TreeNr, TWD_pd))

TWD_lwp_coef_pd = LWP_twd_pd %>% group_by(TreeNr) %>% nest() %>%
  mutate(model = map(data, ~ lm(LWP_mean ~ TWD_pd, data = .x)),
         TWD_lwp_intercept_pd = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
         TWD_lwp_slope_pd = map_dbl(model, ~ coef(.x)[["TWD_pd"]])) %>% ungroup() %>% 
  select(TreeNr, TWD_lwp_intercept_pd, TWD_lwp_slope_pd)

# correlate midday LWP against TWD per tree
LWP_twd_md = LWP %>% filter(pd_md == "md", Method=="unbagged") %>% 
  group_by(Date, TreeNr) %>% summarize(LWP_mean=min(LWP_mean)) %>%
  select(date = Date, TreeNr, LWP_mean) %>%
  inner_join(TWD_norm %>% filter(sensor_loc == "stem") %>%
               mutate(TreeNr = as.numeric(tree_id)) %>%
               select(date, TreeNr, TWD_md))

TWD_lwp_coef_md = LWP_twd_md %>% group_by(TreeNr) %>% nest() %>%
  mutate(model = map(data, ~ lm(LWP_mean ~ TWD_md, data = .x)),
         TWD_lwp_intercept_md = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
         TWD_lwp_slope_md = map_dbl(model, ~ coef(.x)[["TWD_md"]])) %>% ungroup() %>% 
         select(TreeNr, TWD_lwp_intercept_md, TWD_lwp_slope_md)

# combine for plotting
LWP_TWD = rbind(LWP_twd_md %>% rename(TWD=TWD_md) %>% mutate(pd_md="md"),
                LWP_twd_pd %>% rename(TWD=TWD_pd))

ggplot(LWP_TWD, aes(TWD, LWP_mean, color=pd_md))+geom_point()+
  geom_smooth(method="lm")+facet_wrap(~TreeNr)+theme_bw()+
  labs(x="Tree Water Deficit", y = "Leaf Water Potential (MPa)")

# predict predawn and midday LWP from TWD for stem sensors
TWD_norm = TWD_norm %>% mutate(TreeNr = as.numeric(tree_id)) %>%
  left_join(TWD_lwp_coef_pd) %>% left_join(TWD_lwp_coef_md) %>%
  mutate(LWP_md = if_else(sensor_loc == "stem",
                          TWD_lwp_intercept_md + TWD_lwp_slope_md * TWD_md, NA),
         LWP_pd = if_else(sensor_loc == "stem",
                          TWD_lwp_intercept_pd + TWD_lwp_slope_pd * TWD_pd, NA)) %>%
  select(-c(TreeNr, TWD_lwp_intercept_md, TWD_lwp_slope_md, TWD_lwp_intercept_pd, TWD_lwp_slope_pd))

# exclude pre-dawn predictions outside of calibration range
TWD_min = LWP_twd_pd %>% group_by(TreeNr) %>% summarize(TWD_min = min(TWD_pd))

TWD_norm = TWD_norm %>% mutate(TreeNr = as.numeric(tree_id)) %>% left_join(TWD_min) %>% 
  mutate(LWP_pd = if_else(TWD_pd < TWD_min, NA, LWP_pd)) %>% select(-TreeNr, -TWD_min)

# plots
ggplot(TWD_norm, aes(date, TWD_pd, color=tree_id, group=tree_name))+geom_line()+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(date, LWP_pd, color=tree_id, group=tree_name))+geom_line()+
  theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(date, TWD_md, color=tree_id, group=tree_name))+geom_line()+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()+guides(color="none")

ggplot(TWD_norm, aes(date, LWP_md, color=tree_id, group=tree_name))+geom_line()+
  theme_bw()+guides(color="none")

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
  group_by(date, sensor_loc) %>% summarize_all(list(mean), na.rm=T)

ggplot(TWD_daily, aes(date, TWD_pdn))+geom_line(color="green")+
  geom_line(aes(date, MDS_norm), color="black")+
  facet_wrap(~sensor_loc, ncol=1, scales="free")+theme_bw()

# exclude winter months
TWD_daily$year = year(TWD_daily$date)
TWD_daily$month = month(TWD_daily$date)
TWD_daily[TWD_daily$month < 5 | TWD_daily$month > 11, c("TWD_pd", "TWD_pdn", "TWD_md", "MDS", "MDS_norm", "LWP_pd", "LWP_md")] = NA

#write_csv(TWD_daily, "TWD_daily.csv", na="")
#TWD_daily = read_csv("TWD_daily.csv")

# mutual plotting data frames
sap_plot = sap_daily %>% mutate(ddate=as.Date(paste(2000, month, day(date), sep="-"))) %>% 
  filter(date < "2022-11-01", month > 4, month < 11)
TWD_plot = TWD_daily %>% mutate(ddate=as.Date(paste(2000, month(date), day(date), sep="-"))) %>%
  filter(date >= "2021-05-07", date < "2022-11-01", month > 4, month < 11, sensor_loc == "stem")
LWP_pd_plot = LWP %>% rename(LWP=LWP_mean) %>% filter(pd_md=="pd") %>% group_by(Date, year) %>% 
  summarize(LWP_mean=mean(LWP), LWP_sd=sd(LWP)) %>%
  mutate(ddate=as.Date(paste(2000, month(Date), day(Date), sep="-")))
SWP_plot = SWP_mean %>% mutate(ddate=as.Date(paste(2000, month(date), day(date), sep="-"))) %>% 
  filter(date >= "2021-05-07", date < "2022-11-01", month>4, month<11)

p_swp = ggplot(SWP_plot, aes(ddate, SWP/1000, color="Weighted SWP"))+geom_line()+
  geom_line(data=TWD_plot, aes(ddate, LWP_pd, color="TWD-derived LWP"), inherit.aes=F)+
  geom_pointrange(data=LWP_pd_plot, aes(x=ddate, y=LWP_mean, ymin=LWP_mean-LWP_sd, ymax=LWP_mean+LWP_sd, color="Measured LWP"), inherit.aes=F)+
  facet_wrap(~year)+theme_bw()+labs(x="", y="SWP / pre-dawn LWP (MPa)")+
  scale_colour_manual(name="Data", values=c("Weighted SWP"="black","Measured LWP"="red", "TWD-derived LWP"="blue"), 
                      breaks=c("Measured LWP", "TWD-derived LWP", "Weighted SWP"))+
  theme(legend.position="inside", legend.position.inside=c(0.1,0.2))

p_sapn = ggplot(sap_plot, aes(ddate, sapflow_norm*100, color=sensor_loc))+geom_line()+
  labs(x="", y="Normalized daily Sap Flow (%)", color="Location")+theme_bw()+
  facet_wrap(~year)+theme(legend.position="inside", legend.position.inside=c(0.45,0.75))

grid.arrange(p_sapn, p_swp)


# stacked plots
p_sap = ggplot(sap_plot)+geom_line(aes(ddate, sapflow, color=sensor_loc))+
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
  geom_line(data=TWD_plot, aes(ddate, LWP_pd), color="black", inherit.aes=F)+
  geom_pointrange(data=LWP_pd_plot, aes(x=ddate, y=LWP_mean, ymin=LWP_mean-LWP_sd, ymax=LWP_mean+LWP_sd), color="blue", inherit.aes=F)+
  labs(x="", y="SWP / pre-dawn LWP (MPa)", color="Depth (cm)")+theme_bw()+
  scale_x_date(date_breaks="1 month", date_labels="%b")+facet_wrap(~year)+ylim(-4.5, 0)+
  theme(legend.position="inside", legend.position.inside=c(0.05,0.4))

grid.arrange(p_sap, p_twd, p_wp)


# compare TWD against SWP

TWD_SWP = left_join(TWD_daily, SWP_mean)
TWD_SWP = left_join(TWD_SWP, VPD_daily)

ggplot(TWD_SWP, aes(SWP/1000, TWD_pd, color=VPD))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+scale_color_viridis_c(name="VPD (kPa)")+
  theme(legend.position="inside", legend.position.inside=c(.4,.8))+
  labs(x="Weighted Mean SWP (MPa)", y="Pre-dawn TWD")

ggplot(TWD_SWP, aes(SWP/1000, TWD_pdn))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+
  labs(x="Weighted Mean SWP (MPa)", y="Normalized pre-dawn TWD")

ggplot(TWD_SWP, aes(SWP/1000, TWD_md))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+
  labs(x="Weighted Mean SWP (MPa)", y="Midday TWD")

ggplot(TWD_SWP, aes(SWP/1000, MDS_norm))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+
  labs(x="Weighted Mean SWP (MPa)", y="Normalized Maximum Daily Shrinkage")

ggplot(TWD_SWP, aes(SWP/1000, LWP_pd))+geom_point()+
  geom_abline(slope=1, intercept=0)+theme_bw()+
  labs(x="Weighted Mean SWP (MPa)", y="Pre-dawn TWD-derived LWP")

ggplot(TWD_SWP, aes(SWP/1000, LWP_md))+geom_point()+
  geom_abline(slope=1, intercept=0)+theme_bw()+
  labs(x="Weighted Mean SWP (MPa)", y="Midday TWD-derived LWP")

# compare sap against swp

sap_swp = left_join(sap_daily, SWP_mean)
sap_swp = left_join(sap_swp, VPD_daily)

ggplot(filter(sap_swp, month>4, month<11), aes(SWP/1000, sapflow))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+
  labs(x="Weighted Mean SWP (MPa)", y="Daily Sap Flow (kg/day)")

ggplot(filter(sap_swp, month>4, month<11), aes(SWP/1000, sapflow_norm*100, color=VPD))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+scale_color_viridis_c(name="VPD (kPa)")+
  theme(legend.position="inside", legend.position.inside=c(.1,.8))+
  labs(x="Weighted Mean SWP (MPa)", y="Normalized daily Sap Flow (%)")

# sap against twd

sap_twd = left_join(sap_daily, TWD_daily)
sap_twd = left_join(sap_twd, VPD_daily)

ggplot(sap_twd, aes(LWP_pd, sapflow))+geom_point()+theme_bw()+
  labs(x="Pre-dawn TWD-derived LWP", y="Daily Sap Flow (kg/day)")

ggplot(filter(sap_twd, month>4, month<11, year<2023), aes(TWD_pdn, sapflow_norm*100, color=VPD))+geom_point()+
  facet_wrap(~sensor_loc)+theme_bw()+scale_color_viridis_c(name="VPD (kPa)")+
  theme(legend.position="inside", legend.position.inside=c(.4,.8))+
  labs(x="Normalized pre-dawn TWD", y="Normalized Daily Sap Flow (%)")



