# analyze dendrometer data in Pfynwald

library(tidyverse)
library(lubridate)

setwd("../../Data/Pfyn")

# dendrometer data from Zweifel et al. (2020)

den = read_tsv("pfynwald_dendro.tab", skip=43)

# lengthen data frame and separate column names into multiple variables
den_long = den %>% pivot_longer(cols=-1, names_to=c(".value", "tree_id","scenario"), 
                              names_pattern="^(.*?) \\[.*?\\] \\(at tree (\\d+),\\s*([^\\)]+)\\)")

# remove plot from scenario names
den_long$scenario <- den_long$scenario %>% str_replace(" plot$", "") %>% str_trim()

# remove spaces from column names
den_long = den_long %>% rename(SRC='Tree SRC', GRO='Tree GI', datetime='Date/Time')

ggplot(den_long, aes(x=datetime, y=TWD, color=as.factor(tree_id)))+geom_line()+
  facet_wrap(~tree_id)+guides(color="none")

# formatting
den_long = na.omit(den_long)
den_long$date = as.Date(den_long$datetime)
den_long$hour = hour(den_long$datetime)

# aggregate to daily with pre-dawn
den_pd = den_long %>% filter(hour <= 6) %>% 
  group_by(date, tree_id, scenario) %>% summarize(TWD_pd=min(TWD, na.rm=T))

# aggregate to daily with mid-day
den_md = den_long %>% filter(hour >= 12, hour <= 17) %>% 
  group_by(date, tree_id, scenario) %>% summarize(TWD_md=max(TWD, na.rm=T))

den_daily = full_join(den_pd, den_md)

# calculate maximum daily shrinkage
den_daily$MDS = den_daily$TWD_md - den_daily$TWD_pd

# calculate 99th percentile MDS over entire time series
den_mds = den_daily %>% mutate(month=month(date)) %>% 
  filter(month>3, month<11, MDS>0) %>%  # only accept positive shrinkage during summer months
  group_by(tree_id) %>% summarize(mds_max = quantile(MDS, .99))

# join max MDS to daily data frame
den_daily = left_join(den_daily, den_mds)

# normalize pre-dawn TWD and MDS with max MDS
den_daily$TWD_pdn = den_daily$TWD_pd / den_daily$mds_max
den_daily$MDS_norm = den_daily$MDS / den_daily$mds_max

# plots
ggplot(den_daily, aes(x=date, y=TWD_pd, color=as.factor(tree_id)))+geom_line()+
  facet_wrap(~tree_id)+guides(color="none")

ggplot(den_daily, aes(x=date, y=TWD_pdn, color=as.factor(tree_id)))+geom_line()+
  facet_wrap(~tree_id)+guides(color="none")

ggplot(den_daily, aes(TWD_pdn, MDS_norm, color=as.factor(tree_id)))+geom_point()+
  guides(color="none")+theme_bw()

ggplot(filter(den_daily, tree_id=="125"), aes(date, TWD_pdn))+geom_line(color="green")+
  geom_line(aes(date, MDS_norm), color="black")+theme_bw()

# group by scenario
den_meta = den_daily %>% filter(tree_id != "125") %>% na.omit() %>% 
  group_by(date, scenario) %>% summarize_at(vars(TWD_pd, TWD_pdn, MDS_norm), list(mean))
den_meta$year = year(den_meta$date)

ggplot(filter(den_meta, scenario!="irrigation"), aes(date, TWD_pd, color=as.factor(scenario)))+geom_line()+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()+
  labs(x="", y="Pre-dawn TWD", color="Scenario")

ggplot(den_meta, aes(date, TWD_pd, color=as.factor(scenario)))+geom_line()+
  theme_bw()+theme(legend.position="inside",legend.position.inside=c(.9,.85))+
  labs(x="", y="Pre-dawn TWD", color="Scenario")+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))


den_plot = merge(seq(as.Date("2011-01-01"), as.Date("2017-12-31")), c("control","irrigation stop"))
colnames(den_plot) = c("date","scenario")
den_plot = left_join(den_plot, filter(den_meta, scenario!="irrigation"))

ggplot(filter(den_plot, year>2012, year<2016), aes(date, TWD_pd, color=as.factor(scenario)))+geom_line()+
  theme_bw()+theme(legend.position="inside",legend.position.inside=c(.9,.85), 
                   legend.title=element_text(size=12), legend.text=element_text(size=12),
                   axis.text=element_text(size=12), axis.title=element_text(size=14))+
  labs(x="", y=expression("Pre-dawn TWD ("*mu*"m)"), color="Scenario")+
  scale_color_manual(values=c("#E69F00","#009E73"))

# write_csv(den_meta, "Pfyn_twd_2011_17.csv")



# VPDrought 2024 data

TN_dendro = read_csv("TreeNet/tn_timeseries_Pfyn_dendro.csv")
TN_meta = read_csv("TreeNet/tn_metadata_Pfyn_dendro.csv")

# select only VPDrought trees
TN_vpd_meta = TN_meta %>% filter(series_start=="2023-08-01")

# extract scenario names
TN_vpd_meta$scenario = gsub("^pfynwald-([a-z]+_?[a-z]+)_.*", "\\1", TN_vpd_meta$measure_point)

# filter timeseries for VPDrought trees and combine with scenario name
TN_vpden = TN_dendro %>% select(-c(frost, flags, gro_start, gro_end)) %>% 
  filter(series_id %in% TN_vpd_meta$series_id) %>% 
  left_join(select(TN_vpd_meta, series_id, scenario))

ggplot(TN_vpden, aes(ts, twd, color=as.factor(series_id)))+geom_line()+
  facet_wrap(~scenario)+theme_bw()+guides(color="none")

ggplot(TN_vpden, aes(ts, gro_yr, color=as.factor(series_id)))+geom_line()+
  facet_wrap(~scenario)+theme_bw()+guides(color="none")


# formatting
TN_vpden$date = as.Date(TN_vpden$ts)
TN_vpden$hour = hour(TN_vpden$ts)

# aggregate to daily with pre-dawn
vpden_pd = TN_vpden %>% filter(hour <= 6) %>% 
  group_by(date, series_id, scenario) %>% summarize(twd_pd=min(twd, na.rm=T))

# aggregate to daily with mid-day
vpden_md = TN_vpden %>% filter(hour >= 12, hour <= 17) %>% 
  group_by(date, series_id, scenario) %>% summarize(twd_md=max(twd, na.rm=T))

TN_vpden_daily = full_join(vpden_pd, vpden_md)

# calculate maximum daily shrinkage
TN_vpden_daily$MDS = TN_vpden_daily$twd_md - TN_vpden_daily$twd_pd

# calculate 99th percentile MDS over entire time series
vpden_mds = TN_vpden_daily %>% mutate(month=month(date)) %>% 
  filter(month>3, month<11, MDS>0) %>%  # only accept positive shrinkage during summer months
  group_by(series_id) %>% summarize(mds_max = quantile(MDS, .99))

# join max MDS to daily data frame
TN_vpden_daily = left_join(TN_vpden_daily, vpden_mds)

# normalize pre-dawn TWD and MDS with max MDS
TN_vpden_daily$twd_pdn = TN_vpden_daily$twd_pd / TN_vpden_daily$mds_max
TN_vpden_daily$MDS_norm = TN_vpden_daily$MDS / TN_vpden_daily$mds_max

ggplot(TN_vpden_daily, aes(date, twd_pd, color=as.factor(series_id)))+geom_line()+
  facet_wrap(~scenario)+guides(color="none")+theme_bw()

ggplot(filter(TN_vpden_daily, !(series_id %in% c(1580,1582))), aes(date, twd_pdn, color=as.factor(series_id)))+geom_line()+
  facet_wrap(~scenario)+guides(color="none")+theme_bw()


# group by scenario
TN_vpden_meta = TN_vpden_daily %>% filter(!(series_id %in% c(1580,1582))) %>% 
  group_by(date, scenario) %>% summarize_at(vars(twd_pd, twd_pdn, MDS_norm), list(mean))

TN_vpden_meta$scenario = factor(TN_vpden_meta$scenario, 
                            levels=c("control","roof","roof_vpd","irrigation","irrigation_vpd"))

ggplot(TN_vpden_meta, aes(date, twd_pd, color=as.factor(scenario)))+geom_line()+
  theme_bw()+theme(legend.position="inside",legend.position.inside=c(.8,.8))+
  labs(x="", y="Pre-dawn TWD", color="Scenario")

ggplot(TN_vpden_meta, aes(date, twd_pdn, color=as.factor(scenario)))+geom_line()+
  theme_bw()+#theme(legend.position="inside",legend.position.inside=c(.9,.9))+
  labs(x="", y="Normalized Pre-dawn TWD")+guides(color="none")


# older data from Richard

den_con = readRDS("data for Pfynwald/older/dendro/PFY_DR-TWD.RDS")
den_irr = readRDS("data for Pfynwald/older/dendro/PFY_WE-TWD.RDS")

den_con = den_con %>% select(-c(frost, flags, version)) %>% na.omit()
den_irr = den_irr %>% select(-c(frost, flags, version)) %>% na.omit()

