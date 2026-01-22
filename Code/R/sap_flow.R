# analyze sap flow data in Pfynwald

library(tidyverse)
library(lubridate)
library(readxl)
library(zoo)

## common data and functions

source("../../Data/Pfyn/data for Pfynwald/Pfyn_sf_helpers.R") # from Richard Peters for VPDrought

sapwood_bin_size = 300

# inventoried trees
Pfyn_inv = read_excel("../../Data/Pfyn/Pfyn_DBH_2023.xlsx")
Pfyn_inv = filter(Pfyn_inv, Plot!=0)

Pfyn_inv$SWD = mapply(function(dbh, treat) {SWD(dbh, treat)}, Pfyn_inv$DBH, Pfyn_inv$Treatment) # sapwood depth
Pfyn_inv$BDD = BDD(Pfyn_inv$DBH) # bark thickness
Pfyn_inv$sapwood_area = mapply(function(dbh, swd, bdd){calculate_sapwood_area(dbh, swd, bdd)}, Pfyn_inv$DBH, Pfyn_inv$SWD, Pfyn_inv$BDD)

ggplot(Pfyn_inv, aes(DBH))+geom_histogram()+facet_wrap(~Treatment, scales="free_y")
ggplot(Pfyn_inv, aes(sapwood_area))+geom_histogram()+facet_wrap(~Treatment, scales="free_y")
ggplot(Pfyn_inv, aes(sapwood_area))+geom_histogram()+facet_wrap(~Treatment_VPD, scales="free_y")
hist(filter(Pfyn_inv, Treatment_VPD=="Control")$sapwood_area)

# group sapwood area into bins per plot and treatment and calculate proportional area
Pfyn_inv$sapwood_area_bin = cut_width(Pfyn_inv$sapwood_area, sapwood_bin_size, boundary=0)
Pfyn_inv_trt = Pfyn_inv %>% group_by(Treatment, sapwood_area_bin) %>% summarize(n=n()) %>% mutate(prop = n/sum(n))

# sapwood area by plot
Pfyn_inv_sum_VPD = Pfyn_inv %>% group_by(Plot, Treatment_VPD) %>% summarize(sapwood_area_total = sum(sapwood_area)) # in cm^2
Pfyn_inv_sum = Pfyn_inv %>% group_by(Plot, Treatment) %>% summarize(sapwood_area_total = sum(sapwood_area)) # in cm^2

# plot area
# irrigation stop plot areas in m^2: p2=475, p3=450, p6=428, p7=357
# roof area seems to be about 16x16 m = 256 m^2
Pfyn_inv_sum_VPD$plot_area = c(1000, 475, 525, 488, 256, 256, 638, 256, 256, 428, 572, 357, 643, 744, 256, 256)
Pfyn_inv_sum$plot_area = c(1000, 525, 475, 1000, 1000, 572, 428, 643, 357, 1000)

# treatment sapwood area
sap_area_cont_vpd = sum(filter(Pfyn_inv_sum_VPD, Treatment_VPD=="Control")$sapwood_area_total)
sap_area_cont = sum(filter(Pfyn_inv_sum, Treatment=="Control")$sapwood_area_total)
sap_area_irr = sum(filter(Pfyn_inv_sum_VPD, Treatment_VPD=="Irrigation")$sapwood_area_total)
sap_area_drt = sum(filter(Pfyn_inv_sum_VPD, Treatment_VPD=="Roof")$sapwood_area_total)
sap_area_irst = sum(filter(Pfyn_inv_sum, Treatment=="Irrigation stop")$sapwood_area_total)

# treatment area
area_cont_vpd = sum(filter(Pfyn_inv_sum_VPD, Treatment_VPD=="Control")$plot_area)
area_cont = sum(filter(Pfyn_inv_sum, Treatment=="Control")$plot_area)
area_irr = sum(filter(Pfyn_inv_sum_VPD, Treatment_VPD=="Irrigation")$plot_area)
area_drt = sum(filter(Pfyn_inv_sum_VPD, Treatment_VPD=="Roof")$plot_area)
area_irst = sum(filter(Pfyn_inv_sum, Treatment=="Irrigation stop")$plot_area)

# TreeNet metadata
TN_ctr_meta = read_csv("../../Data/Pfyn/TreeNet/tn_metadata_Pfyn_sap_control.csv")
TN_irr_meta = read_csv("../../Data/Pfyn/TreeNet/tn_metadata_Pfyn_sap_irrigation.csv")

# calculate sapwood characteristics of measured trees

sap_profile = function(df, sap_bin) {

  df$SWD = mapply(function(dbh, treat) {SWD(dbh, treat)}, df$tree_dbh, "Control")
  df$BDD = BDD(df$tree_dbh)
  df$sapwood_area = mapply(function(dbh, swd, bdd){calculate_sapwood_area(dbh, swd, bdd)}, df$tree_dbh, df$SWD, df$BDD)
  df$sapwood_area_bin = cut_width(df$sapwood_area, sap_bin, boundary=0)
  df$rad_cor = sapply(df$SWD, rad_cor) # radial profile correction
  
  return(df)
}

# control scenario
TN_sf_ctr = TN_ctr_meta %>% select(tree_name, tree_dbh) %>% 
  filter(tree_name %in% c(109, 110, 124))

TN_sf_ctr = sap_profile(TN_sf_ctr, sapwood_bin_size)

# irrigation stop scenario
TN_sf_irst = TN_irr_meta %>% select(tree_name, tree_dbh) %>% 
  filter(tree_name %in% c(274, 275, 276))

TN_sf_irst = sap_profile(TN_sf_irst, sapwood_bin_size)

# irrigation scenario
TN_sf_irr = TN_irr_meta %>% select(tree_name, tree_dbh) %>% 
  filter(tree_name %in% c(246, 247, 250))

TN_sf_irr = sap_profile(TN_sf_irr, sapwood_bin_size)


## sap flow data from Zweifel et al. (2020)

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
sf_daily = sf_long %>% mutate(datehour=floor_date(datetime, "1 hour")) %>% 
  group_by(datehour, tree_id, scenario) %>% summarize(sf=mean(Tree_SR, na.rm=T), .groups="drop") %>% 
  mutate(date=as.Date(datehour)) %>% select(-datehour) %>% 
  group_by(date, tree_id, scenario) %>% summarize_at(vars(sf), list(sum))
sf_daily = na.omit(sf_daily)

ggplot(sf_daily, aes(date, sf, color=as.factor(tree_id)))+geom_line(alpha=.8)+guides(color="none")+
  facet_wrap(~scenario, ncol=1)+theme_bw()+labs(x="",y="Mean Daily Sap Flux [kg/day]", color="Tree ID")

ggplot(filter(sf_daily, scenario!="irrigation", tree_id!="125"), aes(date, sf, color=as.factor(scenario)))+
  stat_summary(geom="line", fun=mean)+theme_bw()+labs(x="", y="Mean Daily Sap Flow [kg/day]", color="Scenario")

# group by scenario
sf_meta = sf_daily %>% filter(tree_id != "125", tree_id!="274") %>% group_by(date, scenario) %>% summarize(sf=mean(sf, na.rm=T))
sf_meta$year = year(sf_meta$date)

ggplot(filter(sf_meta, scenario!="irrigation"), aes(date, sf, color=as.factor(scenario)))+geom_line()+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()+
  theme(legend.position="inside",legend.position.inside=c(.1,.95))+
  labs(x="", y="Mean Daily Sap Flow [kg/h]", color="Scenario")

# write_csv(sf_meta, "Pfyn_sap_2011_17.csv")

## upscaling procedure

# control scenario
sf_ctr = sf_daily %>% filter(tree_id != "125", scenario=="control")

# combine with measurements for weighting
sf_ctr = left_join(mutate(sf_ctr, tree_id=as.numeric(tree_id)), rename(TN_sf_ctr, tree_id=tree_name))

# sap flow normalized by sapwood area and corrected for radial profile
sf_ctr$sf_norm = sf_ctr$sf / sf_ctr$sapwood_area * sf_ctr$rad_cor

# group by size class and add inventory statistics by treatment
sf_ctr_bin = sf_ctr %>% group_by(date, sapwood_area_bin) %>% summarize_at(vars(sf, sf_norm), list(mean))
sf_ctr_bin = left_join(sf_ctr_bin, filter(Pfyn_inv_trt, Treatment=="Control"))

# upscale by sapwood area
sf_ctr_daily = sf_ctr_bin %>% group_by(date) %>% summarize(Tr=sum(sf_norm*prop))
sf_ctr_daily$Tr = sf_ctr_daily$Tr * sap_area_cont / area_cont
sf_ctr_daily$Tr_rm = rollmean(sf_ctr_daily$Tr, 14, fill=NA) # rolling mean

ggplot(sf_ctr_daily, aes(date, Tr))+geom_line()

sf_ctr_yr = sf_ctr_daily %>% mutate(year=year(date)) %>% 
  group_by(year) %>% summarize(Tr=sum(Tr))

# irrigation stop scenario
sf_irst = sf_daily %>% filter(scenario=="irrigation stop", tree_id != "274")

# combine with measurements for weighting
sf_irst = left_join(mutate(sf_irst, tree_id=as.numeric(tree_id)), rename(TN_sf_irst, tree_id=tree_name))

# sap flow normalized by sapwood area and corrected for radial profile
sf_irst$sf_norm = sf_irst$sf / sf_irst$sapwood_area * sf_irst$rad_cor

# group by size class and add inventory statistics by treatment
sf_irst_bin = sf_irst %>% group_by(date, sapwood_area_bin) %>% summarize_at(vars(sf, sf_norm), list(mean))
sf_irst_bin = left_join(sf_irst_bin, filter(Pfyn_inv_trt, Treatment=="Irrigation stop"))

# upscale by sapwood area
sf_irst_daily = sf_irst_bin %>% group_by(date) %>% summarize(Tr=sum(sf_norm*prop))
sf_irst_daily$Tr = sf_irst_daily$Tr * sap_area_irst / area_irst
sf_irst_daily$Tr_rm = rollmean(sf_irst_daily$Tr, 14, fill=NA) # rolling mean

ggplot(sf_irst_daily, aes(date, Tr))+geom_line()

sf_irst_yr = sf_irst_daily %>% mutate(year=year(date)) %>% 
  group_by(year) %>% summarize(Tr=sum(Tr))

# irrigation scenario
sf_irr = sf_daily %>% filter(scenario=="irrigation")

# combine with measurements for weighting
sf_irr = left_join(mutate(sf_irr, tree_id=as.numeric(tree_id)), rename(TN_sf_irr, tree_id=tree_name))

# sap flow normalized by sapwood area and corrected for radial profile
sf_irr$sf_norm = sf_irr$sf / sf_irr$sapwood_area * sf_irr$rad_cor

# group by size class and add inventory statistics by treatment
sf_irr_bin = sf_irr %>% group_by(date, sapwood_area_bin) %>% summarize_at(vars(sf, sf_norm), list(mean))
sf_irr_bin = left_join(sf_irr_bin, filter(Pfyn_inv_trt, Treatment=="Irrigation"))

# upscale by sapwood area
sf_irr_daily = sf_irr_bin %>% group_by(date) %>% summarize(Tr=sum(sf_norm*prop))
sf_irr_daily$Tr = sf_irr_daily$Tr * sap_area_irr / area_irr
sf_irr_daily$Tr_rm = rollmean(sf_irr_daily$Tr, 14, fill=NA) # rolling mean

ggplot(sf_irr_daily, aes(date, Tr))+geom_line()

sf_irr_yr = sf_irr_daily %>% mutate(year=year(date)) %>% 
  group_by(year) %>% summarize(Tr=sum(Tr))

# combine and output into single file
trans_comp = rbind(mutate(sf_ctr_daily, scen="Control"), mutate(sf_irst_daily, scen="Irrigation stop"))
trans_comp = rbind(trans_comp, mutate(sf_irr_daily, scen="Irrigation"))

ggplot(trans_comp, aes(date, Tr, color=scen))+geom_line()+theme_bw()+
  labs(x="", y="Transpiration (mm/day)", color="Treatment")

tr_plot = merge(seq(as.Date("2011-01-01"), as.Date("2017-12-31"), by=1), c("Control","Irrigation","Irrigation stop"))
colnames(tr_plot) = c("date","scen")
tr_plot = left_join(tr_plot, trans_comp)

ggplot()+geom_rect(data=filter(irr, year %in% c(2011, 2012, 2013, 2016, 2017)), aes(xmin=on, xmax=off, ymin=-Inf, ymax=Inf), alpha=.5, fill="lightblue")+
  geom_point(data=tr_plot, aes(date, Tr, color=scen), inherit.aes=F)+
  geom_line(data=tr_plot, aes(date, Tr_rm, color=scen), linewidth=1.1, inherit.aes=F)+
  labs(x="", y="Transpiration (mm/day)", color="Treatment")+theme_bw()+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))+
  ggtitle("Upscaled Transpiration from Sap Flow Data in Pfynwald")+
  theme(legend.position="inside", legend.position.inside=c(.57,.82), plot.title=element_text(hjust=.5, size=16),
        axis.text.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text.y=element_text(size=12), 
        legend.text=element_text(size=12), legend.title=element_text(size=14))

#write_csv(trans_comp, "../../Data/Pfyn/Pfyn_trans_2011_17.csv")

# difference between control and irrigation stop
ctr_irst_comp = inner_join(rename(sf_ctr_daily, ctr=Tr), rename(sf_irst_daily, irst=Tr))
ctr_irst_comp$year = year(ctr_irst_comp$date)

ctr_irst_comp$diff = ctr_irst_comp$irst - ctr_irst_comp$ctr

ggplot(ctr_irst_comp, aes(date, diff))+geom_line()+geom_hline(yintercept=0)+
  facet_wrap(~year, ncol=1, scales="free_x")+theme_bw()

# 2021 - 2022 data

sap = readRDS("../../Data/Pfyn/data for Pfynwald/2021-2022/PFY_sfd_cleaned.Rds")
sap$tree = factor(sap$tree)

# remove outliers
sap$.annotation = ifelse(is.na(sap$.annotation), "OK", sap$.annotation)
sap = filter(sap, .annotation!="Outlier")

sap = select(sap, -c(.dcrkey, .annotation, selection_count, k))

ggplot(sap, aes(timestamp, sfd, color=tree))+geom_line()+facet_wrap(~tree)+
  guides(color="none")+theme_bw()

# group to daily
sap_daily = sap %>% mutate(datehour=floor_date(timestamp, "1 hour")) %>% 
  group_by(datehour, tree) %>% summarize_all(list(mean), .groups="drop") %>% 
  mutate(date=as.Date(datehour)) %>% select(-c(datehour, timestamp)) %>% 
  group_by(date, tree) %>% summarize_at(vars(sfd,q025,q975), list(sum))

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

# upscaling

sap_ctr = filter(sap_daily_test, meta=="Control")
sap_ctr$tree_id = if_else(sap_ctr$tree=="Pfynwald_Control_02_15_ch6", 109, 124)

# combine with measurements for weighting
sap_ctr = left_join(sap_ctr, rename(TN_sf_ctr, tree_id=tree_name))

# add inventory statistics by treatment
sap_ctr = left_join(sap_ctr, filter(Pfyn_inv_trt, Treatment=="Control"))

# upscale by sapwood area
sap_ctr_daily = sap_ctr %>% group_by(date) %>% summarize(Tr=sum(sfd*rad_cor*prop))
sap_ctr_daily$Tr = sap_ctr_daily$Tr * sap_area_cont / area_cont / 1000

ggplot(sap_ctr_daily, aes(date, Tr))+geom_line()


## VPDrought 2024-25 data

#sap_vpd = readRDS("../../Data/Pfyn/data for Pfynwald/2024/VPDrought_SF_L3_2025-02-02.RDS")
sap_vpd = readRDS("../../Data/Pfyn/data for Pfynwald/2024/VPDrought_sfden_2025-06-18.RDS")
#sap_vpd25 = readRDS("../../Data/Pfyn/data for Pfynwald/2024/VPDrought_SF_L3_2025-12-01.RDS")

sap_vpd$Date.Time = as.POSIXct(sap_vpd$Date.Time, tz="CET", format="%Y-%m-%d %H:%M:%S")
sap_vpd$date = as.Date(sap_vpd$Date.Time, tz="CET")

ggplot(sap_vpd, aes(Date.Time, Total_Sap_Flow_kg_h, color=as.factor(Tree_id)))+geom_line()+
  facet_wrap(~Treatment)+theme_bw()+guides(color="none")

sap_vpd_hourly = sap_vpd %>% filter(Total_Sap_Flow_kg_h > 0) %>% 
  mutate(datehour=floor_date(Date.Time, "1 hour")) %>% 
  group_by(datehour, Tree_id, Treatment) %>% summarize(Total_Sap_Flow=mean(Total_Sap_Flow_kg_h, na.rm=T))
sap_vpd_hourly$date = as.Date(sap_vpd_hourly$datehour, tz="CET")

#
sap_vpd_daily = sap_vpd_hourly %>% filter(Total_Sap_Flow>0) %>% group_by(date, Tree_id, Treatment) %>% 
  summarize(Total_Sap_Flow=sum(Total_Sap_Flow, na.rm=T))

ggplot(sap_vpd_daily, aes(date, Total_Sap_Flow, color=as.factor(Tree_id)))+geom_line()+
  facet_wrap(~Treatment)+theme_bw()+guides(color="none")

sap_vpd_meta = sap_vpd_daily %>% group_by(date, Treatment) %>% 
  summarize(Total_Sap_Flow=mean(Total_Sap_Flow, na.rm=T))

ggplot(sap_vpd_meta, aes(date, Total_Sap_Flow, color=as.factor(Treatment)))+geom_line()+
  theme_bw()+labs(x="",color="Treatment")

# write_csv(sap_vpd_meta, "../../Data/Pfyn/Pfyn_sap_vpd.csv")

## upscale sap flow to transpiration for VPDrought trees

sapwood_bin_size = 150 # higher sample size allows finer resolution
# re-summarize inventory data
Pfyn_inv$sapwood_area_bin = cut_width(Pfyn_inv$sapwood_area, sapwood_bin_size, boundary=0)
Pfyn_inv_VPDtrt = Pfyn_inv %>% group_by(Treatment_VPD, sapwood_area_bin) %>% summarize(n=n()) %>% mutate(prop = n/sum(n))

# first, retrieve dbh for each tree from tree_spec function
tree_specs = do.call(rbind, lapply(unique(sap_vpd$Tree_id), tree_spec))
tree_specs$meta = ifelse(grepl("irrigation",tree_specs$treat),"Irrigation","Control")

# add dummy tree for bins later
tree_specs = rbind(tree_specs, setNames(data.frame("Pinus","-999",1,"Control",0,"Control"), names(tree_specs)))

# add plot # from inventory
tree_specs$TreeNo = as.numeric(tree_specs$tree)
tree_specs = left_join(tree_specs, select(Pfyn_inv, TreeNo, Plot))
tree_specs$Plot = ifelse(is.na(tree_specs$Plot) & tree_specs$TreeNo != -999, 3, ifelse(tree_specs$Plot==0, 8, tree_specs$Plot))

# add sapwood information
tree_specs = sap_profile(rename(tree_specs, tree_dbh=dbh), sapwood_bin_size)

tree_specs = filter(tree_specs, TreeNo != -999) # remove dummy

tree_specs_ctr = filter(tree_specs, treat == "control")
tree_specs_irr = filter(tree_specs, treat == "irrigation")
tree_specs_drt = filter(tree_specs, treat == "roof")

## upscale sap flow

# normalize sap flow by sapwood area and scale with stand sapwood area

# function to scale up sap flow
UPSCALE = function(sap, trees, sap_area, plot_area) {
  
  # append sapwood area bins
  sap = left_join(sap, select(rename(trees, Tree_id=tree), Tree_id, meta, sapwood_area, sapwood_area_bin))
  
  # derive sap flux density per unit sapwood area
  sap$Sap_Flux_Density = sap$Total_Sap_Flow / sap$sapwood_area
  
  # average sap flux density for trees in each bin
  sap_bin = sap %>% group_by(datehour, date, meta, sapwood_area_bin) %>% summarize_at(vars(Sap_Flux_Density), list(mean))
  
  # append sapwood area proportion
  sap_bin = left_join(sap_bin, select(rename(Pfyn_inv_VPDtrt, meta=Treatment_VPD), sapwood_area_bin, prop))
  
  # multiply sap flux density by proportion and add trees together
  sap_scaled = sap_bin %>% group_by(datehour, date) %>% summarize(Tr = sum(Sap_Flux_Density*prop, na.rm=T))
  
  # multiply by total sapwood area and divide by plot area
  sap_scaled$Trans = sap_scaled$Tr * sap_area / plot_area
  
  # aggregate to daily
  sap_daily = sap_scaled %>% group_by(date) %>% summarize(Tr=sum(Trans))
  
  return(sap_daily)
  
}

# control
sap_vpd_ctr = filter(sap_vpd_hourly, Tree_id %in% tree_specs_ctr$tree)
sap_vpd_ctr_daily = UPSCALE(sap_vpd_ctr, tree_specs_ctr, sap_area_cont_vpd, area_cont_vpd)

ggplot(sap_vpd_ctr_daily, aes(date, Tr))+geom_line()
sum(sap_vpd_ctr_daily$Tr)

# irrigation
sap_vpd_irr = filter(sap_vpd_hourly, Tree_id %in% tree_specs_irr$tree)
sap_vpd_irr_daily = UPSCALE(sap_vpd_irr, tree_specs_irr, sap_area_irr, area_irr)

ggplot(sap_vpd_irr_daily, aes(date, Tr))+geom_line()
sum(sap_vpd_irr_daily$Tr)

# roof
sap_vpd_drt = filter(sap_vpd_hourly, Tree_id %in% tree_specs_drt$tree)
sap_vpd_drt_daily = UPSCALE(sap_vpd_ctr, tree_specs_ctr, sap_area_drt, area_drt)

ggplot(sap_vpd_drt_daily, aes(date, Tr))+geom_line()


trans_vpd_comp = rbind(mutate(sap_vpd_ctr_daily, scen="Control"), mutate(sap_vpd_irr_daily, scen="Irrigation"))
#trans_vpd_comp = rbind(trans_vpd_comp, mutate(sap_vpd_drt_daily, scen="Roof"))

ggplot(trans_vpd_comp, aes(date, Tr, color=scen))+geom_line()+theme_bw()+
  labs(x="", y="Transpiration (mm/day)", color="Treatment")

# older data

file_list = list.files("../../Data/Pfyn/data for Pfynwald/older/sap flow/", "*.RDS")

# import individual files into list
sap_list = sapply(file_list, function(x) readRDS(paste0("../../Data/Pfyn/data for Pfynwald/older/sap flow/", x)))

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
sap_daily_df = sap_df %>% mutate(datehour=floor_date(timestamp, "1 hour")) %>% 
  group_by(datehour, TreeID) %>% summarize_all(list(mean), .groups="drop") %>% 
  mutate(date=as.Date(datehour)) %>% select(-c(datehour, timestamp)) %>% 
  group_by(date, TreeID) %>% summarize_at(vars(sfd,q025,q975), list(sum))

ggplot(sap_daily_df, aes(date, sfd, color=as.factor(TreeID)))+geom_line()+
  facet_wrap(~TreeID)+theme_bw()+guides(color="none")

# write_csv(sap_daily_df, "../../Data/Pfyn/data for Pfynwald/older/sap flow/Pfyn_sap_old.csv")

# upscaling

# control scenario
sap_ctr = sap_daily_df %>% filter(TreeID %in% c("110", "124"))

# combine with measurements for weighting
sap_ctr = left_join(mutate(sap_ctr, TreeID=as.numeric(TreeID)), rename(TN_sf_ctr, TreeID=tree_name))
sap_ctr = left_join(sap_ctr, filter(Pfyn_inv_trt, Treatment=="Control"))

# upscale by sapwood area
sap_ctr_daily = sap_ctr %>% group_by(date) %>% summarize(Tr=sum(sfd/sapwood_area*rad_cor*prop))
sap_ctr_daily$Tr = sap_ctr_daily$Tr * sap_area_cont / area_cont

ggplot(sap_ctr_daily, aes(date, Tr))+geom_line()

# irrigation scenario
sap_irr = sap_daily_df %>% filter(TreeID %in% c("246", "247", "250"))

# combine with measurements for weighting
sap_irr = left_join(mutate(sap_irr, TreeID=as.numeric(TreeID)), rename(TN_sf_irr, TreeID=tree_name))
sap_irr = left_join(sap_irr, filter(Pfyn_inv_trt, Treatment=="Irrigation"))

# upscale by sapwood area
sap_irr_daily = sap_irr %>% group_by(date) %>% summarize(Tr=sum(sfd/sapwood_area*rad_cor*prop))
sap_irr_daily$Tr = sap_irr_daily$Tr * sap_area_irr / area_irr

ggplot(sap_irr_daily, aes(date, Tr))+geom_line()



# TreeNet data

TN_ctr = read_csv("../../Data/Pfyn/TreeNet/tn_timeseries_Pfyn_sap_control.csv")
TN_irr = read_csv("../../Data/Pfyn/TreeNet/tn_timeseries_Pfyn_sap_irrigation.csv")

# process raw data with TREX
library(TREX)
library(shiny)
library(plotly)

# one time series at a time
raw = TN_irr %>% filter(series_id==1018) %>% select(-series_id) %>% rename(timestamp=ts)

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
                  start="2024-01-01 00:00:00",
                  end="2024-12-31 23:50:00",
                  time.int=60,
                  max.gap=120,
                  df=TRUE)

# save(input_clip, file="input23.RData")

outlier()

input_list = readRDS("input23_Cleaned.Rds")
#input_clean = input_list$series_input
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
# write_csv(sfd_df, "Pfyn_sap_tree247.csv")