# calculate RWU from soil water content sensors in Pfynwald

library(tidyverse)
library(lubridate)
library(zoo)

declab = read_csv("N:/prj/Soil/Projekte/Pfynwald/Pfynwald_KM/1_Decent_lab_data/filtered/PFY_hh_all_filtered_2025-10-23.csv")

# hourly VWC from decent lab
dec_VWC = declab %>% filter(para == "VWC") %>% rename(VWC=value, datetime=hour) %>% 
  filter(yr > 2014, yr < 2025) %>% select(-c(para, sensor, plot, yr, mo)) %>% 
  mutate(date = as.Date(datetime), hour=hour(datetime))

# correct for gravel
dec_VWC$VWC = dec_VWC$VWC * (1 - 0.62)

# control scenario
VWC_ctr = dec_VWC %>% filter(meta=="control") %>% select(-meta)

# node 83 has least missing values
VWC = VWC_ctr %>% filter(node == "P38") %>% select(-node)

ggplot(VWC, aes(datetime, VWC, color=as.factor(depth)))+geom_line()

# separate by depth
VWC10 = VWC %>% filter(depth == 10)
VWC80 = VWC %>% filter(depth == 80)

VWC10$dVWC = c(0, diff(VWC10$VWC))
VWC80$dVWC = c(0, diff(VWC80$VWC))

ggplot(VWC10, aes(hour, dVWC, color="10"))+geom_point()+
  geom_point(data=VWC80, aes(hour, dVWC, color="80"))

# filtered column
VWC10$dVWC_filt = VWC10$dVWC
VWC80$dVWC_filt = VWC80$dVWC

# filter for daytime values
VWC10$dVWC_filt[VWC10$hour<6 | VWC10$hour>17] = NA
VWC80$dVWC_filt[VWC80$hour<6 | VWC80$hour>17] = NA

ggplot(filter(VWC10, date>="2016-05-01", date<"2016-05-08"), aes(datetime, VWC))+geom_line()
ggplot(filter(VWC10, date>="2016-05-01", date<"2016-05-08"), aes(datetime, dVWC))+geom_line()

# ignore increases in VWC
VWC10$dVWC_filt[VWC10$dVWC_filt > 0] = NA
VWC80$dVWC_filt[VWC80$dVWC_filt > 0] = NA

# also remove large values
VWC10$dVWC_filt[VWC10$dVWC_filt < -0.001] = NA
VWC80$dVWC_filt[VWC80$dVWC_filt < -0.0004] = NA

ggplot(VWC10, aes(datetime, dVWC_filt, color="10"))+geom_point()+
  geom_point(data=VWC80, aes(datetime, dVWC_filt, color="80"))

# convert to RWU
VWC10$RWU = VWC10$dVWC_filt * 200
VWC80$RWU = VWC80$dVWC_filt * 600

ggplot(VWC10, aes(hour, RWU, color="10"))+geom_point()+
  geom_point(data=VWC80, aes(hour, RWU, color="80"))

# summarize to daily
VWC10_day = VWC10 %>% group_by(date) %>% summarize(dVWC_hr=sum(dVWC, na.rm=T), VWC=mean(VWC, na.rm=T), RWU10=sum(RWU, na.rm=T))
VWC80_day = VWC80 %>% group_by(date) %>% summarize(dVWC_hr=sum(dVWC, na.rm=T), VWC=mean(VWC, na.rm=T), RWU80=sum(RWU, na.rm=T))

VWC10_day$dVWC_d = c(0, diff(VWC10_day$VWC))
VWC80_day$dVWC_d = c(0, diff(VWC80_day$VWC))

# ignore increases in VWC
VWC10_day$dVWC_d[VWC10_day$dVWC_d > 0] = NA
VWC80_day$dVWC_d[VWC80_day$dVWC_d > 0] = NA

VWC10_day$RWU_d = VWC10_day$dVWC_d * 200
VWC80_day$RWU_d = VWC80_day$dVWC_d * 600

ggplot(VWC10_day, aes(date, RWU10, color="hour"))+geom_line()+
  geom_line(aes(date, RWU_d, color="day"))

ggplot(VWC80_day, aes(date, RWU80, color="hour"))+geom_line()+
  geom_line(aes(date, RWU_d, color="day"))

RWU = inner_join(select(VWC10_day, date, RWU10), select(VWC80_day, date, RWU80))
RWU$RWU = RWU$RWU10 + RWU$RWU80

ggplot(RWU, aes(date, -RWU))+geom_line()+theme_bw()+labs(x="", y="RWU (mm/day)")

ggplot(RWU, aes(date, RWU10, color="10"))+geom_line()+
  geom_line(aes(date, RWU80, color="80"))

RWU_long = pivot_longer(select(RWU, -RWU), -date)
RWU_long$depth = if_else(RWU_long$name == "RWU10", 10, 80)

ggplot(RWU_long, aes(date, -value, fill=as.factor(depth)))+geom_col()+theme_bw()+
  labs(x="", y="RWU (mm/day)", fill="Depth")

# function to calculate root water uptake from VWC

calc_RWU = function(id, df) {
  
  VWC = df %>% filter(node == id) %>% select(-node)
  
  # separate by depth
  VWC10 = VWC %>% filter(depth == 10)
  VWC80 = VWC %>% filter(depth == 80)
  
  # calculate delta VWC
  VWC10$dVWC = c(0, diff(VWC10$VWC))
  VWC80$dVWC = c(0, diff(VWC80$VWC))
  
  # filter for daytime values
  VWC10$dVWC[VWC10$hour<6 | VWC10$hour>17] = NA
  VWC80$dVWC[VWC80$hour<6 | VWC80$hour>17] = NA
  
  # ignore increases in VWC
  VWC10$dVWC[VWC10$dVWC > 0] = NA
  VWC80$dVWC[VWC80$dVWC > 0] = NA
  
  # also remove large values
  VWC10$dVWC[VWC10$dVWC < -0.001] = NA
  VWC80$dVWC[VWC80$dVWC < -0.0004] = NA
  
  # convert to RWU
  VWC10$RWU = VWC10$dVWC * 200
  VWC80$RWU = VWC80$dVWC * 600
  
  # summarize to daily
  VWC10_day = VWC10 %>% group_by(date) %>% summarize(RWU10=sum(RWU, na.rm=T)*-1)
  VWC80_day = VWC80 %>% group_by(date) %>% summarize(RWU80=sum(RWU, na.rm=T)*-1)
  
  RWU = inner_join(select(VWC10_day, date, RWU10), select(VWC80_day, date, RWU80))
  RWU$RWU = RWU$RWU10 + RWU$RWU80
  RWU$RWU_rm = rollmean(RWU$RWU, 14, na.pad=T)
  RWU$RWU10_rm = rollmean(RWU$RWU10, 14, na.pad=T)
  RWU$RWU80_rm = rollmean(RWU$RWU80, 14, na.pad=T)
  
  RWU$node = id
  
  return(RWU)
}

# apply function over all nodes

RWU_list = lapply(unique(dec_VWC$node), calc_RWU, dec_VWC)
RWU_all = do.call(rbind, RWU_list)

# add treatment back by node
dec_nodes = dec_VWC %>% filter(depth==10, datetime==dec_VWC$datetime[1]) %>% select(node, meta)
RWU_all = left_join(RWU_all, dec_nodes)

ggplot(RWU_all, aes(date, RWU, color=node))+geom_line()+facet_wrap(~meta, ncol=1)+
  theme_bw()+labs(x="", y="RWU (mm/day)")+guides(color="none")

# average over all nodes
RWU_avg = RWU_all %>% select(-node) %>% 
  group_by(date, meta) %>% summarize_all(list(mean))

ggplot()+geom_rect(data=filter(irr, year > 2014, year<2025), aes(xmin=on, xmax=off, ymin=-Inf, ymax=Inf), alpha=.5, fill="lightblue")+
  geom_point(data=RWU_avg, aes(date, RWU, color=meta), inherit.aes=F)+
  geom_line(data=RWU_avg, aes(date, RWU_rm, color=meta), linewidth=1.1, inherit.aes=F)+
  labs(x="", y="Root Water Uptake (mm/day)", color="Treatment")+theme_bw()+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))+
  ggtitle("Root Water Uptake from Soil Water Content Data in Pfynwald")+
  theme(legend.position="inside", legend.position.inside=c(.9,.9), plot.title=element_text(hjust=.5, size=16),
        axis.text.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text.y=element_text(size=12), 
        legend.text=element_text(size=12), legend.title=element_text(size=14))


# separate by depth
RWU_avg$dVWC10_rm = RWU_avg$RWU10_rm / 200
RWU_avg$dVWC80_rm = RWU_avg$RWU80_rm / 600
dVWC_long = pivot_longer(select(RWU_avg, date, meta, dVWC10_rm, dVWC80_rm), -c(date, meta))
dVWC_long$depth = if_else(dVWC_long$name == "dVWC10_rm", 10, 80)

# average over all years
dVWC_long$month = month(dVWC_long$date)
dVWC_long$day = day(dVWC_long$date)

dVWC_yr = dVWC_long %>% select(-c(date, name)) %>% 
  group_by(meta, month, day, depth) %>% summarize(dVWC=mean(value))
dVWC_yr$date = as.Date(paste(2000, dVWC_yr$month, dVWC_yr$day, sep="-"))

ggplot(dVWC_yr, aes(date, dVWC, color=as.factor(depth)))+geom_line(linewidth=1.5)+
  scale_x_date(date_breaks="1 month", date_labels="%b")+
  facet_wrap(~meta, ncol=1)+theme_bw()+labs(x="", y="dVWC", color="Depth")


# compare RWU to transpiration derived from sap flow

trans =  read_csv("../../Data/Pfyn/Pfyn_trans_2011_17.csv")

RWU_comp = RWU_avg %>% select(date, meta, RWU, RWU_rm)
RWU_comp$scen = if_else(RWU_comp$meta=="control", "Control", 
                if_else(RWU_comp$meta=="stop", "Irrigation stop", "Irrigation"))
RWU_comp = select(RWU_comp, -meta)

# combine with overlapping years
date_df = data.frame(date=seq.Date(as.Date("2015-01-01"), as.Date("2017-12-31"), by=1))
RWU_comp = left_join(date_df, RWU_comp)
RWU_comp = left_join(RWU_comp, trans)
RWU_comp$RWU[is.na(RWU_comp$Tr)] = NA

ggplot(RWU_comp, aes(date, RWU, color="RWU"))+geom_line()+
  geom_line(aes(date, Tr, color="Trans"))+facet_wrap(~scen, ncol=1)+
  theme_bw()+labs(x="", y="Water Uptake (mm/day)", color="Source")+
  theme(legend.position="inside", legend.position.inside=c(.2, .5))
