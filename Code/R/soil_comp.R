
library(lubridate)
library(tidyverse)
library(readxl)

VPD = read_csv("../../Data/Pfyn/PFY_VPD_swp_swc.csv")
VPD$depth = VPD$depth * -1

swc = read_csv("../../Data/Pfyn/Pfyn_swat.csv")

swp = read_csv("../../Data/Pfyn/PFY_swpc_corr.csv")

swp_long = swp %>% pivot_longer(-c(meta, date), names_to="depth", values_to="SWP")
swp_long$depth = ifelse(swp_long$depth == "SWP_10cm", .1, .8)

ggplot(filter(swp_long, meta=="control"), aes(date, SWP, color="Long-term Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="control", depth<1), aes(date, SWP_corr, color="VPD Control"))+
  facet_wrap(~depth, ncol=1)+theme_bw()+labs(x="", y="SMP (kPa)", color="Source")+
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12), 
        strip.text = element_text(size=12, face="bold"),
        axis.text=element_text(size=12), axis.title=element_text(size=14))

ggplot(filter(swp_long, meta=="irrigation"), aes(date, SWP, color="Irr"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="irrigation", depth<1), aes(date, SWP_corr, color="VPD Irr"))+
  facet_wrap(~depth, ncol=1)

ggplot(filter(swc, meta=="control", date>="2022-01-01", depth %in% c(0.1,0.8)), aes(date, VWC, color="Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="control", depth<1), aes(date, VWC, color="VPD Cont"))+
  facet_wrap(~depth, ncol=1)+theme_bw()+labs(x="")

ggplot(filter(swc, meta=="irrigation", date>="2022-01-01", depth %in% c(0.1,0.8)), aes(date, VWC, color="Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="irrigation", depth<1), aes(date, VWC, color="VPD Cont"))+
  facet_wrap(~depth, ncol=1)+theme_bw()+labs(x="")

# swp comparison across scenarios
swp_long$scen = if_else(swp_long$meta=="control", "Control", 
                if_else(swp_long$meta=="irrigated", "Irrigation", "Irrigation stop"))
swp_long$depth = swp_long$depth * 100
irr = read_csv("../../Data/Pfyn/irrigation.csv")

ggplot()+geom_rect(data=filter(irr, year > 2014, year<2025), aes(x=dates, width=1, ymin=-Inf, ymax=Inf), alpha=.5, fill="lightblue")+
  geom_line(data=swp_long, aes(date, SWP, color=scen), linewidth=0.8)+
  facet_wrap(~depth, scales="free_y", ncol=1)+
  theme_bw()+labs(x="", y="SMP (kPa)", color="Scenario")+
  ggtitle("Measured Soil Water Potential Comparison across Scenarios") +
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))+
  theme(plot.title=element_text(hjust=.5), strip.text = element_text(size=12, face="bold"),
        legend.text=element_text(size=12), legend.title=element_text(size=12), 
        axis.text=element_text(size=12), axis.title=element_text(size=14))


# compare swp and sap flow

trans =  read_csv("../../Data/Pfyn/Pfyn_trans_2011_17.csv")

dates = seq.Date(as.Date("2015-01-01"), as.Date("2017-12-31"), by=1)
date_df = merge(date, unique(trans$scen))
colnames(date_df) = c("date", "scen")
trans = left_join(date_df, trans)

# control
ggplot()+geom_line(data=filter(trans, scen=="Control"), aes(date, Tr))+
  geom_line(data=filter(swp_long, scen=="Control", date<"2018-01-01"), aes(date, SWP/1000, color=as.factor(depth)))+
  theme_bw()+labs(x="", y="<--  SWP (MPa)   |   Trans (mm/day)  -->", color="Depth")+
  ggtitle("Comparison between Transpiration and Soil Water Potential for Control Scenario in Pfynwald")+
  theme(legend.position="inside", legend.position.inside=c(.1, .3), plot.title=element_text(hjust=0.5))

ggplot()+geom_line(data=filter(trans, scen=="Control"), aes(date, Tr))+
  geom_line(data=filter(swc, meta=="control", date>"2015-01-01", date<"2018-01-01"), aes(date, VWC*-10, color=as.factor(depth)))+
  theme_bw()+labs(x="", y="<-- -10*VWC   |   Trans (mm/day) -->", color="Depth")

# irrigation
ggplot()+geom_line(data=filter(trans, scen=="Irrigation", date>"2016-01-01", ), aes(date, Tr))+
  geom_line(data=filter(swp_long, scen=="Irrigation", date>"2016-01-01", date<"2018-01-01"), aes(date, SWP/1000, color=as.factor(depth)))+
  theme_bw()+labs(x="", y="<-- SWP (MPa)   |   Trans (mm/day) -->", color="Depth")

ggplot()+geom_line(data=filter(trans, scen=="Irrigation", date>"2016-01-01", ), aes(date, Tr))+
  geom_line(data=filter(swc, meta=="irrigation", date>"2016-01-01", date<"2018-01-01"), aes(date, VWC*-10, color=as.factor(depth)))+
  theme_bw()+labs(x="", y="<--  -10*VWC   |   Trans (mm/day) -->", color="Depth")

# irrigation stop
ggplot()+geom_line(data=filter(trans, scen=="Irrigation stop"), aes(date, Tr))+
  geom_line(data=filter(swp_long, scen=="Irrigation stop", date<"2018-01-01"), aes(date, SWP/1000, color=as.factor(depth)))+
  theme_bw()+labs(x="", y="<-- SWP (MPa)   |   Trans (mm/day) -->", color="Depth")+
  ggtitle("Comparison between Transpiration and Soil Water Potential for Irrigation Stop scenario in Pfynwald")+
  theme(legend.position="inside", legend.position.inside=c(.1, .3), plot.title=element_text(hjust=0.5))

ggplot()+geom_line(data=filter(trans, scen=="Irrigation stop"), aes(date, Tr))+
  geom_line(data=filter(swc, meta=="irrigation_stop", date>"2015-01-01", date<"2018-01-01"), aes(date, VWC*-10, color=as.factor(depth)))+
  theme_bw()+labs(x="", y="<-- -10*VWC   |   Trans (mm/day) -->", color="Depth")


# combine VWC records into single file

swc = read_csv("../../Data/Pfyn/PFY_swat.csv")
swc_ext = read_csv("../../Data/Pfyn/PFY_swat_ext.csv")
swc_ext = filter(swc_ext, VWC >= 0.03)
swc_ext$VWC[18834] = NA
swc_ext = na.omit(swc_ext)

swc = swc %>% select(-c(setup, ...1)) %>% filter(date<"2015-01-01")
swc$meta = if_else(swc$meta == "irrigated","irrigation",swc$meta)
swc$VWC = swc$VWC/100 * (1 - 0.62) # correct for gravel content

ggplot(filter(swc, depth %in% c(10, 80)), aes(date, VWC))+geom_line()+geom_line(data=filter(swc_ext, meta!="stop"), aes(date, VWC, color="Ext"))+
  facet_wrap(depth~meta, ncol=1)

swc_comp = rbind(swc, swc_ext)

#write_csv(swc_comp, "../../Data/Pfyn/Pfyn_swat.csv")


# compare SWP against Lorenz data

SWP20_LW = read_excel("../../Data/Pfyn/SWP_Pfynwald_LW.xlsx", sheet=1)
SWP80_LW = read_excel("../../Data/Pfyn/SWP_Pfynwald_LW.xlsx", sheet=2)

SWP20_LW$date = as.Date(SWP20_LW$`TO_CHAR(M.MEAS_DATE,'DD.MM.YYY`, format="%d.%m.%Y %H:%M")
SWP80_LW$date = as.Date(SWP80_LW$`TO_CHAR(M.MEAS_DATE,'DD.MM.YYY`, format="%d.%m.%Y %H:%M")

SWP20_LW$depth = 10
SWP80_LW$depth = 80

SWP_LW = rbind(select(SWP20_LW, date, VALUE, depth), select(SWP80_LW, date, VALUE, depth))

SWP_LW$VALUE[SWP_LW$VALUE == -9999999] = NA

declab = read_csv("N:/prj/Soil/Projekte/Pfynwald/Pfynwald_KM/1_Decent_lab_data/filtered/PFY_hh_all_filtered_2025-10-23.csv")

# hourly SWP from decent lab
dec_SWP = declab %>% filter(para == "SWP", meta=="control") %>% rename(SWP=value, datetime=hour) %>% 
  select(-c(para, sensor, meta, plot, yr, mo)) %>% mutate(date = as.Date(datetime))

dec_SWP = dec_SWP %>% group_by(date, depth) %>% summarize(SWP=mean(SWP, na.rm=T))

ggplot(filter(dec_SWP, date<="2022-10-01"), aes(date, SWP, color="DecentLab"))+geom_line(linewidth=.8)+
  geom_line(data=SWP_LW, aes(date, VALUE, color="LoWa Network"), linewidth=.8)+
  facet_wrap(~depth, ncol=1)+theme_bw()+labs(x="")+
  theme(legend.position="inside", legend.position.inside=c(.1,.2))

