
library(lubridate)
library(tidyverse)

VPD = read_csv("../../Data/Pfyn/PFY_VPD_swp_swc.csv")
VPD$depth = VPD$depth * -1

swc = read_csv("../../Data/Pfyn/Pfyn_swat.csv")
swc$VWC = swc$VWC * (1-0.62) # correct for gravel content

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

ggplot()+geom_rect(data=filter(irr, year > 2014, year<2025), aes(xmin=on, xmax=off, ymin=-Inf, ymax=Inf), alpha=.5, fill="lightblue")+
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

swc_ext = filter(swc_ext, VWC<69, VWC>10)

swc = swc %>% select(-c(setup, ...1)) %>% filter(date<"2014-04-11")
swc$meta = if_else(swc$meta == "irrigated","irrigation",swc$meta)

swc_comp = rbind(swc, swc_ext)
swc_comp$VWC = swc_comp$VWC / 100
swc_comp$depth = swc_comp$depth * 100

#write_csv(swc_comp, "../../Data/Pfyn/Pfyn_swat.csv")



