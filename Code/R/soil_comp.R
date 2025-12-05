
library(lubridate)
library(tidyverse)

VPD = read_csv("../../Data/Pfyn/PFY_VPD_swp_swc.csv")
VPD$depth = VPD$depth * -1

swc = read_csv("../../Data/Pfyn/PFY_swat.csv")
swc$depth = swc$depth / 100

swc_ext = read_csv("../../Data/Pfyn/PFY_swat_ext.csv")
swc_ext$depth = swc_ext$depth / 100
swc_ext$VWC = swc_ext$VWC * 100

swp = read_csv("../../Data/Pfyn/PFY_swpc.csv")

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

ggplot(filter(swc, meta=="control", depth %in% c(0.1,0.8)), aes(date, VWC, color="Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="control", depth<1), aes(date, VWC*100, color="VPD Cont"))+
  geom_line(data=filter(swc_ext, meta=="control"), aes(date, VWC, color="Control Ext"))+
  facet_wrap(~depth, ncol=1)+theme_bw()+labs(x="")

ggplot(filter(swc, meta=="irrigated"), aes(date, VWC, color="Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="irrigation", depth<1), aes(date, VWC*100, color="VPD Cont"))+
  geom_line(data=filter(swc_ext, meta=="irrigation", VWC>10, VWC<50), aes(date, VWC, color="Control Ext"))+
  facet_wrap(~depth, ncol=1)


ggplot(swp_long, aes(date, SWP, color=meta))+geom_line()+
  facet_wrap(~depth, scales="free_y", ncol=1)+
  theme_bw()+labs(x="", y="SMP (kPa)", color="Scenario")+
  ggtitle("Measured Soil Water Potential Comparison acros Scenarios") +
  theme(plot.title=element_text(hjust=.5), strip.text = element_text(size=12, face="bold"),
        legend.text=element_text(size=12), legend.title=element_text(size=12), 
        axis.text=element_text(size=12), axis.title=element_text(size=14))

# combine VWC records into single file

swc_ext = filter(swc_ext, VWC<69, VWC>10)

swc = swc %>% select(-c(setup, ...1)) %>% filter(date<"2014-04-11")
swc$meta = if_else(swc$meta == "irrigated","irrigation",swc$meta)

swc_comp = rbind(swc, swc_ext)
swc_comp$VWC = swc_comp$VWC / 100
swc_comp$depth = swc_comp$depth * 100

#write_csv(swc_comp, "../../Data/Pfyn/Pfyn_swat.csv")
