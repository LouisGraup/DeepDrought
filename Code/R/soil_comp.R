
library(lubridate)
library(tidyverse)

VPD = read_csv("../../Data/Pfyn/PFY_VPD_swp_swc.csv")
VPD$depth = VPD$depth * -1

swc = read_csv("../../Data/Pfyn/PFY_swat.csv")
swc$depth = swc$depth / 100

swp = read_csv("../../Data/Pfyn/PFY_swpc.csv")

swp_long = swp %>% pivot_longer(-c(meta, date), names_to="depth", values_to="SWP")
swp_long$depth = ifelse(swp_long$depth == "SWP_10cm", .1, .8)

ggplot(filter(swp_long, meta=="control"), aes(date, SWP, color="Long-term Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="control", depth<1), aes(date, SWP_corr, color="VPD Control"))+
  facet_wrap(~depth, ncol=1)+theme_bw()+labs(x="", y="SMP (kPa)", color="Source")+
  theme(legend.position="none")

ggplot(filter(swp_long, meta=="irrigated"), aes(date, SWP, color="Irr"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="irrigation", depth<1), aes(date, SWP_corr, color="VPD Irr"))+
  facet_wrap(~depth, ncol=1)

ggplot(filter(swc, meta=="control"), aes(date, VWC, color="Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="control", depth<1), aes(date, VWC*100, color="VPD Cont"))+
  facet_wrap(~depth, ncol=1)+theme_bw()+theme(legend.position="none")+labs(x="")

ggplot(filter(swc, meta=="irrigated"), aes(date, VWC, color="Control"))+geom_line()+
  geom_line(data=filter(VPD, treatment=="irrigation", depth<1), aes(date, VWC*100, color="VPD Cont"))+
  facet_wrap(~depth, ncol=1)



