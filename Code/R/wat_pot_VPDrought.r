# compare above- and belowground measured water potentials
# relate to TWD

library(lubridate)
library(tidyverse)

LWP = read_csv("../../Data/Pfyn/PFY_lwp.csv")
LWP$date = as.Date(LWP$MESSTIME)
LWP$treatment = case_when(LWP$treat2=="irrigated"~"irrigation",
                          LWP$treat2=="irrigated-VPD"~"irrigation_vpd",
                          LWP$treat2=="roof-VPD"~"roof_vpd",
                          .default = LWP$treat2)

SWP = read_csv("../../Data/Pfyn/PFY_VPD_swp_swc.csv")
SWP$depth = SWP$depth * -1

SWP = filter(SWP, date>="2024-01-01", date<"2025-01-01", treatment!="control")

ggplot(LWP, aes(wp_value/10, twd_n2))+geom_point()+
  facet_wrap(~tree)

# compare daily SWP and pre-dawn LWP
LWP2 = na.omit(LWP)
LWP2 = filter(LWP2, wp=="pd")

ggplot(SWP, aes(date, SWP_corr/1000, linetype=as.factor(depth)))+geom_line()+
  geom_point(data=LWP2, aes(date, wp_value/10, color=treatment), inherit.aes=F)+
  facet_wrap(~treatment, scales="free_y")

# compare pre-dawn against mid-day LWP
LWP3 = na.omit(LWP) %>% select(-c(twd_n2, treat2))
LWP3 = LWP3 %>% group_by(tree, wp, date, treatment) %>% summarize(lwp=min(wp_value)/10)
LWP3 = pivot_wider(LWP3, names_from=wp, values_from=lwp)

ggplot(LWP3, aes(pd, md, color=as.factor(tree)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+guides(color="none")+theme_bw()


WP_comp = left_join(LWP3, filter(SWP, depth==.1))

ggplot(WP_comp, aes(SWP_corr/1000, pd))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+theme_bw()
