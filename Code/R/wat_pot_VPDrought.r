# compare above- and belowground measured water potentials
# relate to TWD

library(lubridate)
library(tidyverse)
library(readxl)

# leaf water potential in 2025 needs treatment information
LWP25 = read_excel("../../Data/Pfyn/wp_pfynwald_2025.xlsx")

# use TreeNet metadata and inventory 
TN_meta = read_csv("../../Data/Pfyn/TreeNet/tn_metadata_Pfyn_dendro_2025.csv")
TN_meta = rename(TN_meta, meta=site_subplot, tree=tree_name)
Pfyn_inv = read_excel("../../Data/Pfyn/Pfyn_DBH_2023.xlsx")
Pfyn_inv = rename(Pfyn_inv, tree=TreeNo)

LWP25 = left_join(LWP25, select(TN_meta, tree, meta))
LWP25 = left_join(LWP25, select(Pfyn_inv, tree, Treatment_VPD))

LWP25$meta = ifelse(is.na(LWP25$meta), LWP25$Treatment_VPD, LWP25$meta)
LWP25$meta = ifelse(LWP25$meta == "Buffer", "Control", LWP25$meta)
LWP25$meta = with(LWP25, case_when(tree==382 ~ "Irrigation",
                                   tree==412 ~ "Control",
                                   tree==787 ~ "Irrigation_vpd",
                                   tree==816 ~ "Irrigation_vpd",
                                   tree==837 ~ "Irrigation_vpd",
                                   tree==979 ~ "Control",
                                   .default = meta))
LWP25 = select(LWP25, -Treatment_VPD)

write_csv(LWP25, "../../Data/Pfyn/PFY_lwp25.csv")

# 2024 leaf water potential data
LWP = read_csv("../../Data/Pfyn/PFY_lwp.csv")
LWP$date = as.Date(LWP$MESSTIME)
LWP$treatment = case_when(LWP$treat2=="irrigated"~"irrigation",
                          LWP$treat2=="irrigated-VPD"~"irrigation_vpd",
                          LWP$treat2=="roof-VPD"~"roof_vpd",
                          .default = LWP$treat2)

SWP = read_csv("../../Data/Pfyn/soil_hourly_VPDrought.csv")
SWP$depth = SWP$depth * -1
SWP$date = as.Date(SWP$datetime)
SWP = filter(SWP, date>="2024-04-01", date<"2025-01-01", treatment!="control")

ggplot(LWP, aes(wp_value/10, twd_n2))+geom_point()+
  facet_wrap(~tree)

# compare SWP and pre-dawn LWP
LWP2 = na.omit(LWP)
LWP2 = filter(LWP2, wp=="pd")

SWP_pd = SWP %>% mutate(hour=hour(datetime)) %>% filter(hour==6)

ggplot(SWP, aes(date, SWP_corr/1000, linetype=as.factor(depth), group=interaction(as.factor(scaffold), as.factor(depth))))+geom_line()+
  geom_point(data=LWP2, aes(date, wp_value/10, color=treatment), inherit.aes=F)+
  facet_wrap(~treatment, scales="free_y")


ggplot(filter(SWP_pd, treatment %in% c("roof", "roof_vpd")), aes(date, SWP_corr/1000, color=as.factor(scaffold), linetype=as.factor(depth), group=interaction(as.factor(scaffold), as.factor(depth))))+geom_line()+
  geom_point(data=filter(LWP2, treatment %in% c("roof", "roof_vpd")), aes(date, wp_value/10, color=as.factor(scaffold)), inherit.aes=F)+
  facet_wrap(~treatment, scales="free_y")+theme_bw()+labs(x="", y="SWP, LWP (MPa)", color="Scaffold", linetype="Depth (m)")+
  ggtitle("Comparison between pre-dawn Leaf Water Potential and Soil Water Potential for Roof Treatments")+theme(plot.title=element_text(hjust=.5))


# compare pre-dawn against mid-day LWP
LWP3 = na.omit(LWP) %>% select(-c(twd_n2, treat2))
LWP3 = LWP3 %>% group_by(tree, wp, date, treatment) %>% summarize(lwp=min(wp_value)/10)
LWP3 = pivot_wider(LWP3, names_from=wp, values_from=lwp)

ggplot(LWP3, aes(pd, md, color=as.factor(tree)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+guides(color="none")+theme_bw()

# compare SWP & LWP

LWP4 = na.omit(LWP) %>% select(-c(twd_n2, treat2, date)) %>% rename(datetime=MESSTIME)
WP_comp = left_join(LWP4, filter(SWP, depth==.1))

ggplot(filter(WP_comp, treatment %in% c("roof", "roof_vpd")), aes(SWP_corr/1000, wp_value/10, color=as.factor(scaffold)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+theme_bw()



