# analyze in-situ experiment data from Pfynwald
# leaf water potential and isotope data

library(tidyverse)
library(lubridate)
library(readxl)

# leaf water potential

LWP = read_excel("../../Data/Pfyn/Insitu_lwp.xlsx")
LWP = LWP %>% mutate(date=as.Date(Date)) %>% 
  select(-c(Count, site, country, species, `Time pd.`, `Time md.`, `Sample Nr.`, unit, Date)) %>% 
  rename(predawn=`pre-dawn`, tree_id=`Tree Nr.`, scaffold=`Scafold Number`) %>% 
  filter(treatment %in% c("control", "irrigation", "stop irrigation")) 

ggplot(LWP, aes(Date, predawn/10, color=as.factor(tree_id)))+geom_point()+
  facet_wrap(~treatment)+guides(color="none")

ggplot(LWP, aes(predawn, midday, color=as.factor(tree_id)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+guides(color="none")

#write_csv(LWP, "../../Data/Pfyn/Insitu_lwp.csv")

# leaf water potential from VPDrought

LWP_vpd = read_csv("../../Data/Pfyn/PFY_lwp.csv")
LWP_vpd$date = as.Date(LWP_vpd$MESSTIME)
LWP_vpd = na.omit(LWP_vpd)
LWP_vpd = filter(LWP_vpd, wp=="pd")
LWP_vpd$treatment = ifelse(LWP_vpd$treat2=="irrigated","irrigation",LWP_vpd$treat2)

# compare against soil water potential

SWP = read_csv("../../Data/Pfyn/soil_hourly_VPDrought.csv")
SWP$depth = SWP$depth * -1
SWP$date = as.Date(SWP$datetime)
SWP = filter(SWP, date>="2024-04-01", date<"2025-01-01")

SWP_pd = SWP %>% mutate(hour=hour(datetime)) %>% filter(hour==6)

ggplot(filter(SWP_pd, treatment %in% c("control", "irrigation")), aes(date, SWP_corr/1000, color=as.factor(scaffold), linetype=as.factor(depth), group=interaction(as.factor(scaffold), as.factor(depth))))+geom_line()+
  geom_point(data=filter(LWP, treatment %in% c("control", "irrigation"), date>="2024-04-01", date<="2025-01-01"), aes(date, -predawn/10, color=as.factor(scaffold)), inherit.aes=F)+
  geom_point(data=filter(LWP_vpd, treatment == "irrigation"), aes(date, wp_value/10, color=as.factor(scaffold)), inherit.aes=F)+
  facet_wrap(~treatment)+theme_bw()+labs(x="", y="SWP, LWP (MPa)", color="Scaffold", linetype="Depth (m)")+
  ggtitle("Comparison between pre-dawn Leaf Water Potential and Soil Water Potential")+theme(plot.title=element_text(hjust=.5))


# isotopes

iso = read_excel("../../Data/Pfyn/Isotope_Measurements_PFY.xlsx")
iso = iso %>% rename(d18O=d18O_value, d2H=d2H_value) %>% 
  select(Identifier2, d18O, d2H) %>% filter(d18O < 0, d2H < 0)

iso_info = read_excel("../../Data/Pfyn/Sample_Information_PFY.xlsx")
iso_info$treatment = if_else(iso_info$Treatment == "stop", "irrigation stop", iso_info$Treatment)
iso_info = iso_info %>% mutate(date=as.Date(Sample_Date)) %>% 
  rename(tree_id=`Tree Nr.`, depth=`Soil Depth`) %>% 
  select(Type, Identifier2, date, treatment, tree_id, depth)

iso_meta = inner_join(iso, iso_info)

iso_soil = iso_meta %>% filter(Type=="Soil") %>% select(-c(Type, tree_id, Identifier2))
iso_xy = iso_meta %>% filter(Type=="Xylem_CVD") %>% select(-c(Type, depth, Identifier2))

ggplot(iso_soil, aes(d18O, d2H, shape=depth))+geom_point()+
  geom_point(data=iso_xy, aes(d18O, d2H, color=as.factor(tree_id)), inherit.aes=F)+
  facet_wrap(~date)+guides(color="none")

ggplot(iso_soil, aes(date, d18O, shape=depth))+geom_point()+
  geom_point(data=iso_xy, aes(date, d18O, color=as.factor(tree_id)), inherit.aes=F)+
  guides(color="none")

write_csv(iso_soil, "../../Data/Pfyn/Pfyn_insitu_soil_iso.csv")
write_csv(iso_xy, "../../Data/Pfyn/Pfyn_insitu_xylem_iso.csv")
