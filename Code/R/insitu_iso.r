# analyze in-situ experiment data from Pfynwald
# leaf water potential and isotope data

library(tidyverse)
library(lubridate)
library(gridExtra)
library(readxl)

# leaf water potential from in-situ

LWP = read_excel("../../Data/Pfyn/Insitu_lwp.xlsx")
LWP = LWP %>% mutate(date=as.Date(Date)) %>% 
  select(-c(Count, site, country, species, `Time pd.`, `Time md.`, `Sample Nr.`, unit, Date)) %>% 
  rename(predawn=`pre-dawn`, tree_id=`Tree Nr.`, scaffold=`Scafold Number`) %>% 
  filter(treatment %in% c("control", "irrigation", "stop irrigation")) 

ggplot(LWP, aes(date, predawn/10, color=as.factor(tree_id)))+geom_point()+
  facet_wrap(~treatment)+guides(color="none")

ggplot(LWP, aes(predawn, midday, color=as.factor(tree_id)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+guides(color="none")

#write_csv(LWP, "../../Data/Pfyn/Insitu_lwp.csv")


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

# combine dates
iso_meta$date = with(iso_meta, case_when(date=="2023-08-14" ~ as.Date("2023-08-15"),
                                         date=="2023-09-21" ~ as.Date("2023-09-22"),
                                         date=="2024-05-15" ~ as.Date("2024-05-16"),
                                         date=="2024-07-05" ~ as.Date("2024-07-06"),
                                         date=="2024-08-15" ~ as.Date("2024-08-16"),
                                         date=="2024-10-26" ~ as.Date("2024-10-27"),
                                         .default = date))

iso_soil = iso_meta %>% filter(Type=="Soil") %>% select(-c(Type, tree_id, Identifier2))
iso_xy = iso_meta %>% filter(Type=="Xylem_CVD") %>% select(-c(Type, depth, Identifier2))

iso_soil$depth = factor(iso_soil$depth, levels=c("30-50","10-30","0-10"))
iso_soil2 = iso_soil %>% group_by(date, depth, treatment) %>% summarize_at(vars(d18O, d2H), list(mean=mean, sd=sd))
iso_xy2 = iso_xy %>% group_by(date, treatment) %>% summarize_at(vars(d18O, d2H), list(mean=mean, sd=sd))

ggplot(iso_soil, aes(d18O, d2H, shape=depth))+geom_point()+
  geom_point(data=iso_xy, aes(d18O, d2H, color=as.factor(tree_id)), inherit.aes=F)+
  facet_wrap(~date)+guides(color="none")

ggplot(iso_soil, aes(date, d18O, shape=depth))+geom_point()+
  geom_point(data=iso_xy, aes(date, d18O, color=as.factor(tree_id)), inherit.aes=F)+
  guides(color="none")

ggplot(iso_soil, aes(date, d18O, shape=depth, color=treatment))+geom_point()+
  facet_grid(depth~treatment)+theme_bw()

ggplot(iso_xy, aes(date, d18O, color=as.factor(tree_id)))+geom_point()+
  facet_wrap(~treatment, ncol=1)+theme_bw()+guides(color="none")

ggplot(iso_soil, aes(d18O, depth, color=treatment))+stat_summary(geom="point", fun=mean)+
  geom_vline(data=iso_xy2, aes(xintercept=d18O_mean, color=treatment))+facet_wrap(~date)

ggplot(iso_soil2, aes(d18O_mean, depth, color=treatment))+
  geom_pointrange(aes(xmin=d18O_mean-d18O_sd, xmax=d18O_mean+d18O_sd))+
  geom_vline(data=iso_xy2, aes(xintercept=d18O_mean, color=treatment))+
  geom_rect(data=iso_xy2, aes(xmin=d18O_mean-d18O_sd, xmax=d18O_mean+d18O_sd, ymin=-Inf, ymax=Inf, fill=treatment), alpha=0.3, inherit.aes=F)+
  facet_wrap(~date)+theme_bw()+labs(x="d18O", y="Depth (cm)")

iso_soil_out = iso_soil %>% group_by(date, depth, treatment) %>% summarize_at(vars(d18O, d2H), list(mean))
iso_xy_out = iso_xy %>% group_by(date, treatment) %>% summarize_at(vars(d18O, d2H), list(mean))

#write_csv(iso_soil_out, "../../Data/Pfyn/Pfyn_insitu_soil_iso.csv")
#write_csv(iso_xy_out, "../../Data/Pfyn/Pfyn_insitu_xylem_iso.csv")
