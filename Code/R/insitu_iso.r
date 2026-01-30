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

ggplot(LWP, aes(Date, predawn/10, color=as.factor(tree_id)))+geom_point()+
  facet_wrap(~treatment)+guides(color="none")

ggplot(LWP, aes(predawn, midday, color=as.factor(tree_id)))+geom_point()+
  geom_abline(aes(slope=1, intercept=0))+
  facet_wrap(~treatment)+guides(color="none")

#write_csv(LWP, "../../Data/Pfyn/Insitu_lwp.csv")

# leaf water potential from VPDrought
# 2024 data
LWP24 = read_csv("../../Data/Pfyn/PFY_lwp.csv")
LWP24 = na.omit(LWP24)
LWP24$wp_value = LWP24$wp_value / 10
LWP24$date = as.Date(LWP24$MESSTIME)
LWP24$meta = case_when(LWP24$treat2=="irrigated"~"irrigation",
                          LWP24$treat2=="irrigated-VPD"~"irrigation_vpd",
                          LWP24$treat2=="roof-VPD"~"roof_vpd",
                          .default = LWP24$treat2)
LWP24 = select(LWP24, wp_value, wp, date, meta)
LWP24_pd = LWP24 %>% filter(wp=="pd") %>% group_by(date, meta) %>% 
  summarize_at(vars(wp_value), list(predawn_mean=mean, predawn_sd=sd))
LWP24_md = LWP24 %>% filter(wp=="md") %>% group_by(date, meta) %>% 
  summarize_at(vars(wp_value), list(midday_mean=mean, midday_sd=sd))
LWP24 = full_join(LWP24_pd, LWP24_md)

# 2025 data
LWP25 = read_csv("../../Data/Pfyn/PFY_lwp25.csv")
LWP25$value = LWP25$value / 10
LWP25$date = as.Date(paste(LWP25$year, LWP25$month, LWP25$day, sep="-"))
LWP25 = select(LWP25, value, wp, meta, date)
LWP25_pd = LWP25 %>% filter(wp=="pd") %>% rename(predawn=value) %>% 
  group_by(date, meta) %>% summarize_at(vars(predawn), list(predawn_mean=mean, predawn_sd=sd))
LWP25_md = LWP25 %>% filter(wp=="md") %>% rename(midday=value) %>% 
  group_by(date, meta) %>% summarize_at(vars(midday), list(midday_mean=mean, midday_sd=sd))
LWP25 = full_join(LWP25_pd, LWP25_md)
LWP25$meta = tolower(LWP25$meta)

# combine with in-situ data
LWP = read_csv("../../Data/Pfyn/Insitu_lwp.csv")
LWP$predawn = LWP$predawn / -10
LWP$midday = LWP$midday / -10
LWP = LWP %>% group_by(date, treatment) %>% summarize_at(vars(predawn, midday), list(mean=mean, sd=sd))
LWP$meta = if_else(LWP$treatment == "stop irrigation", "stop", LWP$treatment)
LWP = select(LWP, -c(treatment))
LWP = rbind(LWP, LWP24)
LWP = rbind(LWP, LWP25)

LWP_con = filter(LWP, meta=="control")#, predawn_mean < -0.4, predawn_mean > -1.1) # filter out outliers
LWP_stp = filter(LWP, meta=="stop", predawn_mean > -1.0)
LWP_irr = filter(LWP, meta=="irrigation")

# relate to TWD and VPDrought sap flow
TWD = read_csv("../../Data/Pfyn/TreeNet/Pfyn_twd_2010_25.csv")

sap = read_csv("../../Data/Pfyn/Pfyn_sap_vpd.csv")
# normalize sap flow by treatment-specific max sap flow
sap_max = sap %>% group_by(Treatment) %>% summarize(max_sap=max(Total_Sap_Flow))
sap = left_join(sap, sap_max)
sap$Sap_Flow_norm = sap$Total_Sap_Flow / sap$max_sap
sap = select(sap, -max_sap)

TWD_con = TWD %>% filter(meta=="Control") %>% 
  group_by(date) %>% summarize(twd_pdn = mean(twd_pdn, na.rm=T), twd_md=mean(twd_md, na.rm=T))
TWD_stp = TWD %>% filter(meta=="Irr_Stop") %>% 
  group_by(date) %>% summarize(twd_pdn = mean(twd_pdn, na.rm=T), twd_md=mean(twd_md, na.rm=T))
TWD_irr = TWD %>% filter(meta=="Irrigation") %>% 
  group_by(date) %>% summarize(twd_pdn = mean(twd_pdn, na.rm=T), twd_md=mean(twd_md, na.rm=T))

LWP_con = left_join(LWP_con, TWD_con)
LWP_stp = left_join(LWP_stp, TWD_stp)
LWP_irr = left_join(LWP_irr, TWD_irr)

ggplot(LWP_con, aes(twd_pdn, predawn_mean))+geom_point()+geom_smooth(method="lm")+
  labs(x="Normalized pre-dawn TWD", y="Pre-dawn LWP (MPa)")+theme_bw()
ggplot(LWP_con, aes(twd_md, midday_mean))+geom_point()+geom_smooth(method="lm")+
  labs(x="Midday TWD", y="Midday LWP (MPa)")+theme_bw()

ggplot(LWP_stp, aes(twd_pdn, predawn_mean))+geom_point()+geom_smooth(method="lm")+
  labs(x="Normalized pre-dawn TWD", y="Pre-dawn LWP (MPa)")+theme_bw()
ggplot(LWP_stp, aes(twd_md, midday_mean))+geom_point()+geom_smooth(method="lm")+
  labs(x="Midday TWD", y="Midday LWP (MPa)")+theme_bw()

ggplot(LWP_irr, aes(twd_pdn, predawn_mean))+geom_point()+geom_smooth(method="lm")+
  labs(x="Normalized pre-dawn TWD", y="Pre-dawn LWP (MPa)")+theme_bw()
ggplot(LWP_irr, aes(twd_md, midday_mean))+geom_point()+geom_smooth(method="lm")+
  labs(x="Midday TWD", y="Midday LWP (MPa)")+theme_bw()

# linear regressions
lm_pd_con = lm(predawn_mean ~ twd_pdn, LWP_con)
lm_md_con = lm(midday_mean ~ twd_md, LWP_con)

lm_pd_stp = lm(predawn_mean ~ twd_pdn, LWP_stp)
lm_md_stp = lm(midday_mean ~ twd_md, LWP_stp)

lm_pd_irr = lm(predawn_mean ~ twd_pdn, LWP_irr)
lm_md_irr = lm(midday_mean ~ twd_md, LWP_irr)

TWD_con$LWP_pd = predict(lm_pd_con, TWD_con)
TWD_con$LWP_md = predict(lm_md_con, TWD_con)

TWD_stp$LWP_pd = predict(lm_pd_stp, TWD_stp)
TWD_stp$LWP_md = predict(lm_md_stp, TWD_stp)

TWD_irr$LWP_pd = predict(lm_pd_irr, TWD_irr)
TWD_irr$LWP_md = predict(lm_md_irr, TWD_irr)

ggplot(filter(TWD_con, date>="2022-08-01", date<"2025-10-01"), aes(date, LWP_pd))+geom_line()+
  geom_pointrange(data=filter(LWP, meta=="control", date>="2022-08-01", date<"2025-10-01"), aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  theme_bw()
ggplot(filter(TWD_con, date>="2022-08-01", date<"2025-10-01"), aes(date, LWP_md))+geom_line()+
  geom_pointrange(data=filter(LWP, meta=="control", date>="2022-08-01", date<"2025-10-01"), aes(x=date, y=midday_mean, ymin=midday_mean-midday_sd, ymax=midday_mean+midday_sd), color="red", inherit.aes=F)+
  theme_bw()

ggplot(filter(TWD_stp, date>="2022-08-01", date<"2024-10-01"), aes(date, LWP_pd))+geom_line()+
  geom_pointrange(data=filter(LWP, meta=="stop", date>="2022-08-01", date<"2024-10-01"), aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  theme_bw()
ggplot(filter(TWD_stp, date>="2022-08-01", date<"2024-10-01"), aes(date, LWP_md))+geom_line()+
  geom_pointrange(data=filter(LWP, meta=="control", date>="2022-08-01", date<"2024-10-01"), aes(x=date, y=midday_mean, ymin=midday_mean-midday_sd, ymax=midday_mean+midday_sd), color="red", inherit.aes=F)+
  theme_bw()

ggplot(filter(TWD_irr, date>="2022-08-01", date<"2024-10-01"), aes(date, LWP_pd))+geom_line()+
  geom_pointrange(data=filter(LWP, meta=="irrigation", date>="2022-08-01", date<"2024-10-01"), aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  theme_bw()

TWD_lwp = rbind(mutate(TWD_con, meta="control"), mutate(TWD_stp, meta="stop"))
TWD_lwp = rbind(TWD_lwp, mutate(TWD_irr, meta="irrigation"))


# compare against soil water potential

declab = read_csv("N:/prj/Soil/Projekte/Pfynwald/Pfynwald_KM/1_Decent_lab_data/filtered/PFY_hh_all_filtered_2025-11-04.csv")
declab_test = read_csv("N:/prj/Soil/Projekte/Pfynwald/Pfynwald_KM/1_Decent_lab_data/filtered/PFY_hh_test.csv")

SWP = declab %>% filter(para == "SWP") %>% rename(SWP=value, datetime=hour) %>% 
  select(-c(para, sensor, plot, yr, mo)) %>% mutate(date = as.Date(datetime), hour=hour(datetime))
SWP = declab_test %>% filter(para == "swpt") %>% rename(SWP=value, datetime=hour) %>% 
  select(-c(para, sensor, plot, yr, mo)) %>% mutate(date = as.Date(datetime), hour=hour(datetime))

# average daily SWP
SWP_ctr = SWP %>% filter(meta=="control") %>% 
  group_by(date) %>% summarize(SWP=mean(SWP, na.rm=T))
SWP_stp = SWP %>% filter(meta=="stop") %>% 
  group_by(date) %>% summarize(SWP=mean(SWP, na.rm=T))
SWP_irr = SWP %>% filter(meta=="irrigation") %>% 
  group_by(date) %>% summarize(SWP=mean(SWP, na.rm=T))


SWP_pd = SWP %>% filter(hour > 4, hour < 7) %>% 
  group_by(date, depth, meta, node) %>% summarize(SWP=mean(SWP, na.rm=T))
SWP_md = SWP %>% filter(hour > 11, hour < 15) %>% 
  group_by(date, depth, meta, node) %>% summarize(SWP=mean(SWP, na.rm=T))


## pre-dawn
# soil water potential by depth
ggplot(filter(SWP_pd, date >= "2022-08-01"), aes(date, SWP/1000, color=as.factor(depth), fill=as.factor(depth)))+
         stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_pointrange(data=LWP, aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), inherit.aes=F)+
  facet_wrap(~meta, ncol=1)+theme_bw()+labs(x="", y="Pre-dawn SWP / LWP (MPa)", fill="Depth")+guides(color="none")+
  theme(legend.position="inside", legend.position.inside = c(.2,.5))

# average soil water potential
ggplot(filter(SWP_pd, date >= "2022-08-01"), aes(date, SWP/1000))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_pointrange(data=LWP, aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  facet_wrap(~meta, ncol=1)+theme_bw()+labs(x="", y="Pre-dawn SWP / LWP (MPa)")

# average soil water potential and twd-derived leaf water potential
ggplot(filter(SWP_pd, date >= "2022-08-01", date<"2025-10-01"), aes(date, SWP/1000))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_point(data=filter(TWD_lwp, date>="2022-08-01", date<"2025-10-01"), aes(date, LWP_pd), color="blue")+
  geom_pointrange(data=filter(LWP, date<"2025-10-01"), aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  facet_wrap(~meta, ncol=1)+theme_bw()+labs(x="", y="Pre-dawn SWP / LWP (MPa)")


# stacked plots of TWD, sap flow, and water potentials
p_sap = ggplot(filter(sap, Treatment=="control", date>="2024-05-01", date<"2025-10-01"), aes(date, Total_Sap_Flow))+geom_line()+theme_bw()+labs(x="")
p_twd = ggplot(filter(TWD, meta=="Control", date>="2024-05-01", date<"2025-10-01"), aes(date, twd_pdn))+
  stat_summary(geom="line", fun=mean, color="blue")+stat_summary(geom="ribbon", fill="blue", alpha=0.5)+
  theme_bw()+labs(x="", y="Normalized Pre-dawn TWD")#+scale_x_date(limits=c(as.Date("2024-04-01"), as.Date("2024-09-30")))
p_wp = ggplot(filter(SWP_pd, meta=="control", date>="2024-05-01", date<"2025-10-01"), aes(date, SWP/1000))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_pointrange(data=filter(LWP, meta=="control", date>="2024-05-01", date<"2025-10-01"), aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  theme_bw()+labs(x="", y="Pre-dawn SWP / LWP (MPa)")

grid.arrange(p_sap, p_twd, p_wp)

# only twd and water potentials
p_twd = ggplot(filter(TWD_con, date>="2022-08-01", date<"2025-10-01"), aes(date, twd_pdn))+geom_point()+
  theme_bw()+labs(x="")+scale_x_date(limits=c(as.Date("2022-08-01"), as.Date("2025-10-01")))
p_wp = ggplot(filter(SWP_pd, meta=="control", date>="2022-08-01", date<"2025-10-01"), aes(date, SWP/1000))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_pointrange(data=filter(LWP, meta=="control", date>="2022-08-01", date<"2025-10-01"), aes(x=date, y=predawn_mean, ymin=predawn_mean-predawn_sd, ymax=predawn_mean+predawn_sd), color="red", inherit.aes=F)+
  theme_bw()+labs(x="", y="Pre-dawn SWP / LWP (MPa)")

grid.arrange(p_twd, p_wp)

# direct comparison between pre-dawn twd and swp
# control
twd_swp_ctr = inner_join(TWD_con, SWP_ctr)
twd_swp_ctr$month = month(twd_swp_ctr$date)

ggplot(filter(twd_swp_ctr, month>4, month<10), aes(SWP, twd_pdn))+geom_point()+
  labs(x="SWP (kPa)", y="Normalized pre-dawn TWD")+theme_bw()
ggplot(filter(twd_swp_ctr, month>4, month<10), aes(SWP, twd_md))+geom_point()+
  labs(x="SWP (kPa)", y="Midday TWD")+theme_bw()

# lwp derived from twd against swp
ggplot(filter(twd_swp_ctr, date>="2022-08-01", month>4, month<10), aes(SWP/1000, LWP_pd, color="predawn"))+geom_point()+
  geom_point(aes(SWP/1000, LWP_md, color="midday"))+
  geom_abline(slope=1, intercept=0, color="black")+theme_bw()+labs(x="SWP (MPa)", y="LWP (MPa)")

# irrigation stop
twd_swp_stp = inner_join(TWD_stp, SWP_stp)
twd_swp_stp$month = month(twd_swp_stp$date)

ggplot(filter(twd_swp_stp, month>4, month<10), aes(SWP, twd_pdn))+geom_point()+
  labs(x="SWP (kPa)", y="Normalized pre-dawn TWD")+theme_bw()
ggplot(filter(twd_swp_stp, month>4, month<10), aes(SWP, twd_md))+geom_point()+
  labs(x="SWP (kPa)", y="Midday TWD")+theme_bw()

# lwp derived from twd against swp
ggplot(filter(twd_swp_stp, date>="2022-08-01", month>4, month<10), aes(SWP/1000, LWP_pd, color="predawn"))+geom_point()+
  geom_point(aes(SWP/1000, LWP_md, color="midday"))+
  geom_abline(slope=1, intercept=0, color="black")+theme_bw()+labs(x="SWP (MPa)", y="LWP (MPa)")

# irrigation
twd_swp_irr = inner_join(TWD_irr, SWP_irr)
twd_swp_irr$month = month(twd_swp_irr$date)

ggplot(filter(twd_swp_irr, month>4, month<10), aes(SWP, twd_pdn))+geom_point()+
  labs(x="SWP (kPa)", y="Normalized pre-dawn TWD")+theme_bw()
ggplot(filter(twd_swp_irr, month>4, month<10), aes(SWP, twd_md))+geom_point()+
  labs(x="SWP (kPa)", y="Midday TWD")+theme_bw()


# lwp derived from twd against swp
ggplot(filter(twd_swp_irr, date>="2022-08-01", month>4, month<10), aes(SWP/1000, LWP_pd, color="predawn"))+geom_point()+
  geom_point(aes(SWP/1000, LWP_md, color="midday"))+
  geom_abline(slope=1, intercept=0, color="black")+theme_bw()+labs(x="SWP (MPa)", y="LWP (MPa)")


## mid-day
ggplot(SWP_md, aes(date, SWP/1000, color=as.factor(depth), fill=as.factor(depth)))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_pointrange(data=LWP, aes(x=date, y=midday_mean, ymin=midday_mean-midday_sd, ymax=midday_mean+midday_sd), inherit.aes=F)+
  facet_wrap(~meta, ncol=1)+theme_bw()+labs(x="", y="Midday SWP / LWP (MPa)", fill="Depth")+guides(color="none")

ggplot(SWP_md, aes(date, SWP/1000))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=0.5)+
  geom_pointrange(data=LWP, aes(x=date, y=midday_mean, ymin=midday_mean-midday_sd, ymax=midday_mean+midday_sd), color="red", inherit.aes=F)+
  facet_wrap(~meta, ncol=1)+theme_bw()+labs(x="", y="Midday SWP / LWP (MPa)")


# sap and water potentials
sap$month = month(sap$date)

# average daily soil water potential
sap_swp = left_join(filter(sap, Treatment=="control"), SWP_ctr)

ggplot(filter(sap_swp, month>4, month<10), aes(SWP, Sap_Flow_norm))+geom_point()+
  labs(x="SWP (kPa)")+theme_bw()

# TWD-derived leaf water potential
sap_lwp = left_join(filter(sap, Treatment=="Control"), TWD_con)

ggplot(filter(sap_lwp, month>4, month<10), aes(LWP_pd, Sap_Flow_norm))+geom_point()+
  labs(x="Pre-dawn LWP (MPa)")+theme_bw()

ggplot(filter(sap_lwp, month>4, month<10), aes(LWP_md, Sap_Flow_norm))+geom_point()+
  labs(x="Midday LWP (MPa)")+theme_bw()

# measured leaf water potential

sap_lwp2 = inner_join(LWP, rename(sap, meta=Treatment))

ggplot(sap_lwp2, aes(predawn_mean, Sap_Flow_norm))+geom_point()+
  facet_wrap(~meta)+theme_bw()+labs(x="Pre-dawn LWP (MPa)")+ylim(0, 1)

ggplot(sap_lwp2, aes(midday_mean, Sap_Flow_norm))+geom_point()+
  facet_wrap(~meta)+theme_bw()+labs(x="Midday LWP (MPa)")+ylim(0, 1)



# direct comparison
SWP_pd_comp = SWP_pd %>% group_by(date, meta) %>% summarize(SWP=mean(SWP/1000, na.rm=T))
SWP_md_comp = SWP_md %>% group_by(date, meta) %>% summarize(SWP=mean(SWP/1000, na.rm=T))

SWP_pd_comp = inner_join(SWP_pd_comp, select(LWP, date, predawn_mean, meta))
SWP_md_comp = inner_join(SWP_md_comp, select(LWP, date, midday_mean, meta))

ggplot(SWP_pd_comp, aes(SWP, predawn_mean, color="Pre-dawn"))+geom_point(size=3)+
  geom_point(data=SWP_md_comp, aes(SWP, midday_mean, color="Midday"), size=3)+
  geom_abline(slope=1, intercept=0)+facet_wrap(~meta)+theme_bw()+
  labs(x="SWP (MPa)", y="LWP (MPa)")+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title=element_text(size=14),
        strip.text=element_text(size=12))


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
