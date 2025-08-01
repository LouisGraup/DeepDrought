## retrieve and organize Pfynwald soil water data for VPDrought experiment

library(lubridate)
library(tidyverse)
library(rLWFpg)
library(humidity)


start_date <- as.POSIXct("2023-08-01 00:00:00", tz="GMT")
end_date <- as.POSIXct("2025-07-01 00:00:00", tz="GMT")

conn <- db_connect('grauplou', rstudioapi::askForPassword("Database password"), db_host = 'pgdblwf', db_name = 'lwf')


messvar.df = db_tbl(conn, 'ada', 'v_messvar', retrieve = T) %>%
  filter(project_name %in% "Pfynwald irrigation")

soilvar.tbl = db_tbl(conn, 'ada', 'v_messvar', retrieve = F) %>%
  filter(project_name == "Pfynwald irrigation", vartable_name == "pfyn_messdat",
         varname_name %in% c("Soil water content", "Soil water matric potential", "Soil temperature", "Soil electrical conductivity"),
         varfreq_name == "10 Minutes", grepl("PFY_SC", installation_name))
soilvar.df = collect(soilvar.tbl)

pfyndata.tbl = db_tbl(conn, table="pfyn_messdat", retrieve = F)

# code from Katrin

data.df = pfyndata.tbl %>% 
  filter(messvar_id %in% soilvar.df$messvar_id,
         messtime >= start_date, messtime < end_date) %>% 
  inner_join(soilvar.tbl, by="messvar_id") %>% 
  select(messtime, messval, messvar_id, messvar_name, varname_name, varvpos, type, installation_name, scaffold, treatment, descr) %>% collect()

VPD.df <- data.df %>%
  rename(    
    datetime = messtime,
    site = installation_name,
    sensor = type,
    depth = varvpos,
    value = messval
  ) %>% 
  mutate(
    date = as.Date(datetime),
    para = recode(varname_name,
                  "Soil water content" = "VWC",
                  "Soil water matric potential" = "SWP",
                  "Soil electrical conductivity" = "ECS",
                  "Soil temperature" = "TEM")
  ) %>%
  mutate(
    value = ifelse(para == "SWP" & value == 0, -0.1, value)
  ) %>%
  filter(
    depth != "-0.01", sensor != "MPS 6",
    !(para == "SWP" & value > 0)  # Exclude SWP values > 0
  )

write_csv(VPD.df, file="../../Data/Pfyn/soil_VPDrought.csv")

## need to correct for temperature at finer measurement resolution (sub-daily)

# Convert to wide format with separate columns for SWP and TEM
VPD_wide <- VPD.df %>%
  select(datetime, date, site, treatment, depth, para, descr, value) %>%
  group_by(datetime, date, site, treatment, depth, para, descr) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = para,
    values_from = value
  )
summary(VPD_wide)

# Apply temperature correction
VPD_wide <- VPD_wide %>%
  mutate(
    # Convert SWP from kPa to hPa, then calculate pF
    pfraw = ifelse(!is.na(SWP) & SWP < 0, log10(10 * abs(SWP)), NA),
    # Temperature deviation from 22°C
    delT = ifelse(!is.na(TEM), TEM - 22, NA),
    # Ensure delT is within reasonable limits
    delT = ifelse(delT > 10, 10, ifelse(delT < -10, -10, delT)),
    # Calculate delta pF correction (with upper bound)
    delpF = ifelse(
      !is.na(SWP) & !is.na(delT),
      pmin((6.206 * delT + 0.1137 * (delT)^2) * exp(-22.76 / (pfraw + 0.1913)), 2),  # Cap at 2
      NA),
    # Compute corrected pF
    pFcr = ifelse(!is.na(SWP) & !is.na(delpF), pfraw + delpF, NA),
    # Convert corrected pF back to SWP (apply a lower bound)
    SWP_corr = ifelse(!is.na(SWP) & !is.na(pFcr), pmax(((10^pFcr) * -1) / 10, -3000), NA)  # Cap at -1500 kPa
  )
head(VPD_wide)

# summarize daily values
daily.df <- VPD_wide %>% select(-c(datetime, pfraw, delT, delpF, pFcr)) %>% 
  group_by(date, site, treatment, depth, descr) %>%
  summarize_all(list(mean))

#write_csv(daily.df, file="../../Data/Pfyn/soil_daily_VPDrought.csv")

# summarize across treatments
#scen.df = daily.df %>% group_by(date, treatment, depth) %>% summarize_at(vars(SWP_corr, VWC), list(mean))
#write_csv(scen.df, file="../../Data/Pfyn/PFY_VPD_swp_swc.csv")

## start here
daily.df = read_csv("../../Data/Pfyn/soil_daily_VPDrought.csv")

daily.df$year = year(daily.df$date)
daily.df$depth = factor(daily.df$depth, levels=c(-0.01, -0.10, -0.80, -1.20), labels=c("1 cm", "10 cm", "80 cm", "120 cm"))
daily.df$treatment = factor(daily.df$treatment, 
                            levels=c("control","roof","roof_vpd","irrigation","irrigation_vpd"))


# plot mean and standard deviation per treatment and depth
ggplot(daily.df, aes(date, VWC, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~depth, ncol=1)+theme_bw()

ggplot(daily.df, aes(date, ECS, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~depth, ncol=1)+theme_bw()

ggplot(filter(daily.df, date>="2024-04-01"), aes(date, SWP_corr/1000, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~depth, scales="free_y", ncol=1)+theme_bw()+labs(x="", y="SMP (MPa)")+
  theme(legend.position="bottom", legend.text=element_text(size=12),
        legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12), axis.title=element_text(size=14))+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

ggplot(filter(daily.df, date>="2025-03-01"), aes(date, SWP_corr/1000, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~depth, ncol=1, scales="free")+theme_bw()+labs(x="", y="SMP (MPa)")+
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12), strip.text=element_text(size=12, face="bold"),
        axis.text=element_text(size=12), axis.title=element_text(size=14))


# plot mean and standard deviation per depth for roof treatments
ggplot(filter(daily.df, date>="2024-04-01", treatment %in% c("roof_vpd", "roof")), aes(date, SWP_corr/1000, color=depth, fill=depth))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~treatment, ncol=1)+theme_bw()+labs(x="", y="SMP (MPa)")


ggplot(filter(daily.df, date>="2024-04-01", treatment %in% c("roof_vpd", "roof")), aes(date, SWP_corr/1000, color=treatment, fill=treatment))+
  stat_summary(geom="line", fun=mean)+stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~depth, scales="free_y", ncol=1)+theme_bw()+labs(x="", y="SMP (MPa)")+
  theme(legend.position=c(.1,.85), legend.text=element_text(size=12),
        legend.title=element_text(size=12), strip.text = element_text(size=12, face="bold"),
        axis.text=element_text(size=12), axis.title=element_text(size=14))


ggplot(filter(daily.df, treatment %in% c("roof_vpd", "roof"), year==2025), aes(date, SWP_corr/1000, color=depth, fill=depth))+
  stat_summary(geom="line", fun=mean)+
  stat_summary(geom="ribbon", alpha=.3)+
  facet_wrap(~treatment)+theme_bw()+labs(x="", y="SMP (MPa)")


# plot roof treatments for all sensors
ggplot(filter(daily.df, !is.na(descr)), aes(date, SWP_corr/1000, group=interaction(depth, site)))+geom_line(aes(linetype=depth, color=as.factor(site)))+
  facet_grid(descr~treatment)+theme_bw()+labs(x="", y="SMP (MPa)")

# by depth, treatment and position
ggplot(filter(daily.df, !is.na(descr)), aes(date, SWP_corr/1000, group=interaction(descr, site)))+geom_line(aes(color=as.factor(descr)))+
  facet_grid(depth~treatment, scales="free_y")+theme_bw()+labs(x="", y="SMP (MPa)")


# plot roof treatments from recent sensors only
ggplot(filter(daily.df, date>="2025-04-25", !is.na(descr)), aes(date, SWP_corr/1000, group=interaction(depth, site)))+
  geom_line(aes(linetype=depth, color=as.factor(site)))+
  facet_grid(descr~treatment)+theme_bw()+labs(x="", y="SMP (MPa)")

# across plots by position relative to panel
ggplot(filter(daily.df, date>="2025-04-25", !is.na(descr), depth=="10 cm"), aes(date, SWP_corr/1000, group=as.factor(descr)))+
  geom_line(aes(color=as.factor(descr)))+
  facet_wrap(~site, ncol=2)+theme_bw()+labs(x="", y="SMP (MPa)")+
  ggtitle("Comparison of SMP sensors at 10 cm for roof and roof_vpd treatments")
