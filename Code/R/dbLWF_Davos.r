# retrieve Davos LWF data

library(rLWFpg)
library(readxl)
library(humidity)
library(tidyverse)
library(lubridate)

con <- db_connect(username = "grauplou", password = rstudioapi::askForPassword("enter password"), db_host = "pgdblwf", db_name = "lwf")

# all data installations
install_tbl <- db_tbl(con, schema = 'ada', table = 'installation', col = c("installation_id", "installation_name"), retrieve = TRUE)
davos_tbl = filter(install_tbl, grepl("Davos", installation_name)) # installations in davos

# Information about measurements
messvar.tbl <- db_tbl(con, schema = "ada", table = "v_messvar")
messvar.df = collect(messvar.tbl)

# all measurements in davos
davos.df = messvar.df %>% dplyr::filter(installation_name %in% davos_tbl$installation_name, vartable_name=="lwf_messdat")

# soil water measurements in davos
soilvar.tbl <- messvar.tbl %>% dplyr::filter(installation_name %in% davos_tbl$installation_name, vartable_name=="lwf_messdat",
                                             varname_name %in% c("Soil water matric potential", "Soil water content"),
                                             from_sensor != as.POSIXct("2024-09-16 13:00:00"),
                                             !(varunit_symbol %in% c("number", "ps")))
soilvar.df = collect(soilvar.tbl)

lwfdata.tbl = db_tbl(con, table="lwf_messdat")

# combine measurements with metadata
data.df = lwfdata.tbl %>% 
  filter(messvar_id %in% !!soilvar.df$messvar_id) %>% 
  inner_join(soilvar.tbl, by="messvar_id") %>% 
  select(messtime, messval, messvar_name, varname_name, varvpos) %>% collect()

# analyze and visualize soil data
soil.df <- data.df %>% rename(datetime = messtime, depth = varvpos, value = messval)
soil.df$depth = factor(soil.df$depth, levels=c(-0.05, -0.15, -0.5, -0.8, -1.5))

# soil water content
SWC = soil.df %>% filter(varname_name == "Soil water content", grepl("EC", messvar_name))
SWC$value[SWC$value < .1 | SWC$value > .52] = NA

ggplot(SWC, aes(datetime, value, color=messvar_name))+geom_line()+facet_wrap(~depth)+
  guides(color="none")

SWC = SWC %>% rename(SWC=value) %>% select(-varname_name)

SWC_hourly = SWC %>% mutate(datetime=floor_date(datetime, "1 hour")) %>% 
  group_by(datetime, depth, messvar_name) %>% summarize_all(list(mean), na.rm=T)
# write_csv(SWC_hourly, "Davos_SWC_hourly.csv")

ggplot(SWC_hourly, aes(datetime, SWC, color=messvar_name))+geom_line()+
  facet_wrap(~depth)+guides(color="none")

SWC_daily = SWC %>% mutate(date=as.Date(datetime)) %>% select(-datetime) %>% 
  group_by(date, depth, messvar_name) %>% summarize_all(list(mean), na.rm=T)
# write_csv(SWC_daily, "Davos_SWC_daily.csv")

ggplot(SWC_daily, aes(date, SWC, color=messvar_name))+geom_line()+
  facet_wrap(~depth)+guides(color="none")+theme_bw()+labs(x="")

# choose individual sensors for each depth
SWC_filt = SWC_daily %>% mutate(depth = depth * -1000) %>% 
  filter(messvar_name %in% c("EC1_15cm","EC3_50cm","EC2_80cm")) #,"EC2_150cm"))

ggplot(SWC_filt, aes(date, SWC, color=messvar_name))+geom_line()+
  theme_bw()+labs(x="")

SWC_cal = select(SWC_filt, -messvar_name)

# calibration ready output
# write_csv(SWC_cal, "../../Data/Davos/Davos_SWC_cal.csv")

# soil water potential
SWP = soil.df %>% filter(varname_name == "Soil water matric potential") %>% 
  select(-varname_name) %>% rename(SWP=value) %>% 
  filter(!is.na(depth), depth != -0.05, depth != -0.45,
         messvar_name != "PsiSoil15_kPa_Avg")
SWP$SWP[SWP$SWP < -2000 | SWP$SWP > 0] = NA

ggplot(SWP, aes(datetime, SWP, color=messvar_name))+geom_line()+facet_wrap(~depth)+
  guides(color="none")

SWP_hourly = SWP %>% mutate(datetime=floor_date(datetime, "1 hour")) %>% 
  group_by(datetime, depth, messvar_name) %>% summarize(SWP=mean(SWP, na.rm=T))
# write_csv(SWP_hourly, "Davos_SWP_hourly.csv")

ggplot(SWP_hourly, aes(datetime, SWP, color=messvar_name))+geom_line()+
  facet_wrap(~depth)+guides(color="none")

SWP_daily = SWP %>% mutate(date=as.Date(datetime)) %>% select(-datetime) %>% 
  group_by(date, depth, messvar_name) %>% summarize(SWP=mean(SWP, na.rm=T))
# write_csv(SWP_daily, "Davos_SWP_daily.csv")

ggplot(SWP_daily, aes(date, -SWP, color=messvar_name))+geom_line()+
  scale_y_log10()+facet_wrap(~depth)+guides(color="none")+theme_bw()

# average of sensors for each depth
SWP_avg = SWP_daily %>% select(-messvar_name) %>% filter(depth != -1.5) %>% 
  group_by(date, depth) %>% summarize(SWP=mean(SWP, na.rm=T)) %>% 
  mutate(depth = depth * -1000)

ggplot(SWP_avg, aes(date, SWP, color=as.factor(depth)))+geom_line()+facet_wrap(~depth)+
  guides(color="none")+theme_bw()

# calibration ready output
# write_csv(SWP_avg, "../../Data/Davos/Davos_SWP_cal.csv")


# meteo data from LWF database

metvar.tbl <- messvar.tbl %>% dplyr::filter(installation_name == "Davos Wetterdaten EMPA", vartable_name=="lwf_messdat",
                                            sensor_id %in% c(5784, 5787, 5790, 5799, 5808))
metvar.df = collect(metvar.tbl)

# combine measurements with metadata
metdata.df = lwfdata.tbl %>% 
  filter(messvar_id %in% !!metvar.df$messvar_id) %>% 
  inner_join(metvar.tbl, by="messvar_id") %>% 
  select(messtime, messval, messvar_name) %>% collect()

metwide = metdata.df %>% pivot_wider(names_from=messvar_name, values_from=messval)

metwide$year = year(metwide$messtime)
metwide = filter(metwide, year>2007, year<2025)

# derive actual vapor pressure
metwide$Es = SVP.ClaCla(metwide$TEMP+273.15) # saturation vapor pressure hPa
metwide$vappres = WVP2(metwide$FEUCHTE, metwide$Es) / 1000 # actual vapor pressure kPa

metwide$date = as.Date(metwide$messtime)
metdaily = metwide %>% group_by(date) %>% 
  summarize(globrad=mean(STRGLO)*.0864, # convert W/m2 to MJ/day/m2
            tmax=max(TEMP), tmin=min(TEMP), tmean=mean(TEMP), 
            vappres=mean(vappres), windspeed=mean(WIGE), prec=sum(REGEN))

# Disconnect once finished
DBI::dbDisconnect(con)

# compare against data from flux site

fluxdata = read_excel("../../Data/Davos/meteo/Fluxes_half_hourly.xlsx")
fluxdata$Date = as.Date(fluxdata$Date)
fluxdata$Precip = as.numeric(fluxdata$Precip)

fluxdata$ET_f[fluxdata$ET_f < 0] = 0 # ignore negative ET

# derive actual vapor pressure
fluxdata$Es = SVP.ClaCla(fluxdata$Tair_f+273.15) # saturation vapor pressure hPa
fluxdata$vappres = WVP2(fluxdata$RH, fluxdata$Es) / 1000 # actual vapor pressure kPa

fluxdaily = fluxdata %>% group_by(Date) %>% 
  summarize(ET=sum(ET_f), tmin=min(Tair_f), tmax=max(Tair_f), tmean=mean(Tair_f),
            globrad=mean(Rg_f)*.0864, vappres=mean(vappres), prec=sum(Precip))

# fill in missing precip from lwf measurements
for (i in which(is.na(fluxdaily$prec))) {
  rep_date = fluxdaily$Date[i]
  if (rep_date %in% metdaily$date)
    fluxdaily$prec[i] = metdaily$prec[which(metdaily$date == rep_date)]
  else
    fluxdaily$prec[i] = 0
}

# add windspeed from LWF measurements
fluxdaily = left_join(rename(fluxdaily, date=Date), select(metdaily, date, windspeed))
fluxdaily$windspeed[1:365] = fluxdaily$windspeed[366:730]
#fluxdaily$windspeed[is.na(fluxdaily$windspeed)] = 2

# combine lwf and flux site meteo data
meteo_daily = rbind(select(fluxdaily, -ET), filter(metdaily, date >= "2023-01-01"))

# fill in NA values
fill_na <- function(x) {
  n <- length(x)
  
  if (!anyNA(x)) return(x)
  
  i <- 1
  while (i <= n) {
    
    if (is.na(x[i])) {
      
      # Beginning of NA run
      start <- i
      
      while (i <= n && is.na(x[i])) {
        i <- i + 1
      }
      
      # End of NA run
      end <- i - 1
      run_length <- end - start + 1
      
      prev_idx <- start - 1
      next_idx <- end + 1
      
      # Only fill if both neighbours exist
      if (prev_idx >= 1 &&
          next_idx <= n &&
          !is.na(x[prev_idx]) &&
          !is.na(x[next_idx])) {
        
        # Linear interpolation
        for (k in seq_len(run_length)) {
          x[start + k - 1] <-
            x[prev_idx] +
            (x[next_idx] - x[prev_idx]) * k / (run_length + 1)
        }
      }
      
    } else {
      i <- i + 1
    }
  }
  
  x
}

# Apply to every column except date
meteo_daily[] <- lapply(meteo_daily, function(col) {
  if (is.numeric(col)) fill_na(col) else col
})

# write_csv(meteo_daily, "Davos_meteo.csv")

# also export ET as separate file
# write_csv(fluxdaily[, c("date", "ET")], "Davos_ET.csv")


# meteo comparisons
ggplot(metdaily, aes(date, globrad, color="LWF"))+geom_line()+
  geom_line(data=fluxdaily, aes(date, globrad, color="Flux site"))+
  theme_bw()+labs(x="", y="Global radiation (MJ/m2/day)", color="Source")

ggplot(metdaily, aes(date, tmax, color="LWF"))+geom_line()+
  geom_line(data=fluxdaily, aes(date, tmax, color="Flux site"))+
  theme_bw()+labs(x="", y="Max temperature (deg C)", color="Source")

ggplot(metdaily, aes(date, tmin, color="LWF"))+geom_line()+
  geom_line(data=fluxdaily, aes(date, tmin, color="Flux site"))+
  theme_bw()+labs(x="", y="Min temperature (deg C)", color="Source")

ggplot(metdaily, aes(date, vappres, color="LWF"))+geom_line()+
  geom_line(data=fluxdaily, aes(date, vappres, color="Flux site"))+
  theme_bw()+labs(x="", y="Vapor pressure (kPa)", color="Source")

ggplot(metdaily, aes(date, prec, color="LWF"))+geom_line()+
  geom_line(data=fluxdaily, aes(date, prec, color="Flux site"))+
  theme_bw()+labs(x="", y="Precip (mm)", color="Source")