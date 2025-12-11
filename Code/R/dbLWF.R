library(rLWFpg)
library(tidyverse)
library(lubridate)

con <- db_connect(username = "grauplou", password = rstudioapi::askForPassword("enter password"), db_host = "pgdblwf", db_name = "lwf")

# all data installations
install_tbl <- db_tbl(con, schema = 'ada', table = 'installation', col = c("installation_id", "installation_name"), retrieve = TRUE)

# Information about measurements
messvar.df <- db_tbl(con, schema = "ada", table = "v_messvar", retrieve = TRUE) %>%
  dplyr::filter( project_name %in% 'Pfynwald irrigation', varfreq_name != "1 Minutes")

messvar.tbl <- db_tbl(con, schema = "ada", table = "v_messvar") %>%
  dplyr::filter( project_name %in% 'Pfynwald irrigation')

pfyndata.tbl = db_tbl(con, table="pfyn_messdat")

# irrigation
irr.df = messvar.tbl %>% inner_join(pfyndata.tbl, by = 'messvar_id') %>% 
  filter(variable_id %in% c(263)) %>% 
  select(messvar_id, messtime, messval, varname_name, messvar_name, varunit_symbol, varvpos, treatment)

# precip
pcp.tbl <- ada_get_data(
  conn = con,
  messvar = c(134, 135),
  messtime_from = '2021-06-06 12:00:00',
  messtime_to = '2025-04-01 24:00:00'
)

pcp.df = pcp.tbl %>% inner_join(messvar.tbl) %>% 
  select(messvar_id, messtime, messval, varname_name, messvar_name, varunit_symbol, varvpos, treatment) %>% 
  collect()

# test difference between data collection methods
pcp.wide = pcp.df %>% pivot_wider(id_cols=messtime, names_from=messvar_name, values_from=messval)
pcp.wide$diff = pcp.wide$RaineH3_amount_dif_Tot - pcp.wide$WS100_amount_dif_Tot
ggplot(pcp.wide, aes(RaineH3_amount_dif_Tot, WS100_amount_dif_Tot))+geom_point()


# soil water potential from in-situ experiments

soilvar.tbl = messvar.tbl %>% filter(varname_name %in% c("Soil water matric potential", "Soil temperature"),
                                   vartable_name == "pfyn_messdat", varfreq_name != "1 Minutes",
                                   grepl("PFY_SS", installation_name))
soilvar.df = collect(soilvar.tbl)

data.df = pfyndata.tbl %>% 
  filter(messvar_id %in% !!soilvar.df$messvar_id) %>% 
  inner_join(soilvar.tbl, by="messvar_id") %>% 
  select(messtime, messval, messvar_id, messvar_name, varname_name, varvpos, type, treatment) %>% collect()

soil.df <- data.df %>%
  rename(    
    datetime = messtime,
    sensor = type,
    depth = varvpos,
    value = messval
  ) %>% 
  mutate(
    date = as.Date(datetime),
    para = recode(varname_name,
                  "Soil water matric potential" = "SWP",
                  "Soil temperature" = "TEM")
  ) %>%
  filter(
    !(para == "SWP" & value > 0)  # Exclude SWP values > 0
  )

soil_wide <- soil.df %>%
  select(datetime, date, treatment, depth, para, value) %>%
  group_by(datetime, date, treatment, depth, para) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = para,
    values_from = value
  )

# remove values less than -4000 kPa
soil_wide$SWP[soil_wide$SWP < -4000] = NA

# apply SWP corrections
soil_wide <- soil_wide %>%
  mutate(
    # Convert SWP from kPa to hPa, then calculate pF
    pfraw = ifelse(!is.na(SWP) & SWP < 0, log10(10 * abs(SWP)), NA),
    # Temperature deviation from 22°C
    delT = ifelse(!is.na(TEM), TEM - 22, NA),
    # Ensure delT is within reasonable limits
    delT = ifelse(delT > 10, 10, ifelse(delT < -10, -10, delT)),
    # Calculate delta pF correction (with upper bound)
    delpF = ifelse(!is.na(SWP) & !is.na(delT),
                   pmin((6.206 * delT + 0.1137 * (delT)^2) * exp(-22.76 / (pfraw + 0.1913)), 2),  # Cap at 2
                   NA),
    # Compute temperature-corrected pF
    pFtc = ifelse(!is.na(SWP) & !is.na(delpF), pfraw + delpF, NA),
    # temperature-corrected SWP
    SWP_tcorr = ifelse(!is.na(SWP) & !is.na(pFtc), pmax(((10^pFtc) * -1) / 10, -2000), NA),  # Cap at -2000 kPa
    # Correct low water potentials
    pFcr = ifelse(depth==10, ifelse(pFtc >= 3.811, 15.002 - 7.825 * pFtc + 1.283 * pFtc ^ 2, pFtc),  # soil at 10 cm uses sand-silt-humus coefficients
                  ifelse(pFtc >= 3.909, 14.508 - 7.413 * pFtc + 1.203 * pFtc ^ 2, pFtc)), # soil at >80 cm uses sand-silt coefficients
    # Convert corrected pF back to SWP (apply a lower bound)
    SWP_corr = ifelse(!is.na(SWP) & !is.na(pFcr), pmax(((10^pFcr) * -1) / 10, -4000), NA)  # Cap at -4000 kPa
  )

ggplot(soil_wide, aes(datetime, SWP_corr, color=as.factor(depth)))+geom_line()+
  facet_wrap(~treatment, ncol=1)+labs(x="",y="SWP (kPa)", color="Depth (m)")+theme_bw()

# aggregate to hourly
hourly.df = soil_wide %>% mutate(datetime=floor_date(datetime, "1 hour")) %>% select(-c(date, pfraw, delT, delpF, pFtc, pFcr)) %>% 
  group_by(datetime, treatment, depth) %>% summarize_all(list(mean))
#write_csv(hourly.df, file="../../Data/Pfyn/soil_hourly_insitu.csv")

ggplot(hourly.df, aes(datetime, SWP_corr, color=as.factor(depth)))+geom_line()+
  facet_wrap(~treatment, ncol=1)+labs(x="",y="SWP (kPa)", color="Depth (m)")+theme_bw()

# summarize daily values
daily.df <- soil_wide %>% select(-c(datetime, pfraw, delT, delpF, pFtc, pFcr)) %>% 
  group_by(date, treatment, depth) %>% summarize_all(list(mean))
#write_csv(daily.df, file="../../Data/Pfyn/soil_daily_insitu.csv")

ggplot(filter(daily.df, date>="2023-01-01"), aes(date, SWP_corr, color=as.factor(depth)))+geom_line()+
  facet_wrap(~treatment, scales="free_y", ncol=1)+
  labs(x="",y="SWP (kPa)", color="Depth (m)")+theme_bw()

ggplot(filter(daily.df, date>="2023-01-01"), aes(date, SWP_corr, color=as.factor(treatment)))+geom_line()+
  facet_wrap(~depth, scales="free_y")+
  labs(x="",y="SWP (kPa)", color="Depth (m)")+theme_bw()


# Disconnect once finished
DBI::dbDisconnect(con)


## without rLWFpg

conn = DBI::dbConnect(
  drv = RPostgres::Postgres(),
  dbname = "lwf",
  host = "pgdblwf.wsl.ch",
  port = 5432,
  user = "grauplou",
  password = rstudioapi::askForPassword("enter password"))

var_tbl = dplyr::tbl(conn, dbplyr::in_schema("ada", "v_messvar")) %>% collect()

data_tbl = dplyr::tbl(conn, dbplyr::in_schema("ada", "pfyn_messdat"))

data_df = data_tbl %>% inner_join(var_tbl) %>% filter(variable_id %in% c(263)) %>% 
  select(messvar_id, messtime, messval, varname_name, messvar_name, varunit_symbol, varvpos, treatment) %>% 
  collect()
