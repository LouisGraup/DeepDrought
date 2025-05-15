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
