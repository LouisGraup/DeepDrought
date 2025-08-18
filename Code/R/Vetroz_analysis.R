# combine and analyze Vetroz data

library(tidyverse)

setwd("~/DeepDrought/Data/Daten_Vetroz_Lorenz")

# TreeNet data
TN = read_csv("TreeNet_dendro/tn_timeseries.csv")
meta = read_csv("TreeNet_dendro/tn_metadata.csv")

# long-term dendro data for 2 stems
TWD12 = read_csv("TWD_hourly_2014_2022_Nodes_1_2.csv")
Stem_rad12 = read_csv("Stem_radius_hourly_2014_2022_Nodes_1_2.csv")

# soil water potential data
SWP = read_csv("SWP_20_80_110_160cm_hourly_2014_2022.csv")
SWP1cm = read_csv("SWP_1cm_Vetroz.csv")

# sap flow
sap_files = list.files(path = ".", pattern="Sapflow_node")
sap = lapply(sap_files, read_csv)
