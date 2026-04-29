
## use helper functions to create input files for LWFBrook90.jl

library(tidyverse)
library(lubridate)
library(LWFBrook90R)

setwd("../../Data/Daten_Vetroz_Lorenz")

## read in meteo data

meteo = read.csv("meteo.csv")
meteo$dates = as.Date(meteo$date)
meteo = rename(meteo, prec=prec_Sion)

## soil and root data

# derive soil hydraulic parameters from soil texture data
soil = read_csv("soil_texture.csv")
soil$c_org = soil$c_org / 10 # convert to kg/kg
soil$gravel = soil$gravel / 100 # convert to %
soil$upper = soil$upper / 100 # convert to m
soil$thickness = soil$thickness / 100 # convert to m
soil = mutate(soil, lower=upper - thickness, .before=thickness)

soil_hyd = cbind(soil, hydpar_wessolek_tab(texture=soil$texture))
soil_ext = select(soil_hyd, !c(upper, thickness))
#soil_hyd2 = cbind(soil, hydpar_puh2(clay=soil$clay, silt=soil$silt, sand=soil$sand, bd=soil$bd, oc.pct=soil$c_org))

# Mualem van Genuchten - derived variables
MvG_SWC = function(psi, alpha, n, ths, thr) {
  # psi in hPa
  m = 1 - 1/n
  
  wetness = 1/((1+(alpha*psi)^n)^m)
  
  theta = wetness * (ths - thr) + thr
  
  return(theta)
}

soil_hyd$sat_vol = with(soil_hyd, MvG_SWC(1, alpha/100, npar, ths, thr))
soil_hyd$fc_vol = with(soil_hyd, MvG_SWC(63, alpha/100, npar, ths, thr))
soil_hyd$pwp_vol = with(soil_hyd, MvG_SWC(15848, alpha/100, npar, ths, thr))
soil_hyd$res_vol = with(soil_hyd, MvG_SWC(10000000, alpha/100, npar, ths, thr))
soil_hyd$awc_vol = soil_hyd$fc_vol - soil_hyd$pwp_vol
soil_hyd$awc_mm = soil_hyd$awc_vol * (1 - soil_hyd$gravel) * soil_hyd$thickness * 1000
sum(soil_hyd$awc_mm)

# analyze fine root fraction
rootden_per_fineearthvol = c(4.736, 3.169, 1.59, 1.569, 0.68, 1.205)
rootfrac_per_fineearthvol = rootden_per_fineearthvol / sum(rootden_per_fineearthvol)

# actual depths measured
#root_upper = c(0, 15, 45, 75, 115, 145)
#root_lower = c(10, 25, 55, 85, 125, 155)
root_upper = c(0, 15, 35, 65, 100, 140) / -100
root_lower = c(15, 35, 65, 100, 140, 165) / -100

root_df = data.frame(upper = root_upper, lower=root_lower, rootden = rootden_per_fineearthvol)

# extend soil with sensor nodes
layer_nodes = c(-0.19, -0.79, -1.09, -1.1, -1.59, -1.6) # depths of additional layers needed
row_ids = c(2, 5, 6, 6, 7, 7) # ids correspond to rows to be copied
for (i in 1:length(layer_nodes))  {
  
  d = layer_nodes[i] # depth to add
  r = row_ids[i] # row to copy
  
  row = soil_ext[r,]
  row$lower = d
 
  soil_ext = rbind(soil_ext, row) 
  
}
soil_ext = sort_by(soil_ext, soil_ext$lower, decreasing=T)
soil_ext = mutate(soil_ext, upper = c(0,soil_ext$lower[1:(nrow(soil_ext)-1)]), .before=lower)

# combine soil with root data
soil_nodes = c(max(soil_ext$upper), soil_ext$lower)
root_table = make_rootden(soilnodes = soil_nodes, method="table", rootdat=root_df)
soil_ext$rootden = root_table

## use LWFBrook90R to prepare input files for julia version

source("../../../LWFBrook90.jl/src/generate_LWFBrook90jl_Input_mod.R", echo=F)

opt = set_optionsLWFB90(startdate=as.Date("2010-01-01"), enddate=as.Date("2024-12-31"), 
                        root_method="soilvar", budburst_method="Menzel", 
                        leaffall_method="vonWilpert", lai_method="Coupmodel")

# coordinates 46.222265 N, 7.254621 E
par = set_paramLWFB90(maxlai=4.12, winlaifrac=0, height=12, height_ini=12, age_ini=82, bypar=1,
                      coords_y=46.222265, eslope=20, aspect=139, budburst_species="Quercus robur")
generate_LWFBrook90jl_Input("Vetroz","vetroz",".", options_b90=opt, param_b90=par, climate=meteo, soil=soil_ext)

