# helper functions to apply soil water potential correction in dry soils
# making tables from (Walthert, Cobos & Schleppi 2023) easily accessible

corr_coef = function(type) {
  # returns coefficients and minimum pF values from Table 2 
  # needed to apply the equation of the form pFc = a + b*pF + c*pF^2
  
  corr_coef = switch(type,
                     sand = list(a=17.219, b=-8.98, c=1.428, pFmin=3.881, psi_max=-760),
                     silt = list(a=1.328, b=-0.492, c=0.292, pFmin=3.954, psi_max=-899),
                     clay = list(a=63.963, b=-33.381, c=4.603, pFmin=3.961, psi_max=-914),
                     loam = list(a=14.508, b=-7.413, c=1.203, pFmin=3.909, psi_max=-811), # silt-sand
                     loam_humus = list(a=15.002, b=-7.825, c=1.283, pFmin=3.811, psi_max=-647), # silt-sand-humus
                     humus = list(a=5.610, b=-2.274, c=0.474, pFmin=3.755, psi_max=-569),
                     stop('Invalid soil type. Options are sand, silt, clay, loam, loam_humus, or humus.'))
  return(corr_coef)
}


soil_type = function(sand, silt, clay, OC=1) {
  # provides simplified soil type based on soil texture
  # as rough approximation of table 2
  # inputs are in %
  
  soil_type = ifelse(OC >= 10, 'humus',
              ifelse(OC >= 2, 'loam_humus',
              ifelse(clay >= 40, 'clay',
              ifelse(silt >= 60, 'silt',
              ifelse(sand >= 60, 'sand', 
                      'loam')))))
  return(soil_type)
}


pF_corr = function(pF, type) {
  # implements correction equation for dry soils based on soil texture type
  
  coef = corr_coef(type)
  
  a = coef$a
  b = coef$b
  c = coef$c
  pFmin = coef$pFmin
  
  pFc = ifelse(pF >= pFmin, a + b * pF + c * pF ^ 2, pF)
  
  SWP_corr = (10^pFc) / -10
  SWP_corr[SWP_corr < -4000] = NA
  
  return(SWP_corr)
  
}