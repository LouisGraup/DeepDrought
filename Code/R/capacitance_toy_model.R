# toy model of plant water storage capacitance

library(tidyverse)

# Mualem van Genuchten functions
SWP_MvG = function(SWC) {
  ths = .42 * (1 - 0.62)
  alpha = 5
  n = 1.23
  m = 1 - 1 / n
  
  # calculates SWP in kPa
  SWP = 9.81 * (-1 / alpha) * ((SWC / ths) ^ (-1 / m) - 1) ^ (1 / n)
  return(SWP) 
}

VWC_MvG = function(SWP) {
  # requires SWP in kPa
  ths = .42 * (1 - 0.62)
  alpha = 5
  n = 1.23
  m = 1 - 1 / n
  
  VWC = ths/((1 + (-alpha * SWP / 9.81) ^ n) ^ m)
  return(VWC)
}

# model parameters
n = 50 # time steps

# soil parameters
SWATMAX = 100 # total soil water storage (mm)
TH_S = 0.42 * 0.38 # field capacity * gravel content (%)
dz = SWATMAX / TH_S # effective soil depth (mm)
PSI0 = -0.1 # initial soil and plant potential

# plant parameters
C = 2 # plant storage capacitance (mm/MPa)
Vc_max = 20 # maximum storage capacity (mm)
Ks = 2.4 # plant storage conductance (mm/day/MPa)
Kp = 7 # total plant conductance (mm/day/MPa)
PSI_CR = -1.2 # MPa

# derived parameters
fX = 0.5 # xylem fraction (-) (only used to split Kp into above- and below-ground conductances)
Kx = Kp * fX
Kr = Kp * (1 - fX)
FC = VWC_MvG(PSI0 * 1000) # "field capacity" or initial water content
WP = VWC_MvG(PSI_CR * 1000) # wilting point
AWC = FC - WP
PAWSC = AWC * dz

# forcing
PET = rep(1, n) # potential ET (mm/day)

# states
PL_ST_pd = rep(0, n) # plant storage volume at predawn (mm)
PL_ST_md = rep(0, n) # plant storage volume at midday (mm)
PL_WP_pd = rep(0, n) # plant storage potential at predawn (MPa)
PL_WP_md = rep(0, n) # plant storage potential at midday (MPa)
SO_WP = rep(0, n) # soil storage potential (MPa)
SO_WC = rep(0, n) # soil volumetric water content (%)
SWAT  = rep(0, n) # soil water storage (mm)

# initialize states
PL_WP_pd[1] = PSI0
PL_WP_md[1] = PSI0
PL_ST_pd[1] = Vc_max + C * PSI0
PL_ST_md[1] = Vc_max + C * PSI0
SO_WP[1] = PSI0
SO_WC[1] = FC
SWAT[1]  = SWATMAX * (FC / TH_S)

# fluxes
RWU = rep(0, n) # root water uptake (mm/day)
PLFL = rep(0, n) # plant flow (mm/day)
PLRF = rep(0, n) # plant storage refilling (mm/day)

for (i in 2:n) {
  
  ## solve the series of resistances
  
  KT = Kr + Ks # total conductance
  
  # weighted total water potential
  PSIT = PL_WP_pd[i-1]*(Ks/KT) + SO_WP[i-1]*(Kr/KT)
  
  # soil and plant water available for transpiration
  SUPPLY = (PSIT - PSI_CR) * (KT + Kx)
  
  # actual ET
  AET = min(SUPPLY, PET[i])
  
  # distribute uptake according to proportion of conductance in each store
  #RWU[i] = (SO_WP[i-1] - PSIT + AET / KT) * Kr
  #PLFL[i] = (PL_WP[i-1] - PSIT + AET / KT) * Ks
  
  ## from ChatGPT to resolve sub-daily behavior without sub-daily timesteps
  
  D = 0.5 # daylight fraction
  omega = pi / D
  
  # calculate maximum rate of transpiration within a day
  # assuming sinusoidal forcing (aligned with LWFBrook90)
  # of the form T(t) = Tmax * sin(pi*t/D)
  Tmax = (omega/2) * AET
  
  # analytical solution to ODE of the form
  # C * dPSI/dt = Kr * (PSI_s - PSI_p) - T(t)
  # for calculating midday water potential
  psi_midday = SO_WP[i-1] - Tmax / sqrt(Kr^2 + (omega*C)^2)
  
  # storage release
  PLFL[i] = C * (PL_WP_pd[i-1] - psi_midday)
  
  # update midday states
  PL_ST_md[i] = PL_ST_pd[i-1] - PLFL[i] 
  PL_WP_md[i] = psi_midday
  
  # calculate RWU removing PLFL
  #RWU[i] = (SO_WP[i-1] - PSIT + (AET - PLFL[i]) / KT) * Kr
  RWU[i] = AET - PLFL[i]
  
  # update soil water
  SWAT[i]  = SWAT[i-1] - RWU[i]
  SO_WC[i] = (SWAT[i] * TH_S) / SWATMAX
  SO_WP[i] = SWP_MvG(SO_WC[i]) / 1000
  
  # calculate next day's predawn states based on today's midday
  psi_predawn = SO_WP[i] + (psi_midday - SO_WP[i]) * exp(-Kr*(1-D)/C)
  
  PLRF[i] = C * (psi_predawn - psi_midday)
  
  # remove plant refill from soil water
  SWAT[i]  = SWAT[i] - PLRF[i]
  SO_WC[i] = (SWAT[i] * TH_S) / SWATMAX
  SO_WP[i] = SWP_MvG(SO_WC[i]) / 1000
  
  ## update pre-dawn states of following day
  PL_ST_pd[i] = PL_ST_md[i] + PLRF[i]
  PL_WP_pd[i] = psi_predawn
  
}

solution = data.frame(n=1:n, RWU, PLFL, PLRF, PL_ST_pd, PL_ST_md, 
                      PL_WP_pd, PL_WP_md, SWAT, SO_WC, SO_WP)

# derived quantities
solution$AET = solution$RWU + solution$PLFL
solution$f_store = solution$PLFL / solution$AET

# fluxes
ggplot(solution, aes(n, RWU, color="RWU"))+geom_line()+
  geom_line(aes(n, PLFL, color="PLFL"))+geom_line(aes(n, PLRF, color="PLRF"))+
  geom_line(aes(n, AET, color="AET"))

# stores
ggplot(solution, aes(n, PL_ST_pd, color="PL_ST_pd"))+geom_line()+
  geom_line(aes(n, PL_ST_md, color="PL_ST_md"))+
  geom_line(aes(n, SWAT, color="SWAT"))

# potentials
ggplot(solution, aes(n, PL_WP_pd, color="PL_WP_pd"))+geom_line()+
  geom_line(aes(n, PL_WP_md, color="PL_WP_md"))+
  geom_line(aes(n, SO_WP, color="SO_WP"))



## toy model w/o capacitance
# to compare time to soil water depletion

SO_WP = rep(0, n) # soil storage potential (MPa)
SO_WC = rep(0, n) # soil volumetric water content (%)
SWAT  = rep(0, n) # soil water storage (mm)

# initialize states
SO_WP[1] = PSI0
SO_WC[1] = VWC_MvG(SO_WP[1]*1000) # use MvG to convert initial water potential to vwc
SWAT[1]  = SWATMAX * (SO_WC[1] / TH_S)

# fluxes
RWU = rep(0, n) # root water uptake (mm/day)

for (i in 2:n) {
  
  # soil water available for transpiration
  SUPPLY = (SO_WP[i-1] - PSI_CR) * (Kr + Kx)
  
  # actual ET
  AET = min(SUPPLY, PET[i])
  
  # distribute uptake according to conductance in each store
  RWU[i] = AET
  
  # update soil water
  SWAT[i]  = SWAT[i-1] - RWU[i]
  SO_WC[i] = (SWAT[i] * TH_S) / SWATMAX
  SO_WP[i] = SWP_MvG(SO_WC[i]) / 1000
  
}

solu = data.frame(n=1:n, RWU, SWAT, SO_WC, SO_WP)

ggplot(solu, aes(n, RWU, color="RWU"))+geom_line()

ggplot(solu, aes(n, SWAT, color="SWAT"))+geom_line()

ggplot(solu, aes(n, SO_WP, color="SO_WP"))+geom_line()


## ChatGPT sub-daily model for comparison

# --- PARAMETERS ---

dt      <- 3600                 # timestep (seconds, 1 hour)
nsteps  <- 24                   # 24 hours
time    <- seq(0, 23, by=1)

Kr <- 2e-5                      # root conductance (kg s-1 MPa-1)
Kx <- 5e-5                      # xylem conductance
Cp <- 5                         # plant capacitance (kg MPa-1)

psi_soil <- -0.5                # soil water potential (MPa)

# transpiration parameters
Tmax <- 2e-4                    # max transpiration rate (kg s-1)

# --- STORAGE VECTORS ---

psi_p  <- numeric(nsteps)       # plant potential
psi_l  <- numeric(nsteps)       # leaf potential
Fr     <- numeric(nsteps)       # root uptake
Fx     <- numeric(nsteps)       # xylem flow
Tplant <- numeric(nsteps)

# initial condition (predawn equilibrium)
psi_p[1] <- psi_soil
psi_l[1] <- psi_soil

# --- DIURNAL TRANSPIRATION FUNCTION ---

transpiration_demand <- function(hour){
  # simple sinusoidal daytime demand (6–18 h)
  if(hour >= 6 & hour <= 18){
    return(Tmax * sin(pi*(hour-6)/12))
  } else {
    return(0)
  }
}

# --- TIME LOOP ---

for(i in 2:nsteps){
  
  Tplant[i] <- transpiration_demand(time[i])
  
  # xylem flow equals transpiration demand
  Fx[i] <- Tplant[i]
  
  # root flow based on potential gradient
  Fr[i] <- Kr * (psi_soil - psi_p[i-1])
  
  # plant capacitance ODE
  dpsi_dt <- (Fr[i] - Fx[i]) / Cp
  
  psi_p[i] <- psi_p[i-1] + dpsi_dt * dt
  
  # leaf potential from xylem conductance
  psi_l[i] <- psi_p[i] - Fx[i]/Kx
}

# --- DAILY INTEGRALS ---

total_T  <- sum(Tplant) * dt
total_Fr <- sum(Fr) * dt

storage_change <- Cp * (psi_p[nsteps] - psi_p[1])

cat("Total transpiration:", total_T, "\n")
cat("Total root uptake:", total_Fr, "\n")
cat("Storage change:", storage_change, "\n")
cat("Mass balance check (T - Fr - dW):",
    total_T - total_Fr - storage_change, "\n")

psi_pd_sub  <- psi_p[1]
psi_md_sub  <- min(psi_p)

storage_sub <- Cp * (psi_pd_sub - psi_md_sub)

# ===============================
# DAILY ANALYTICAL SOLUTION
# ===============================

D <- 12 * 3600                       # daylight seconds
omega <- pi / D

Ktot <- Kr + Kx
tau  <- Cp / Ktot

Tmean <- total_T / D
A     <- Tmax

# analytical midday ψ drop
psi_md_daily <- psi_soil -
  (Tmean / Kr) -
  (A / sqrt((Ktot)^2 + (Cp*omega)^2))

storage_daily <- Cp * (psi_soil - psi_md_daily)

# ===============================
# OUTPUT COMPARISON
# ===============================

cat("\n--- SUB-DAILY ---\n")
cat("Midday psi:", psi_md_sub, "\n")
cat("Storage:", storage_sub, "\n")

cat("\n--- DAILY ANALYTICAL ---\n")
cat("Midday psi:", psi_md_daily, "\n")
cat("Storage:", storage_daily, "\n")

cat("\n--- DIFFERENCE ---\n")
cat("Midday psi diff:", psi_md_sub - psi_md_daily, "\n")
cat("Storage diff:", storage_sub - storage_daily, "\n")