# analyze Mualem-van Genuchten parameters derived from Puhlmann & von Wilpert PTF

library(tidyverse)
library(lubridate)

swc =  read_csv("../../Data/Pfyn/Pfyn_swat.csv")
swc$VWC = swc$VWC * (1 - 0.62)

swp = read_csv("../../Data/Pfyn/PFY_swpc_corr.csv")
swp$meta = if_else(swp$meta=="irrigated", "irrigation", swp$meta)
swp_long = swp %>% pivot_longer(-c(meta, date), names_to="depth", values_to="SWP")
swp_long$depth = ifelse(swp_long$depth == "SWP_10cm", 10, 80)

sw_comp = inner_join(swp_long, swc)

ggplot(sw_comp, aes(log10(-10*SWP), VWC))+geom_point()+facet_grid(meta~depth)
ggplot(sw_comp, aes(VWC, SWP))+geom_point()+facet_grid(meta~depth)


# PTF constants
ths10 = .4138*(1-0.62) # correct for gravel content
alpha10 = 8.1650
n10 = 1.1769
m10 = 1 - 1 / n10

ths80 = .3839*(1-0.62) # correct for gravel content
alpha80 = 6.005
n80 = 1.2223
m80 = 1 - 1 / n80

# apply MvG formula to test soil water potential transformation
# from LWFBrook90.jl (g converts m to kPa): 
# 9.81 .* (-1 ./ p_MvGα) .* (max.(u_aux_WETNES, eps) .^ (-1 ./ (p_MvGm)) .-1) .^ (1 ./ p_MvGn)
sw_comp$SWP_MvG = if_else(sw_comp$depth == 10,
                          9.81 * (-1 / alpha10) * ((sw_comp$VWC / ths10) ^ (-1 / m10) - 1) ^ (1 / n10),
                          9.81 * (-1 / alpha80) * ((sw_comp$VWC / ths80) ^ (-1 / m80) - 1) ^ (1 / n80))

# compare predicted SWP to measured
ggplot(sw_comp, aes(SWP, SWP_MvG))+geom_point()+facet_grid(meta~depth)+
  geom_abline(slope=1, intercept=0)

ggplot(sw_comp, aes(log10(-10*SWP), VWC))+geom_point()+
  geom_point(aes(log10(-10*SWP_MvG), VWC), color="red")+facet_grid(meta~depth)

ggplot(sw_comp, aes(VWC, SWP))+geom_point()+
  geom_point(aes(VWC, SWP_MvG), color="red")+facet_grid(meta~depth)


# test parameters
thst10 = .2
alphat10 = 6
nt10 = 1.2
mt10 = 1 - 1 / nt10

thst80 = .2
alphat80 = 5
nt80 = 1.22
mt80 = 1 - 1 / nt80

sw_comp$SWP_test = if_else(sw_comp$depth == 10,
                          9.81 * (-1 / alphat10) * ((sw_comp$VWC / thst10) ^ (-1 / mt10) - 1) ^ (1 / nt10),
                          9.81 * (-1 / alphat80) * ((sw_comp$VWC / thst80) ^ (-1 / mt80) - 1) ^ (1 / nt80))

ggplot(sw_comp, aes(log10(-10*SWP), VWC))+geom_point()+
  geom_point(aes(log10(-10*SWP_MvG), VWC), color="red")+
  geom_point(aes(log10(-10*SWP_test), VWC), color="blue")+
  facet_grid(meta~depth)+xlim(1, 5)

ggplot(sw_comp, aes(VWC, SWP))+geom_point()+
  geom_point(aes(VWC, SWP_MvG), color="red")+
  geom_point(aes(VWC, SWP_test), color="blue")+
  facet_grid(meta~depth)



# basic sensitivity
vwc = seq(.2, 0.5, .01)

n = seq(1.16, 1.3, .02)
for (i in 1:length(n)) {
  
  m = 1 - 1 / n[i]
  swp_n = 9.81 * (-1 / alpha10) * (vwc ^ (-1 / m) - 1) ^ (1 / n[i])
  
  if (i==1)
    out_n = data.frame(vwc=vwc, swp=swp_n, n=n[i])
  else
    out_n = rbind(out_n, data.frame(vwc=vwc, swp=swp_n, n=n[i]))
      
}

ggplot(out_n, aes(log10(-10*swp), vwc, group=n, color=n))+geom_line()
ggplot(out_n, aes(vwc, swp, group=n, color=n))+geom_line()


alp = seq(6, 12, .2)
for (i in 1:length(alp)) {
  
  swp_a = 9.81 * (-1 / alp[i]) * (vwc ^ (-1 / mt10) - 1) ^ (1 / nt10)
  
  if (i==1)
    out_a = data.frame(vwc=vwc, swp=swp_a, a=alp[i])
  else
    out_a = rbind(out_a, data.frame(vwc=vwc, swp=swp_a, a=alp[i]))
  
}

ggplot(out_a, aes(log10(-10*swp), vwc, group=a, color=a))+geom_line()
ggplot(out_a, aes(vwc, swp, group=a, color=a))+geom_line()

