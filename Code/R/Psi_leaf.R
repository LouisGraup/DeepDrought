# exploratory analysis of leaf water potential measurements

library(tidyverse)
library(lubridate)

LWP_Bhutan = read_csv("../../Data/Bhutan/LWP_Bhutan.csv")
LWP_Bhutan$Date = if_else(LWP_Bhutan$Date == "10.07.2026", "09.07.2026", LWP_Bhutan$Date) # shift date for predawn-midday comparison

LWP_Visp = read_csv("../../Data/Visp/LWP_Visp.csv")
LWP_Visp$Date = if_else(LWP_Visp$Date == "11.07.2026", "10.07.2026", LWP_Visp$Date) # shift date for predawn-midday comparison
LWP_Visp$TreeNr = LWP_Visp$Remarks # use shorthand notation instead of tree #

LWP_comp = rbind(LWP_Bhutan, select(LWP_Visp, -Remarks))
LWP_comp$Date = as.Date(LWP_comp$Date, format="%d.%m.%Y")
LWP_comp$Species = case_when(LWP_comp$Species == "Eiche" ~ "Oak",
                             LWP_comp$Species == "Esche" ~ "Ash",
                             LWP_comp$Species == "Fichte" ~ "Spruce",
                             LWP_comp$Species == "Tanne" ~ "Fir",
                             LWP_comp$Species == "Foehre" ~ "Pine",
                             .default = "Other")

ggplot(LWP_comp, aes(Species, LWP_mean, group=interaction(Type, Species), fill=Species, alpha=Type))+geom_boxplot()+
  facet_wrap(~Site, nrow=1, scales="free_x")+theme_bw()+scale_alpha_manual(values=c(0.2, 1.0))+
  guides(alpha=guide_legend(override.aes = list(fill="black")))+
  labs(y="LWP (bar)")+ggtitle("Leaf water potential (LWP) measurements comparison")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        strip.text.x=element_text(size=14), plot.title=element_text(size=16, hjust=.5))

ggplot(LWP_comp, aes(TreeNr, LWP_mean, group=interaction(Type, Species), color=Species, shape=Type))+geom_point(size=2.5)+
  facet_wrap(~Site, nrow=1, scales="free_x")+theme_bw()+
  labs(x="Tree ID", y="LWP (bar)")+ggtitle("Leaf water potential (LWP) measurements comparison")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        strip.text.x=element_text(size=14), plot.title=element_text(size=16, hjust=.5))

# predawn - midday difference
LWP_pdmd = LWP_comp %>% select(Site, Date, Type, TreeNr, Species, LWP_mean) %>% 
  pivot_wider(names_from=Type, values_from=LWP_mean) %>% mutate(LWP_diff = Midday - Predawn)

ggplot(LWP_pdmd, aes(Species, LWP_diff, fill=Species))+geom_boxplot()+
  facet_wrap(~Site, nrow=1, scales="free_x")+theme_bw()+
  labs(y="LWP difference (Midday - Predawn) (bar)")+ggtitle("Leaf water potential (LWP) measurements comparison")+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), 
        strip.text.x=element_text(size=14), plot.title=element_text(size=16, hjust=.5))