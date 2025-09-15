#TDP----
#topic: functions for calculating sap flow for the pfynwald project
#date: 24-04-2024
#person: Richard L. Peters

#sapwood depth from the diameter----
SWD<-function(dbh=20,class="Control"){
#error checking
if(is.numeric(dbh)==F)stop("dbh is not numerical.")
if(is.character(class)==F)stop("class is not a character.")

#functions
if(class=="Control"){swt<-1.09125+0.16795*dbh} #fitted model parameters with a mixed-effect model
if(class=="Irrigation"){swt<-0.07907+1.09125+0.16795*dbh}
if(class=="Irrigation stop"){swt<-1.09125-0.13645+0.16795*dbh}

#output
return(swt)
}

#bark thickness from the diameter----
BDD<-function(dbh=20){
#error checking
if(is.numeric(dbh)==F)stop("dbh is not numerical.")

#functions
bt<-  -0.20093+0.11157*dbh  

#output
return(bt)
}

#set fixed parameters----
fixed_param<-function(){
#(-) Beta Correction Value for Wounding [Correction parameter for wounding.]
B_param            <-data.frame(wound_width_cm=c(0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.26,0.28,0.30),
                                B=c(1.7023,1.7585,1.8265,1.8905,1.9572,2.0267,2.0991,2.1482,2.2817,2.4467, 2.5985)) #B parameters reported by Burgess et al. 2001 (https://watermark.silverchair.com/21-9-589.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA3AwggNsBgkqhkiG9w0BBwagggNdMIIDWQIBADCCA1IGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMNYrIfCXiF4aMpi6ZAgEQgIIDI3QBF8VKeVn766NI7QVmJixE014uZoTkgGDMIrYliOw3uURln5Nh_mqFRcgFKVvGHSqm2hM8j5f_DDtm1lcfwYdmLg5buz_mePEa5ljeW6wDcfiF_oPPlI6tYhfaE0MTv0SQ-tBGikn_2ge4T90i9PAOFS5gSRYPElpxpcGYDRSiAdzivnliW4DVEXukQgutEB_MqcrcU323WsPFz23C4Bmdm3g4cyxzpFkiTHPubfd13LLfg3gN2pWlaShpeh3KNsFEgbNktMZkNTOGL9usqaespJTS3GYh2WWBE6EuwhBPNJM0VmQ-23T4JRwpc72SHNbXxlmwyF2clB3XMP0r5AVMtW95TpFndJVN8KUuzhUQhMyE1rhmJ_ESKdOtF-0moRDklVmjYZIY7iLENCqHWMvgPvditUuTgRnklFFySrPVx-2RDW06u59wJY8Rxa8K_ReLp1HKioy4ZP13veMJDQUognqLqLB3zU-JdGwh26WCYn9sQPEgYDjGHsDs631_UiRKbBwxyxt_DTbG1EP43c7p0t4aXdtDlDSl2MFcr0GDCIk9Ojsx_iR7FL-Hw7h2JWaAgmKulOeHjopGrdvJVOwMVODH8kfjqo0VCp-r62wCMKMycmeKGjV3c3R5Y-Awyyi8pIA3Ffo5vhby9DzR9Ik5RzxjpSMT0z197Bsmn_ECkdH0PY2A1C91K-6I3-0xsEebG8zc88TgJS2JoeGQdIF3SSH2Hm_3pbGDpioHPiSeHqPJuC00IdAvXby7dkVC1zkZLa_St_7n1wmRIAmIXUIB9D_7ROVS5SRIHFKqXyvKR-GfKpsHWbaxNNZzsd0Ic0lVfseGZV0nZ1jZZ0eVnHUUjjXad1UoAZFnTPPXrz2QY0HBd7PI9UGcVyvotoSnrWi8XoyTZdhDS26q-o4tXKlMxfqdJEnL24iou9ce-iWTdesN0i1J7hNgKA3rKT_nWjQ64huzDg1hU1kQxtjdstGzBy9xHKaIBv7bQVpwmUCMPhx2mAc5LZ08n6U9-g_h0_ZpD6wjx2QcHhoSVaFMjIC5oV_kmJsHfKmhpW16uOxwf9SR)

#we do not correct for the wound diameter yet
d_w   <- 0.21 #(cm) Wound Diameter (2-dec. place ONLY) [The wound developed in the sapwood due to insertion of needles.]

sel_row            <-which(B_param$wound_width_cm-round(d_w,2)==0)
if(length(sel_row)==0){print("No wound correction factor available for this value (using default)")
B                  <- 1.8905 }else{B<- B_param[sel_row,"B"]}
c_d                <- 1200 #(J/kg/°C at 20 °C) Wood Matrix Specific Heat Capacity [Wood Matrix Specific Heat Capacity.]
c_w                <- 4182 #(J/kg/°C at 20 °C) Specific Heat Capacity of Sap [Sap Specific Heat Capacity.]
ρ_w                <- 1000 #(kg/m3) Density of Water [Density of Water.]
ρ_cw               <- 1530 #(kg/m3) Cell Wall Density [Density of Cell Wall.]
K_w                <- 0.5984 #(J/kg/°C at 20 °C) Thermal Conductivity of Water [Thermal conductivity of water (Lide 1992).]

#output
out<-data.frame(B=B,c_d=c_d,c_w=c_w,ρ_w=ρ_w,ρ_cw=ρ_cw,K_w=K_w)
return(out)
}

#measured parameters----
meas_param<-function(class="Control",TTD=1,BDD=1,SWD=1){
#error checking
if(is.character(class)==F)stop("class is not a character.")
if(is.numeric(TTD)==F)stop("TTD is not numerical.")
if(is.numeric(BDD)==F)stop("BDD is not numerical.")
if(is.numeric(SWD)==F)stop("SWD is not numerical.")
  
#function
if(class=="Control"){ρ_d<-476.99} #fitted model parameters with a mixed-effect model
if(class=="Irrigation"){ρ_d<-476.99+47.57}
if(class=="Irrigation stop"){ρ_d<-476.99+47.57}
ρ<-1019.393
c<-2765.90
mc<-1.12051

#adding all other parameters
fp                 <-data.frame(fixed_param())

B                  <-fp$B
c_d                <-fp$c_d
c_w                <-fp$c_w
ρ_w                <-fp$ρ_w
ρ_cw               <-fp$ρ_cw
K_w                <-fp$K_w

#algebraic parameters----
ρ_d                <- ρ_d #(kg/m3) Basic Density of Dry Wood [Basic Density of Dry Wood.]
ρ                  <- ρ #(kg/m3) Basic Density of Fresh Wood [Basic Density of Fresh Wood.]
mc                 <- mc #(kg/kg) Gravimetric Water Content of Sapwood [Gravimetric Water Content of Sapwood.]
mc_v               <- mc*(ρ_d/1000) #(m3/m3) Volumetric Water Content of Sapwood [Volumetric Water Content of Sapwood]
mc_FSP             <- 0.2*(ρ_d*(1/ρ_w))^-0.5 #(-) Water Content at Fibre Saturation Point [Water Content of Sapwood at Fibre Saturation Point.]
Fv                 <- 1-ρ_d*((0.6536+mc)/1000) #(-) Void Fraction of Wood [Void fraction of wood (Swanson 1983) indicats the volume fraction of air in the wood.]
K_d                <- 0.04182*(21-(20*Fv)) #(J/kg/°C at 20 °C) Thermal Conductivity of Dry Wood Matrix [Thermal Conductivity of Dry Wood Matrix (Swanson 1983).]
Fv_FSP             <- 1-(ρ_d/ρ_w)*((ρ_w/ρ_cw)+mc_FSP) #(-) Void Fraction of Wood at Fibre Saturation Point [Void fraction of wood at fibre saturation point from Vandegehuchte & Steppe (2012)  indicats the volume fraction of air in the wood.]
c                  <- c #(J/kg/°C at 20 °C) Specific Heat Capacity of Fresh Wood [Specific heat capacity of fresh wood (Edwards and Warwick 1984).]
ρc                 <- (ρ)*(c) #(J/m3/°C) Volumetric Specific Heat Capacity of Fresh Wood [Volumetric Specific Heat Capacity of Fresh Wood.]
K_Van              <- K_w*(mc-mc_FSP)*(ρ_d/ρ_w)+0.04186*(21-20*Fv_FSP) #(J/kg/°C at 20 °C) Thermal Conductivity (K_Van) [Stem thermal conductivity from Vandegehuchte & Steppe (2012).]
HWA                <- pi*((TTD/2)-BDD-if(SWD<1.51){SWD}else{2})^2 #(cm2) Heartwood Area [Heartwood area.]
ISWA               <- if(SWD<1.51){0}else{pi*((TTD/2)-BDD-1)^2-HWA} #(cm2) Inner Sapwood Area [Inner Sapwood Area equals trunk radius (TDD/2) minus bark depth (BDD) minus heartwood area. When sapwood depth < 1.51cm, inner sapwood area equals zero.]
OSWA               <- pi*((TTD/2)-BDD)^2-HWA-ISWA #(cm2) Outer Sapwood Area [Outer Sapwood Area equals trunk radius (TDD/2) minus bark depth (BDD) minus inner sapwood area minus heartwood area.]
k_Van              <- ((K_Van/ρc))*10000 #(cm2/s) Thermal Diffusivity (k_Van) [The rate of diffusion of heat in the wood and sap matrix as calculated from Vandegehuchte & Steppe (2012). (Multiplied by 10000 to convert m2/s to cm2/s)]

#output  
out<-data.frame(
  B=B,c_d=c_d,c_w=c_w,ρ_w=ρ_w,ρ_cw=ρ_cw,K_w=K_w,ρ_d=ρ_d,ρ=ρ,mc=mc,mc_v=mc_v,mc_FSP=mc_FSP,Fv=Fv,K_d=K_d,Fv_FSP=Fv_FSP,c=c,ρc=ρc,K_Van=K_Van,HWA=HWA,ISWA,OSWA=OSWA,k_Van=k_Van,
  TTD=TTD,BDD=BDD,SWD=SWD,class=class)
return(out)
}

#Vc output----
vc_cal<-function(raw_data,param,x_d_o=0.6,x_d_i=0.6,x_u_o=0.6,x_u_i=0.6){

#algebraic parameters  
tm_outer           <- round((x_d_o^2)/(4*param$k_Van),1) #(seconds) Outer Tmax under Zero Flow [Time to maximum temperature rise under zero sap flow conditions in the outer position.]
tm_inner           <- round((x_d_i^2)/(4*param$k_Van),1) #(seconds) Inner Tmax under Zero Flow [Time to maximum temperature rise under zero sap flow conditions in the inner position.]
tm_Peclét_outer    <- round(0.9696*tm_outer-3.2363,1) #(seconds) Outer Peclét Tmax [The Tmax value to transition between slow and fast flow for the outer position.]
tm_Peclét_inner    <- round(0.9696*tm_inner-3.2363,1) #(seconds) Inner Peclét Tmax [The Tmax value to transition between slow and fast flow for the inner position.]

#measured sap flux density expressed in cm/h
Vc_HRM_Outer <- (((2*param$k_Van)/(x_d_o+x_u_o))*raw_data$alpha_outer+((x_d_o-x_u_o)/(2*58.5)))*param$B*3600 #(cm/hr)
Vc_HRM_Inner <- (((2*param$k_Van)/(x_d_i+x_u_i))*raw_data$alpha_inner+((x_d_i-x_u_i)/(2*58.5)))*param$B*3600 #(cm/hr)
#issue with NA values as the fluxes can sometime be to low as it takes to long to reach the maximum pulse height due to slow fluxes
Vc_Tmax_Outer<- suppressWarnings((sqrt(((4*param$k_Van)/3)*(log(1-(3/raw_data$tMaxT_outer)))+((x_d_o^2)/(raw_data$tMaxT_outer*(raw_data$tMaxT_outer-3)))))*3600)*param$B #(cm/hr)
Vc_Tmax_Inner<- suppressWarnings((sqrt(((4*param$k_Van)/3)*(log(1-(3/raw_data$tMaxT_inner)))+((x_d_i^2)/(raw_data$tMaxT_inner*(raw_data$tMaxT_inner-3)))))*3600)*param$B #(cm/hr)

DMA_Outer    <- NA
for(i in c(1:length(tMaxT_outer))){if(tMaxT_outer[i]<tm_Peclét_outer){DMA_Outer[i]<-Vc_Tmax_Outer[i]}else{DMA_Outer[i]<-Vc_HRM_Outer[i]}} #(cm/hr)
DMA_Inner    <- NA
for(i in c(1:length(tMaxT_inner))){if(tMaxT_inner[i]<tm_Peclét_inner){DMA_Inner[i]<-Vc_Tmax_Inner[i]}else{DMA_Inner[i]<-Vc_HRM_Inner[i]}} #(cm/hr)

#merge all data for the output
Vc<-data.frame(Date.Time=raw_data$Date.Time,Vc_HRM_Outer,Vc_HRM_Inner,Vc_Tmax_Outer,Vc_Tmax_Inner,DMA_Outer,DMA_Inner)
return(Vc)
}

#misalignment----
align<-function(Vc,param,hours=c(2,3),dates=as.Date(c(as.Date("2024-01-28"):as.Date("2024-02-28"))),method="HRM",auto_correct=T){
  
#error checking
if(is.numeric(hours)==F)stop("hours provided are not numeric.")
if(length(which(method%in%c("HRM","Tmax","DMA")))==0)stop("select either HRM, Tmax, or DMA.")

#prepare input data
if(method=="HRM"){time_series<-Vc[,c(1,2,3)]}
if(method=="Tmax"){time_series<-Vc[,c(1,4,5)]}
if(method=="DMA"){time_series<-Vc[,c(1,6,7)]}

#relevant functions
right = function (string, char) {substr(string,nchar(string)-(char-1),nchar(string))}
left = function (string,char) {substr(string,1,char)}

#isolate the values
time_series$hours<-as.numeric(left(right(format(time_series$Date.Time),8),2))
time_series$date <-as.Date(left(time_series$Date.Time,10)) 
sel_time         <-time_series[which(as.character(time_series$hours)%in%as.character(hours)),]
sel_time         <-sel_time[which(as.character(sel_time$date)%in%as.character(dates)),]

#inner probe
if(nrow(na.omit(sel_time[,c(1,3,4,5)]))==0){

if(auto_correct==T){
sel_time         <-time_series[which(as.character(time_series$hours)%in%as.character(hours)),]
sel_time         <-na.omit(sel_time[,c(1,3,4,5)])  
agg_time<-aggregate(sel_time[,-1],by=list(sel_time$date),FUN=mean,na.rm=T)  
agg_sel <-agg_time[which(as.character(agg_time$date)%in%as.character(as.Date(c(as.Date((agg_time)[1,"date"]):as.Date(na.omit(agg_time)[1,"date"]+30))))),]
cor_inner<-as.numeric(quantile(agg_sel$Vc_HRM_Inner,probs=0.1,na.rm=T))

if(method=="HRM"){
  Vc[,3]<-Vc[,3]-cor_inner
}
if(method=="Tmax"){
  Vc[,5]<-Vc[,5]-cor_inner
}
if(method=="DMA"){
  Vc[,7]<-Vc[,7]-cor_inner  
}  
  
}else{
  #output generation
  if(method=="HRM"){
    Vc[,3]<-Vc[,3]
  }
  if(method=="Tmax"){
    Vc[,5]<-Vc[,5]
  }
  if(method=="DMA"){
    Vc[,7]<-Vc[,7]  
  }
}  
}else{

#calculate the daily means
agg_time<-aggregate(sel_time[,-1],by=list(sel_time$date),FUN=mean,na.rm=T)
cor_inner<-as.numeric(quantile(agg_time$Vc_HRM_Inner,probs=0.1,na.rm=T))#min(agg_time$Vc_HRM_Inner,na.rm=T)

#output generation
if(method=="HRM"){
  Vc[,3]<-Vc[,3]-cor_inner
}
if(method=="Tmax"){
  Vc[,5]<-Vc[,5]-cor_inner
}
if(method=="DMA"){
  Vc[,7]<-Vc[,7]-cor_inner  
}
}

#outer
sel_time         <-time_series[which(as.character(time_series$hours)%in%as.character(hours)),]
sel_time         <-sel_time[which(as.character(sel_time$date)%in%as.character(dates)),]

if(nrow(na.omit(sel_time[,c(1,2,4,5)]))==0){
  
  if(auto_correct==T){
    sel_time         <-time_series[which(as.character(time_series$hours)%in%as.character(hours)),]
    sel_time         <-na.omit(sel_time[,c(1,2,4,5)])  
    agg_time<-aggregate(sel_time[,-1],by=list(sel_time$date),FUN=mean,na.rm=T)  
    agg_sel <-agg_time[which(as.character(agg_time$date)%in%as.character(as.Date(c(as.Date((agg_time)[1,"date"]):as.Date(na.omit(agg_time)[1,"date"]+30))))),]
    cor_outer<-as.numeric(quantile(agg_sel$Vc_HRM_Outer,probs=0.1,na.rm=T))
    
    if(method=="HRM"){
      Vc[,2]<-Vc[,2]-cor_outer
    }
    if(method=="Tmax"){
      Vc[,4]<-Vc[,4]-cor_outer
    }
    if(method=="DMA"){
      Vc[,6]<-Vc[,6]-cor_outer
    }
    
  }else{
    #output generation
    if(method=="HRM"){
      Vc[,2]<-Vc[,2]
    }
    if(method=="Tmax"){
      Vc[,4]<-Vc[,4]
    }
    if(method=="DMA"){
      Vc[,6]<-Vc[,6]
    }
  }  
}else{
  
  #calculate the daily means
  agg_time<-aggregate(sel_time[,-1],by=list(sel_time$date),FUN=mean,na.rm=T)
  cor_outer<-as.numeric(quantile(agg_time$Vc_HRM_Outer,probs=0.1,na.rm=T))#min(agg_time$Vc_HRM_Inner,na.rm=T)
  
  #output generation
  if(method=="HRM"){
    Vc[,2]<-Vc[,2]-cor_outer
  }
  if(method=="Tmax"){
    Vc[,4]<-Vc[,4]-cor_outer
  }
  if(method=="DMA"){
    Vc[,6]<-Vc[,6]-cor_outer
  }
}

return(Vc)

#TRASH CODE
#CORRECTION DOES NOT WORK ALWAYS FOR THE TMAX METHOD DUE TO LOW FLUX
# x_d Correction with TMax Equation
# TMax Equation: (4*k_van/3 * ln(1-3*tMaxT) + (x_d^2/(tMaxT(tMaxT-3))))^0.5 *B*3600
# X_d = (k_van*4/3*ln(1-3/tMaxT)*(-1)*tMaxT(tMaxT-3))^0.5
#x_d_oCorr <- (((param$K_Van*4/3)*log(1-3/agg_time$tMaxT_outer)*(-1)*agg_time$tMaxT_outer*(agg_time$tMaxT_outer-3)))^0.5
#x_d_iCorr <- (((param$K_Van*4/3)*log(1-3/agg_time$tMaxT_inner)*(-1)*agg_time$tMaxT_inner*(agg_time$tMaxT_inner-3)))^0.5

# x_u Correction with HRM Equation and corrected x_d 
# HRM Equation: Vc_HRM_Outer (cm/hr): (((2*k_van)/(x_d_o+x_u_o))*alphaOuter+((x_d_o-x_u_o)/(2*58,5)))*B*3600
# x_u = ((4*k_van**alpha*58.5 - x_d^2)*(-1))^0.5
#x_u_oCorr <- (((4*param$K_Van*agg_time$alpha_outer*58.5)+(x_d_oCorr)^2))^0.5
#x_u_iCorr <- (((4*param$K_Van*agg_time$alpha_inner*58.5)+(x_d_iCorr)^2))^0.5
}

#J output----
j_cal<-function(Vc,param,method="HRM"){

#error checking
if(length(which(method%in%c("HRM","Tmax","DMA")))==0)stop("select either HRM, Tmax, or DMA.")

#prepare input data
if(method=="HRM"){sel_input<-Vc[,c(1,2,3)]}
if(method=="Tmax"){sel_input<-Vc[,c(1,4,5)]}
if(method=="DMA"){sel_input<-Vc[,c(1,6,7)]}

Sap_Flux_Density_Outer <- (sel_input[,2]*param$ρ_d*(param$c_d+(param$mc*param$c_w)))/(param$ρ_w*param$c_w) #(cm3/cm2/hr)
Sap_Flux_Density_Inner <- (sel_input[,3]*param$ρ_d*(param$c_d+(param$mc*param$c_w)))/(param$ρ_w*param$c_w) #(cm3/cm2/hr)

#generate output
J                       <-data.frame(Date.Time=sel_input$Date.Time,Sap_Flux_Density_Outer,Sap_Flux_Density_Inner)
J                       <-J[order(J$Date.Time),]
return(J)
}

#Q output----
q_cal<-function(J,param){
#input
input_data=J
param=param
  
#functions
Total_Sap_Flow <- ((input_data$Sap_Flux_Density_Outer*param$OSWA)+(input_data$Sap_Flux_Density_Inner*param$ISWA))/1000 #(L/hr or kg/hr)

#output
Q              <-data.frame(Date.Time=input_data$Date.Time,Total_Sap_Flow)
Q              <-Q[order(Q$Date.Time),]
return(Q)
}

#tree_spec----
tree_spec<-function(tree_id=796){
#process
tree_id<-as.character(tree_id)

base=data.frame(
  species="Pinus sylvestris",
  tree   =c("795"           ,  "796"           ,  "797"           ,  "38"     ,  "40"     ,   "41"     ,   "597"    ,  "600"    ,  "601"    ,  "235"       ,  "237"       ,  "241"       ,  "389"       ,  "391"       ,  "392"       ,  "755"       ,  "756"       ,  "757"       ,  "466"     ,  "469"     ,  "472"     ,  "647"     ,  "677"     ,  "679"     , "913"     ,  "914"     ,  "919"     ,  "254"           ,  "263"          ,  "332"           ,  "819"           ,  "820"           ,  "937"           ,  "518" ,  "522" ,  "572" ,   "669",   "705",  "706",  "703" ,  "1060", "1086", "1037", "1038", "1039","1085","671"),
  dbh    =c(25.8            ,  20.8            ,  21.7            ,   35.1    ,   18.9    ,   16.5     ,      28.3  ,     21.3  ,     16.3  ,     22.2     ,     21.1     ,     28       ,     31.2     ,    31.5      ,    31        ,     23.9     ,     27.8     ,     33.5     ,     20.3   ,     24.4   ,     29.1   ,     18.8   ,     19.6   ,     20     ,    26     ,     18.4   ,     24.8   ,     27.3         ,     37.1        ,     35.3         ,     32.8         ,     30.4         ,     33.5         ,    22.2,    17.4, 25.8   ,    36.4,    32.2,   24.9,   22.5 ,    20.5,   18.7,  23.4 ,   23.6,   28.4,  32 , 30),
  treat  =c("irrigation_vpd",  "irrigation_vpd",  "irrigation_vpd",  "control",  "control",   "control",   "control",  "control",  "control",  "irrigation",  "irrigation",  "irrigation",  "irrigation",  "irrigation",  "irrigation",  "irrigation",  "irrigation",  "irrigation",  "roof_vpd",  "roof_vpd",  "roof_vpd",  "roof_vpd",  "roof_vpd",  "roof_vpd", "roof_vpd",  "roof_vpd",  "roof_vpd",  "irrigation_vpd", "irrigation_vpd",  "irrigation_vpd",  "irrigation_vpd",  "irrigation_vpd",  "irrigation_vpd",  "roof",  "roof",  "roof",  "roof",  "roof", "roof",  "roof", "roof" , "roof", "control", "control", "control","roof","roof"),
  scaf   =c(14              ,  14              ,  14              ,  1        ,  1        ,   1        ,   701      ,  701      ,  701      ,  3           ,  3           ,  3           ,  4           ,  4           ,  4           ,  901         ,  901         ,  901         ,  6         ,  6         ,  6         ,  12        ,  12        ,  12        , 16        ,  16        ,  16        ,  10              ,  10             ,  10              ,  15              ,  15              ,             15   ,  11    ,  11    ,  11    ,  13    ,    13  , 13    ,  13    ,  17    , 17    , 18       , 18       , 18 ,  17, 13)
  )

#error
if(length(which(base$tree==tree_id))==0)stop("no such tree in the database.")

#output
sel<-base[which(base$tree==tree_id),]
return(sel)
}

#corrections----
corrections<-function(Vc,tree_id=796){

#add the corrections to the table
cortab<-data.frame(matrix(
c(5,40	,"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T13:00:00Z",	"remove",	NA,
  5,40	,"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T13:00:00Z",	"remove",	NA,
  9,601,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T11:40:00Z",	"remove",	NA,
  9,601,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T11:30:00Z",	"remove",	NA,
  9,601,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T11:00:00Z",	"shift",	-13.9609376,
  9,601,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T11:00:00Z",	"scale",	14.47443275,
  12,	241,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T13:10:00Z",	"remove",	NA,
  12,	241,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T13:10:00Z",	"remove",	NA,
  15,	755,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T14:10:00Z",	"remove",	NA,
  15,	755,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T14:10:00Z",	"remove",	NA,
  23,	677,	"Vc_HRM_Inner",	"2024-04-07T13:20:00Z",	"2024-04-07T13:20:00Z",	"remove",	NA,
  34,	518,	"Vc_HRM_Outer",	"2024-03-18T12:20:00Z",	"2024-03-18T12:20:00Z",	"remove",	NA,
  34,	518,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T12:00:00Z",	"shift",	-9.768614547,
  34,	518,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T12:00:00Z",	"scale",	36.93692951,
  34,	518,	"Vc_HRM_Inner",	"2024-03-18T12:20:00Z",	"2024-03-18T12:20:00Z",	"remove",	NA,
  34,	518,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T12:00:00Z",	"shift",	-17.36687796,
  34,	518,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T12:00:00Z",	"scale",	45.66227816,
  36,	572,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-17T06:50:00Z",	"remove",	NA,
  37,	669,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T12:20:00Z",	"remove",	NA,
  37,	669,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T12:20:00Z",	"remove",	NA,
  44,	1038,	"Vc_HRM_Outer",	"2024-01-01T00:00:00Z",	"2024-03-18T12:20:00Z",	"remove",	NA,
  44,	1038,	"Vc_HRM_Inner",	"2024-01-01T00:00:00Z",	"2024-03-18T12:20:00Z",	"remove",	NA),ncol=7,byrow=T))
colnames(cortab)<-c("tree","tree_id","sensor","start","end","action","value")
cortab$value<-as.numeric(cortab$value)
#

#check if there are corrections for the tree
rows<-which(as.character(cortab$tree_id)==as.character(tree_id))
if(length(rows)==0){
return(Vc)  
break
}

#select the rows and build a loop which can handle the actions
sel_cor<-cortab[rows,]

#loop for each item
for(i in c(1:nrow(sel_cor))){
ssel<-sel_cor[i,]

#removing values
if(ssel$action=="remove"){
start<-which(Vc$Date.Time==as.POSIXct(ssel$start,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"))
end<-which(Vc$Date.Time==as.POSIXct(ssel$end,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"))  
Vc[c(start:end),ssel$sensor]<-as.numeric(ssel$value)  
}

#shifting values
if(ssel$action=="shift"){
start<-which(Vc$Date.Time==as.POSIXct(ssel$start,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"))
end<-which(Vc$Date.Time==as.POSIXct(ssel$end,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"))  
Vc[c(start:end),ssel$sensor]<-Vc[c(start:end),ssel$sensor]+as.numeric(ssel$value)
}

#scale the values
if(ssel$action=="scale"){
start<-which(Vc$Date.Time==as.POSIXct(ssel$start,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"))
end<-which(Vc$Date.Time==as.POSIXct(ssel$end,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC"))  
start_3<-which(Vc$Date.Time==as.POSIXct(ssel$end,format="%Y-%m-%dT%H:%M:%SZ",tz="UTC")-1*60*60*24*3) #3 days prior to the end we select the maximum
cor_val<-ssel$value/max(Vc[c(start_3:end),ssel$sensor])
Vc[c(start:end),ssel$sensor]<-Vc[c(start:end),ssel$sensor]*cor_val 
}
}

#create output
return(Vc)
}

#radial profile correction----
rad_cor<-function(swd=4,probe.length=2){
  #probe.length=2
  #swd=4
  
  library(pracma)
  beta<- 0.75#0.66#coef(fit)[1]
  gamma1<-1#0.98#coef(fit)[2]
  gamma2<-4.5#2.7#coef(fit)[3]
  
  radial_profile <- function(x, beta, gamma1, gamma2) {
    (x / beta)^(gamma1 - 1) * exp(-(x / beta)^gamma2)
  }
  
  # Integrate over full sapwood depth (0 to 4 cm)
  SFD_full <- integral(function(x) radial_profile(x, beta, gamma1, gamma2), 0, swd) / swd
  
  # Integrate over measured depth (0 to 2 cm)
  SFD_measured <- integral(function(x) radial_profile(x, beta, gamma1, gamma2), 0, probe.length) / probe.length
  
  # Compute correction factor
  correction_factor <- SFD_full / SFD_measured
  return(correction_factor)
}

#sapwood area calculation
calculate_sapwood_area <- function(dbh, swd, bark_thickness) {
  # Compute the wood radius (excluding bark)
  R_wood <- (dbh - 2 * bark_thickness) / 2
  
  # Compute heartwood radius
  R_heartwood <- R_wood - swd
  
  # Ensure heartwood radius is not negative
  R_heartwood <- ifelse(R_heartwood < 0, 0, R_heartwood)
  
  # Compute sapwood area
  sapwood_area <- pi * (R_wood^2 - R_heartwood^2)
  
  return(sapwood_area)
}

#gs calculation
gs_eImp <- function(e,t,vpd,elv){
  #Calculates surface conductance based on imposed evaporation (simplified version of Penman-Monteith for aerodynamically #well-coupled conditions).
  e<-sap_flow_hourly$Sap_Flow_LA_g_m.2_h.1
  t<-sap_flow_hourly$TEMP
  vpd<-sap_flow_hourly$VPD
  elv<-615
  
  #Input: # vpd: vapor pressure deficit [kPa]
  #t: air temperature [°C]
  #elv: elevation [m asl]
  #e: transpiration [g m-2 h-1]
  # Output:
  #gs: canopy conductance [m s-1]
  #Set constants 
  cp <- 1.01                                                                                #Heat capacity of air [J g-1 K-1]
  #Helper variables
  zeroK <- 273.15                                                                      #0 °C in K
  hundredK <- 373.15	                                                #100 °C in K
  gammaT <- 0.005			                                #Temperature lapse rate [k m-1]
  rAir <- 287.04			                                #Gas constant for dry air [J kg-1 K-1]
  rWat <- 461.5	                                 #Gas constant for water vapour
  pZero <- 1013.25		#Standard air pressure
  g <- 9.81			                                #Gravitational acceleration [m s-2]
  #Functions of temperature (adapted from PREVAH subroutine EPMONDT / FORHYCS subroutine EVAP)
  tk <- t+zeroK
  tr <- 1.0-(hundredK/tk)
  tp <- tk+gammaT*elv/2.0		                                #Potential temperature			             
  pt <- pZero * exp(-g*elv/rAir/tp)	                                #Atmospheric pressure
  lambda <- 2500.8 - 2.36*t + 0.0016*(t^2) - 0.00006*(t^3)	#Latent heat of vaporization of water [J g-1]
  gamma <- (pt*cp/lambda/(rAir/rWat))*100.0	                #Psychrometric "constant" [Pa K-1]
  rhoAir <- (pt/(rAir*tk)) * 100000.0	                                #Density of air [g m-3]
  #Unit conversion
  e <- e/3600                                                                            # hours to seconds
  vpd <- vpd*1000  #kPa to Pa
  #vpd[which(vpd<100)]<-NA                                                  #Remove values below 100 Pa
  gs <- (lambda*e*gamma)/(rhoAir*cp*vpd)
  return(gs)}

#relevant function
gs_e <- function(gs,gs_max,t,vpd,elv){
  # Calculates surface conductance based on imposed evaporation (simplified version of Penman-Monteith
  # for aerodynamically well-coupled conditions).
  #
  # Input: 
  # vpd: vapor pressure deficit [kPa]
  # t: air temperature [ C]
  # elv: elevation [m asl]
  # e: transpiration [g m-2 h-1]
  #
  # Output:
  # gs [m s-1]
  # Set constants 
  cp <- 1.01          #Heat capacity of air [J g-1 K-1]
  #Helper variables
  zeroK <- 273.15			#0  C in K
  hundredK <- 373.15	#100  C in K
  gammaT <- 0.005			#Temperature lapse rate [k m-1]
  rAir <- 287.04			#Gas constant for dry air [J kg-1 K-1]
  rWat <- 461.5				#Gas constant for water vapor
  pZero <- 1013.25		#Standard air pressure
  g <- 9.81					  #Gravitational acceleration [m s-2]
  
  
  #Functions of temperature (adapted from PREVAH subroutine EPMONDT / FORHYCS subroutine EVAP)
  tk <- t+zeroK
  tr <- 1.0-(hundredK/tk)
  tp <- tk+gammaT*elv/2.0											                  #Potential temperature
  pt <- pZero * exp(-g*elv/rAir/tp)									            #Atmospheric pressure
  lambda <- 2500.8 - 2.36*t + 0.0016*(t^2) - 0.00006*(t^3)			#Latent heat of vaporization of water [J g-1]
  gamma <- (pt*cp/lambda/(rAir/rWat))*100.0								      #Psychrometric "constant" [Pa K-1]
  rhoAir <- (pt/(rAir*tk)) * 100000.0								            #Density of air [g m-3]
  
  #Unit conversion
  #e <- e/3600  # hours to seconds
  vpd <- vpd*1000  #kPa to Pa
  #vpd[which(vpd<0.01)]<-0.01
  #vpd[which(vpd<100)]<-NA #100 This is the one I used
  #vpd[which(vpd<1)]<-1 #100
  gc <-gs*gs_max
  e<-(rhoAir*cp*vpd*gc)/(lambda*gamma) 
  e<-e*3600#g m-2 h-1
  return(e)
}


#add this function
#data_corr <- data_wide %>%
#  mutate(
#    # Convert SWP from kPa to hPa, then calculate pF
#    pfraw = ifelse(!is.na(SWP) & SWP < 0, log10(10 * abs(SWP)), NA),
#    
#    # Temperature deviation from 22°C
#    delT = ifelse(!is.na(TEM), TEM - 22, NA),
#    
#    # Calculate delta pF correction
#    delpF = ifelse(
#      !is.na(SWP) & !is.na(delT),
#      (6.206 * delT + 0.1137 * (delT)^2) * exp(-22.76 / (pfraw + 0.1913)),
#      NA
#    ),
#    
#    # Compute corrected pF
#    pFcr = ifelse(!is.na(SWP) & !is.na(delpF), pfraw + delpF, NA),
#    
#    # Convert corrected pF back to SWP
#    SWP_corr = ifelse(!is.na(SWP) & !is.na(pFcr), ((10^pFcr) * -1) / 10, NA)
#  )

#calculate gs and PET
calculate_pet_and_gs <- function(e, t, vpd, rad, wind, wind_height = 2, elv = 0) {
  # e    = transpiration rate (g m^-2 h^-1)
  # t    = air temperature (°C)
  # vpd  = vapor pressure deficit (kPa)
  # rad  = net radiation (W/m^2)
  # wind = wind speed (m/s)
  # elv  = elevation (m)
  
  # Constants
  R <- 8.314           # Universal gas constant [J/mol/K]
  Mw <- 0.01801528     # kg/mol
  cp <- 1013           # J/kg/°C
  r_air <- 287.05      # J/kg/K (gas constant for dry air)
  lambda <- 2.45e6     # J/kg
  gamma <- 0.066       # Psychrometric constant [kPa/°C]
  g <- 9.81            # m/s^2
  
  # Pressure at elevation (Pa)
  P <- 101325 * exp(-elv * g / (r_air * (t + 273.15)))  # Pa
  tK <- t + 273.15  # Kelvin
  
  # Adjust wind to 2 m if needed
  if (wind_height != 2) {
    wind <- wind * log(67.8 * 2 - 5.42) / log(67.8 * wind_height - 5.42)
  }
  
  # Slope of saturation vapor pressure curve (kPa/°C)
  Delta <- 4098 * (0.6108 * exp((17.27 * t)/(t + 237.3))) / ((t + 237.3)^2)
  
  # Air density (kg/m^3)
  rho_air <- P / (r_air * tK)  # Ideal gas law
  
  # VPD in Pa
  VPD_Pa <- vpd * 1000
  
  # Aerodynamic resistance (s/m) — approximate
  ra <- 208 / wind  # [s/m], FAO estimate
  
  # PET in kg/m^2/s
  PET_kg_m2_s <- (Delta * rad + rho_air * cp * VPD_Pa / ra) /
    (lambda * (Delta + gamma * (1 + 0.34 * wind)))
  
  # PET in mol/m^2/s
  PET_mol_m2_s <- PET_kg_m2_s / Mw
  
  # Measured transpiration (g/m²/h) → mol/m²/s
  E_kg_m2_s <- e / 1000 / 3600
  E_mol_m2_s <- E_kg_m2_s / Mw
  
  # Stomatal conductance (mol/m²/s)
  gs <- (E_mol_m2_s * P) / (VPD_Pa * R * tK)
  gs[!is.finite(gs)] <- NA
  
  return(list(
    PET_mol_m2_s = PET_mol_m2_s,
    gs_mol_m2_s = gs
  ))
}

#model version 2 (the correct one)
calculate_pet_and_gs <- function(e, t, vpd, rad, wind, wind_height = 2, elv = 0) {
  # Inputs:
  # e    = transpiration rate (g/m²/h)
  # t    = air temperature (°C)
  # vpd  = vapor pressure deficit (kPa)
  # rad  = net radiation (W/m²)
  # wind = wind speed (m/s)
  # elv  = elevation (m)
  
  # ==== CONSTANTS ====
  R      <- 8.314        # J/mol/K
  Mw     <- 0.01801528   # kg/mol
  cp_g   <- 1.01         # J/g/K
  cp_kg  <- 1013         # J/kg/K
  r_air  <- 287.05       # J/kg/K
  r_wat  <- 461.5        # J/kg/K
  g      <- 9.81         # m/s²
  p0     <- 101325       # Pa (sea-level)
  
  # ==== TEMPERATURE ====
  tK <- t + 273.15
  
  # ==== PRESSURE ADJUSTED FOR ELEVATION ====
  tp <- tK + 0.005 * elv / 2
  P_Pa <- p0 * exp(-g * elv / (r_air * tp))
  
  # ==== LATENT HEAT OF VAPORIZATION (J/g) ====
  lambda_g <- 2500.8 - 2.36 * t + 0.0016 * t^2 - 0.00006 * t^3  # J/g
  lambda_kg <- lambda_g * 1000  # J/kg
  
  # ==== PSYCHROMETRIC CONSTANT (Pa/K) ====
  gamma_pa <- (P_Pa / 100) * cp_g / lambda_g / (r_air / r_wat) * 100
  
  # ==== AIR DENSITY ====
  rho_air_kg <- P_Pa / (r_air * tK)
  rho_air_g <- rho_air_kg * 1000
  
  # ==== VPD ====
  VPD_Pa <- vpd * 1000
  
  # ==== TRANSPIRATION RATE ====
  e_g_m2_s <- e / 3600
  
  # ==== STOMATAL CONDUCTANCE (m/s) ====
  gs_m_s <- (lambda_g * e_g_m2_s * gamma_pa) / (rho_air_g * cp_g * VPD_Pa)
  
  # ==== CONVERT TO mol/m²/s ====
  gs_mol_m2_s <- gs_m_s * P_Pa / (R * tK)
  
  # ==== PENMAN–MONTEITH PET ====
  # Adjust wind to 2 m height if needed
  if (wind_height != 2) {
    wind <- wind * log(67.8 * 2 - 5.42) / log(67.8 * wind_height - 5.42)
  }
  
  # Slope of saturation vapor pressure curve (Delta, kPa/°C)
  Delta <- 4098 * (0.6108 * exp((17.27 * t) / (t + 237.3))) / ((t + 237.3)^2)
  
  # Aerodynamic resistance
  ra <- 208 / wind
  
  # PET (kg/m²/s)
  PET_kg_m2_s <- (Delta * rad + rho_air_kg * cp_kg * VPD_Pa / ra) /
    (lambda_kg * (Delta + gamma_pa / 1000 * (1 + 0.34 * wind)))  # convert gamma to kPa/°C
  
  # PET (mol/m²/s)
  PET_mol_m2_s <- PET_kg_m2_s / Mw
  
  return(list(
    PET_mol_m2_s = PET_mol_m2_s,
    gs_mol_m2_s = gs_mol_m2_s
  ))
}

#conversion for transpiration
predict_transpiration_pm_gc <- function(gc, t, vpd, rad, ws, elv = 0, wind_height = 2) {
  # ---- CONSTANTS ----
  R        <- 8.314        # J/mol/K
  cp       <- 1013         # J/kg/K
  r_air    <- 287.05       # J/kg/K
  r_wat    <- 461.5        # J/kg/K
  p0       <- 101325       # Pa
  g        <- 9.81         # m/s²
  
  # ---- Convert temperature to Kelvin ----
  tK <- t + 273.15
  
  # ---- Adjust pressure at elevation ----
  tp <- tK + 0.005 * elv / 2
  P <- p0 * exp(-g * elv / (r_air * tp))  # Pa
  
  # ---- Latent heat of vaporization (temperature-dependent) ----
  lambda_g <- 2500.8 - 2.36 * t + 0.0016 * t^2 - 0.00006 * t^3  # J/g
  lambda <- lambda_g * 1000  # J/kg
  
  # ---- Psychrometric constant (Pa/K) ----
  gamma <- (P / 100) * cp / lambda_g / (r_air / r_wat)
  
  # ---- Convert VPD to Pa ----
  VPD_Pa <- vpd * 1000
  
  # ---- Slope of saturation vapor pressure curve (kPa/°C) ----
  Delta <- 4098 * (0.6108 * exp((17.27 * t) / (t + 237.3))) / ((t + 237.3)^2)
  
  # ---- Air density (kg/m³) ----
  rho_air <- P / (r_air * tK)
  
  # ---- Adjust wind speed to 2 m if needed ----
  if (wind_height != 2) {
    ws <- ws * log(67.8 * 2 - 5.42) / log(67.8 * wind_height - 5.42)
  }
  
  # ---- Aerodynamic resistance ----
  ra <- 208 / ws  # s/m
  
  # ---- Stomatal resistance ----
  rs <- (1 / gc) * (R * tK) / P  # s/m
  
  # ---- Penman–Monteith equation ----
  numerator <- Delta * rad + rho_air * cp * VPD_Pa / ra
  denominator <- lambda * (Delta * 1000 + gamma * (1 + rs / ra))  # Delta converted to Pa
  
  E_kg_m2_s <- numerator / denominator
  
  # ---- Convert to g/m²/h ----
  E_g_m2_h <- E_kg_m2_s * 1000 * 3600
  
  return(E_g_m2_h)
}