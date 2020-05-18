#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#--- Daily carbon gain calculations for GREAT experiments
#--- This script read the  pre-estimated physiological parameters and met data for GREAT experiment
#--- and used to model leaf photosynthesis (~carbon gain)


# source("R/functions_for_analysis.R")
# source("R/generic_functions.R")
#devtools::install_bitbucket("remkoduursma/plantecophys") 
#library(plantecophys)
#library(dplyr)
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#--- get required functions and packages

path<-getwd()
# source("R/loadLibraries.R")
# library(plantecophys)
# library(doBy)
# library(plotBy)
# library(RColorBrewer)
# library(kableExtra)

source("R/functions/loadLibraries.R")
source("R/functions/functions_for_analysis.R")
source("R/functions/generic_functions.R")
source("R/functions/GREAT_functions.R")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------


#--- read met data and fitted parameters

#--- read GREAT met data 15 min averages

met<-read.csv("Data/MET_DATA/s39_climate_15_min_averages.csv")
met$Date<-as.Date(met$Date)


# met$Tair<-ifelse(met$Room %in% c(4,5),met$Tair-2,met$Tair)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# read leaf, stem and root Rdark values for each room  (values at 25C measured at final harvest)

rday<-read.csv("Parameters/great_Resp_leaf_shoot_root_25C.csv")

met<-merge(met,rday,by="Room")

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# add Rday estimated form ACi curves

rday_aci<-read.csv("Parameters/rday_aci_fits.csv")

met<-merge(met,rday_aci,by="Room")



#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#-- add daily stem/root/leaf mass data (units g(C))
#-- predicted leaf/stem/root mass data from allometric models. smooth curve fitted to get the time course

# srl_mass<-read.csv("Parameters/great_leaf_area_predicted.csv")
# names(srl_mass)[3]<-"LA"
# srl_mass$Date<-as.Date(srl_mass$Date,format="%Y-%m-%d")
# met<-merge(met,srl_mass,by=c("Date","Room"))
# #-- add daily stem/root/leaf mass data (units g(C))
# #-- predicted leaf/stem/root mass data from allometric models. smooth curve fitted to get the time course
# 
# # srl_mass<-read.csv("Parameters/great_leaf_area_predicted.csv")

srl_mass<-read.csv("Parameters/modelled_data.csv")
# names(srl_mass)[3]<-"LA"
srl_mass$Date<-as.Date(srl_mass$Date,format="%Y-%m-%d")

datesout <- seq(as.Date("2016/1/8"), as.Date("2016/2/29"), by="1 day")

interp_var <- function(d, datesout,var){
  
  # sp <- spline(x=d$Date, y=d[[var]], xout=datesout) #smooth curve
  sp <- approx(x=d$Date, y=d[[var]], xout=datesout) # linear interpollation
  
  out<-data.frame(Date=datesout, Room=unique(d$Room), var=sp$y)
  return(out)
}


leaf_area <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Leafarea"))
names(leaf_area)[3]<-"LA"

leaf_area.se <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Leafarea_SE"))
names(leaf_area.se)[3]<-"LA_SE"

Leafmass <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Leafmass"))
names(Leafmass)[3]<-"Leafmass"

Leafmass.se <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Leafmass_SE"))
names(Leafmass.se)[3]<-"Leafmass_SE"


Stemmass <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Stemmass"))
names(Stemmass)[3]<-"Stemmass"

Stemmass.se <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Stemmass_SE"))
names(Stemmass.se)[3]<-"Stemmass_SE"


Rootmass <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Rootmass"))
names(Rootmass)[3]<-"Rootmass"

Rootmass.se <- do.call(rbind, lapply(split(srl_mass, srl_mass$Room), interp_var, datesout=datesout,var="Rootmass_SE"))
names(Rootmass.se)[3]<-"Rootmass_SE"

seedling_mass<-cbind(leaf_area,leaf_area.se[3],Leafmass[3],Leafmass.se[3],Stemmass[3],Stemmass.se[3],Rootmass[3],Rootmass.se[3])


met<-merge(met,seedling_mass,by=c("Date","Room"))

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# read dataframe with g1,alpha,theta and temperature response parameters Ea, dS 

phy.data<-read.csv("Parameters/great_alpha_g1_ea_dels.csv")
# phy.data$alpha<-ifelse(phy.data$Room==4,0.3522958,phy.data$alpha)

phy.data$Date<-NULL

#--- add physiological parameters to med data

metwithphy<-merge(met,phy.data,by="Room")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#--- read vcmax and jmax at 25C
#--- this data contains daily Vcmax and Jmax values from 2016-02-03 to 2016-02-26 

bioch<-read.csv("Parameters/great_vcmax_jmax_25C.csv")
bioch$Date<-as.Date(bioch$Date)

metwithphy<-merge(metwithphy,bioch,by=c("Date","Room"),all=T)


# set Vcmax and Jmax from 26/2/2016 onwards (similar to values at growth temperatures on 26/2/2016)
for(i in 1:length(unique(metwithphy$Room))){
  
  metwithphy$Vcmax25[which(metwithphy$Date>as.Date("2016-02-26") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-26" & bioch$Room==i)$Vcmax25
  metwithphy$Jmax25[which(metwithphy$Date>as.Date("2016-02-26") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-26" & bioch$Room==i)$Jmax25
  
}

# set Vcmax and Jmax from the begining to 02/2/2016  (similar to values at growth temperatures on 02/2/2016)

for(i in 1:length(unique(metwithphy$Room))){
  
  metwithphy$Vcmax25[which(metwithphy$Date<as.Date("2016-02-03") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-03" & bioch$Room==i)$Vcmax25
  metwithphy$Jmax25[which(metwithphy$Date<as.Date("2016-02-03") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-03" & bioch$Room==i)$Jmax25
  
}

# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------


#make negatine PAR values (recorded in night) to zero. No need but sometimes it gave an error

metwithphy$PAR<-ifelse(metwithphy$PAR<0,0,metwithphy$PAR) 
metwithphy<-subset(metwithphy,!is.na(PAR))  
metwithphy <- metwithphy[order(metwithphy$DateTime_hr), ]

metwithphy$period <- ifelse(metwithphy$PAR>2,"Day","Night")



# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------

# set up photosynthesis model 

startdate<-as.Date("2016-01-07")
enddate<-as.Date("2016-02-24")


# get gross GPP, Rday set to 70% Rdark

# photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha,
#                                    theta = theta, Jmax =Jmax25, Vcmax =Vcmax25,delsC=delsV*10^3,delsJ=delsJ*10^3,EaV=EaV*10^3,EaJ=EaJ*10^3,
#                                    Rd0=ifelse(metwithphy$period=="Day",metwithphy$Rd.mean,metwithphy$R_leaf_a),
#                                    Q10=ifelse(metwithphy$period=="Day",0.7,2.1)))


photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha=alpha,
                                    theta = theta, Jmax =Vcmax25*JVr, Vcmax =Vcmax25,delsC=delsV*10^3,delsJ=delsJ*10^3,
                                   EaV=EaV*10^3,EaJ=EaJ*10^3,
                                   Rd0=ifelse(metwithphy$period=="Day",metwithphy$Rd.mean,metwithphy$R_leaf_a),
                                   Q10=ifelse(metwithphy$period=="Day",0.7,2.1)))


photom <- photo_gr[,c(1:6, 8)]
photom_15min <- cbind(photom,metwithphy) # add met data 


# convert GPP to g (Carbon) m-2 15 min

photom_15min$Date <- as.Date(photom_15min$Date)
photom_15min$gCarbon <- with(photom_15min, ALEAF*15*60*10^-6*12.0107) #g Carbon m-2 15 min
photom_15min$gCarbon_la <- with(photom_15min, gCarbon*LA) #total carbon gain (g Carbon plant-1 15 min)


#- convert Rd to g (Carbon) m-2 15 min
photom_15min$gRd <- with(photom_15min, Rd*15*60*10^-6*12.0107) #g Carbon m-2 15 min
photom_15min$gRd_la <- with(photom_15min, gRd*LA) #total Rd (g Carbon plant-1 15 min)



# photom_15min$period <- ifelse(photom_15min$PAR>2,"Day","Night")

# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------


#- calculate leaf respiration rates (Rm Leaf) (mass basis)
#- assume Rday=0.7*Rdark 
#- units in gC gC-1 s-1
#- use Rdark measured at 25C with Q10=2 to calculate Rdark at growth temperatures in 15 min intervals
# Rday not included in GPP calculations
# q10<-2.13
# photom_15min$Rdark<-ifelse(photom_15min$period=="Day",with(photom_15min,R_leaf*q10^((Tair-25)/10))*0.7, with(photom_15min,R_leaf*q10^((Tair-25)/10))) # units gC gC-1 s-1
# photom_15min$R_leaf_se<-with(photom_15min,R_leaf_se*15*60)# units gC gC-1 per 15 min
# 
# 
# # photom_15min$Rdark<-photom_15min$Rdark*10^-9*15*60*12.0107/.48 # unit conversion to gC gC-1(leaf mass) 15 min
# photom_15min$Rdark<-with(photom_15min,Rdark*15*60*Leafmass) # total leaf respiration per plant per 15 min (gC)
# 


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# calculate stem and root respiration

photom_15min$Rdark_stem<-with(photom_15min, R_stem*2.1^((Tair-25)/10)) # units: gC gC (biomass) s-1
photom_15min$Rdark_root<- with(photom_15min, R_root*2.1^((Tair-25)/10))# units: gC gC (biomass) s-1


# unit conversion of stem and root respiration rates
photom_15min$Rdark_stem<-with(photom_15min,Rdark_stem*15*60*Stemmass)  # unit conversion to gC plant-1 15 min-1
photom_15min$R_stem_se<-with(photom_15min,R_stem_se*15*60)# units gC gC-1 per 15 min


photom_15min$Rdark_root<- with(photom_15min,Rdark_root*15*60*Rootmass )# unit conversion to gC plant-1 15 min-1
photom_15min$R_root_se<-with(photom_15min,R_root_se*15*60)# units gC gC-1 per 15 min

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------


#-- get daily total carbon gain for each temperature treatment
#-- both GPP rate and total C gain given here are net of leaf respiration

dCarbon <- doBy::summaryBy(gCarbon+gCarbon_la~Date+Room, data=photom_15min, FUN=sum, keep.names=TRUE ) #gC day-1

# daily total stem and root respiration
dCarbon_resp <- doBy::summaryBy(Rdark_stem+Rdark_root~Date+Room, data=photom_15min, FUN=c(sum), keep.names=TRUE ) #gC day-1


# daily total leaf respiration

# during night
dCarbon_leaf_resp_n <- (doBy::summaryBy(gRd_la~Date+Room, data=subset(photom_15min,photom_15min$period=="Night"), FUN=c(sum), keep.names=TRUE )) #gC day-1
names(dCarbon_leaf_resp_n)[3]<-"R_leaf_n"
# dCarbon_leaf_resp_n$R_leaf_n<-with(dCarbon_leaf_resp_n,abs(R_leaf_n))

# during day

dCarbon_leaf_resp_d <- (doBy::summaryBy(gRd_la~Date+Room, data=subset(photom_15min,photom_15min$period=="Day"), FUN=c(sum), keep.names=TRUE )) #gC day-1
names(dCarbon_leaf_resp_d)[3]<-"R_leaf_d"

dCarbon_leaf_resp<-merge(dCarbon_leaf_resp_n,dCarbon_leaf_resp_d,by=c("Date","Room"))
dCarbon_leaf_resp$Rleaf<-with(dCarbon_leaf_resp,R_leaf_n+R_leaf_d) # daily total leaf respiration




dCarbon_resp_rs<-merge(dCarbon,dCarbon_resp,by=c("Date","Room"))
dCarbon<-merge(dCarbon_resp_rs,dCarbon_leaf_resp,by=c("Date","Room"))

# names(dCarbon_resp)[c(3:8)]<-c("R_leaf","R_stem","R_root","R_leaf_se","R_stem_se","R_root_se")
names(dCarbon)[c(5,6,9)]<-c("R_stem","R_root","R_leaf")

# dCarbon$R_leaf<-with(dCarbon,R_leaf*-1)

# dCarbon<-merge(dCarbon,dCarbon_resp,by=c("Date","Room"))
dCarbon <-subset(dCarbon,Date >startdate & Date <= enddate) 

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# get parameter output for attribution analysis

daily_parameters<-doBy::summaryBy(Vcmax25+Jmax25+EaV+delsV+EaJ+delsJ+JVr+alpha+g1+LA ~Date+Room,data=photom_15min, FUN=mean, keep.names=TRUE )
daily_parameters<-merge(daily_parameters,dCarbon,by=c("Date","Room"))

write.csv(daily_parameters,"Parameters/daily_parameters_for_sensitivity.csv",sep=",")


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# calculate standard error of leaf, stem and root respiration
sum_sqrt<-function(x,na.rm=FALSE){
  
  # source: Hogan, 2006 (http://www.met.rdg.ac.uk/~swrhgnrj/combining_errors.pdf)
  sqrt(sum(x^2,na.rm=na.rm))
}

dCarbon_resp_se <- doBy::summaryBy(R_leaf_se+R_stem_se+R_root_se~Date+Room, data=photom_15min, FUN=c(sum_sqrt), keep.names=TRUE) #gC day-1
dCarbon<-merge(dCarbon,dCarbon_resp_se,by=c("Date","Room"))


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#-- Gap filling in GPP data
#-- 2016-1-20 and 2016-1-21 are the rotation dates (rooms changed). these two days, met data cannot be used  
#-- assumed 2016-1-19 GPP values for these two days (PAR data ~matched with missing dates)

gpp_na<-subset(dCarbon,dCarbon$Date %in% c(as.Date("2016-1-18"),as.Date("2016-1-19")))
gpp_20<-subset(dCarbon,dCarbon$Date %in% as.Date("2016-1-19"))
gpp_20$Date<-as.Date("2016-01-20")
gpp_21<-subset(dCarbon,dCarbon$Date %in% as.Date("2016-1-19"))
gpp_21$Date<-as.Date("2016-01-21")

# gpp_na$Date<-ifelse(gpp_na$Date=="2016-1-18","2016-1-20","2016-1-21")  
dCarbon<-rbind(dCarbon,rbind(gpp_20,gpp_21)) 


# #------------------------------------------------------------------------------------------------------------------
# #------------------------------------------------------------------------------------------------------------------


#--- multiply GPP by self shading factor

# add daily leaf area to the dataframe

dCarbon_with_la<-merge(dCarbon,seedling_mass[c(1:3)],by=c("Date","Room"))


#--- read self shading factor (used mean sigma for each room)
sigma_s<-read.csv("Parameters/leaf_area_vs_self_shading.csv")

dCarbon_with_la<-merge(dCarbon_with_la,sigma_s,by="Room")

dCarbon_with_la$GPP<-with(dCarbon_with_la,gCarbon_la*(Intercept+Slope*LA))  # total GPP per day with self shading/sigma as a function of LA
# dCarbon_with_la$GPP<-with(dCarbon_with_la,gCarbon_la*0.9)  # total GPP per day with self shading/sigma fixed to mean across rooms

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# remove unwanted columns

dCarbon_with_la<-dCarbon_with_la[c("Room","Date","GPP","R_leaf","R_stem","R_root","R_leaf_se","R_stem_se","R_root_se"
)]

# add leaf stem and root mass data

dCarbon_with_la<-merge(dCarbon_with_la,seedling_mass,by=c("Date","Room"))
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------


write.csv(dCarbon_with_la,file = "Parameters/great_daily_carbon_gain_LA.csv",row.names=FALSE,sep=",") # this is the input file for data assimilation


# copy to DA repo

write.csv(dCarbon_with_la,file = sprintf("processed_data/great_daily_carbon_gain_LA_%s.csv",enddate),row.names=FALSE,sep=",")



#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------


#-- get daily total Rm (gC)

# dCarbon_with_la$R_leaf<-with(dCarbon_with_la,R_leaf+(R_leaf*.7)) # roughly add day respiration

dCarbon_with_la$Rm<-with(dCarbon_with_la,R_leaf+R_stem+R_root)



dCarbon_with_la$Rm_percent<-with(dCarbon_with_la,Rm/GPP*100)
dCarbon_with_la$Rm_sr<-with(dCarbon_with_la,R_stem+R_root) # total root and stem respiration

# carb_with_met_la$Rm<-with(carb_with_met_la,Rm_leaf+Rm_stem+Rm_root)


#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------

# get toral GPP, Rm

dat5<-summaryBy(GPP+Rm+Rm_sr~Room,data=dCarbon_with_la,FUN=sum)

# correct GPP for leaf respiration
# leaf day respiration included in photosynthesis model
# leaf night respiration was calculated by the photosynthesis model 
# add leaf respiration to GPP for carbon balanc

dat5$GPP.sum<-with(dat5,GPP.sum+(Rm.sum-Rm_sr.sum))
dat5$Tair.mean<-c(18,21.5,25,28.5,32.5,35.5)


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

#-- get final mass data and add to GPP dataframe

f_mass<-subset(seedling_mass,seedling_mass$Date==enddate)
# f_mass<-subset(srl_mass,srl_mass$Date==as.Date("2016-02-24"))
f_mass<-transform(f_mass, sum=rowSums(f_mass[c("Leafmass","Stemmass","Rootmass")])) # get total plant weight in each treatments (g(C))
names(f_mass)[11]<-"totdm_f"


# dat5<-merge(dat5,f_mass[c(2:7)],by="Room")
dat5<-merge(dat5,f_mass,by="Room")

# add initial mass to the dataframe (mass data on 8/1/2016)

dat5$totdm_i<-with(srl_mass[1,],sum(Leafmass,Stemmass,Rootmass)) #units in g(C) and assume similar initial mass for all treatments


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

#-- get total C use 

dat5$fCarbon_calc<-((dat5$totdm_f-dat5$totdm_i)*1.3)+dat5$Rm.sum  #sum of GPP= sum of Rm (gC)+(1.3*Biomass change)gC

# dat5$fCarbon_calc<-((dat5$totdm_f_test-dat5$totdm_i)*1.3)+dat5$Rm.sum  #test with John's final biomass estimates


dat5$Rg<-((dat5$totdm_f-dat5$totdm_i)*.3)  # calculate growth respiration rate

dat5$storage<-with(dat5,GPP.sum-fCarbon_calc) # calculate balance C (storage??)
dat5$st_perc<-dat5$storage/dat5$totdm_f
dat5$Rm_perc<-with(dat5,Rm.sum/GPP.sum*100)



write.csv(dat5,file = "Parameters/great_final_carbon_estimates.csv",row.names=FALSE,sep=",")


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

# get plots

# 1. plot vcmax and Jmax

# source("R/plot_vcmax_and_jmax.R")
# 
# # 2. plot leaf, stem and root Rm
# 
# source("R/plot_leaf_stem_root_respiration.R")
# 
# 
# # 3. plot daily GPP
# 
# source("R/plot_daily_carbon_gain.R")


# 4. plot Carbon balance

source("R/plot_carbon_balance.R")


# 5. plot diurnal variation of daily GPP

# source("R/plot_diurnal_cestimates.R")


