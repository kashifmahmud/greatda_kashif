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
met$DateTime_hr <- as.POSIXct(met$DateTime_hr,format="%Y-%m-%d %T",tz="UTC")

gpp_na<-subset(met,met$Date %in% c(as.Date("2016-1-18"),as.Date("2016-1-19")))
gpp_20<-subset(met,met$Date %in% as.Date("2016-1-19"))
gpp_20$Date<-as.Date("2016-01-20")
gpp_20$DateTime_hr<-rep(seq.POSIXt(as.POSIXct("2016-1-20 00:00:00",tz="UTC"), as.POSIXct("2016-01-20 23:45:00",tz="UTC"), by = "15 min"),each=6)


gpp_21<-subset(met,met$Date %in% as.Date("2016-1-22"))
gpp_21$Date<-as.Date("2016-01-21")
gpp_21$DateTime_hr<-rep(seq.POSIXt(as.POSIXct("2016-1-21 00:00:00",tz="UTC"), as.POSIXct("2016-01-21 23:45:00",tz="UTC"), by = "15 min"),each=6)

# gpp_na$Date<-ifelse(gpp_na$Date=="2016-1-18","2016-1-20","2016-1-21")  

met_na<-rbind(gpp_20,gpp_21)
met<-rbind(met,met_na) 

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
                                   Rd=0))


photom <- photo_gr[,c(1:6, 8)]
photom_15min <- cbind(photom,metwithphy) # add met data 


# convert GPP to g (Carbon) m-2 15 min

photom_15min$Date <- as.Date(photom_15min$Date)
photom_15min$gCarbon <- with(photom_15min, ALEAF*15*60*10^-6*12.0107) #g Carbon m-2 15 min
photom_15min$gCarbon_la <- with(photom_15min, gCarbon*LA) #total carbon gain (g Carbon plant-1 15 min)


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# calculate stem and root respiration (units gC gC 15 min-1)

photom_15min$Rdark_stem<-with(photom_15min, R_stem*2.1^((Tair-25)/10)) # units: gC gC (biomass) s-1
photom_15min$Rdark_root<- with(photom_15min, R_root*2.1^((Tair-25)/10))# units: gC gC (biomass) s-1


# unit conversion of stem and root respiration rates
photom_15min$Rdark_stem<-with(photom_15min,Rdark_stem*15*60)  # unit conversion to gC plant-1 15 min-1
photom_15min$R_stem_se<-with(photom_15min,R_stem_se*15*60)# units gC gC-1 per 15 min


photom_15min$Rdark_root<- with(photom_15min,Rdark_root*15*60)# unit conversion to gC plant-1 15 min-1
photom_15min$R_root_se<-with(photom_15min,R_root_se*15*60)# units gC gC-1 per 15 min

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# calculate leaf respiration rates (units gC gC 15 min-1)

photom_15min$Rdark_leaf<-ifelse(photom_15min$period=="Day",with(photom_15min,Rd.mean*0.7^((Tair-25)/10)*(LA/Leafmass)*10^-6*12.0107*15*60),
                                with(photom_15min,R_leaf*2.1^((Tair-25)/10))*15*60)


# photom_15min$Rdark_leaf<-ifelse(photom_15min$period=="Day",with(photom_15min,R_leaf*2.1^((Tair-25)/10))*15*60*0.7,
#                                 with(photom_15min,R_leaf*2.1^((Tair-25)/10))*15*60)

photom_15min$R_leaf_se<-with(photom_15min,R_leaf_se*15*60)# units gC gC-1 per 15 min

photom_15min$SLA<-with(photom_15min,LA/Leafmass)
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

photom_15min<-photom_15min[c(8,9,10,12,14,15,24,25,27,28,35:50,52:55,57:61)]
names(photom_15min)[c(31:34)]<-c("GPP","R_stem","R_root","R_leaf")



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

# # add DA parameters
# # has to run DA first
# 
# paramslist <- vector(mode = "list", length = nlevels(treat.group))
# for (i in 1:nlevels(treat.group)) {
#   paramslist[[i]] <- data.frame(result[[i]][[2]])
# }
# 
# params = do.call("rbind", paramslist)
# params<-params[c(1:3,5)]
# names(params)[4]<-"Room"
# 
# params_new<-tidyr::spread(params, variable, Parameter, convert = FALSE)
# 
# metwithphy_n<-merge(photom_15min,params_new,by=c("Date","Room"))
# 
# write.csv(params_new,"output/DA_parameters_final.csv")
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

# add DA fitted data
# has to run DA first

# plot modelled mass data
# 
# 
# paramslist <- vector(mode = "list", length = nlevels(treat.group))
# for (i in 1:nlevels(treat.group)) {
#   paramslist[[i]] <- data.frame(result[[i]][[4]])
# }
# 
# meas_dat <- do.call("rbind", paramslist)
# meas_dat<-subset(meas_dat,!is.na(value))
# 
# 
# data_new<-tidyr::spread(meas_dat, variable, value, convert = FALSE)
# names(data_new)[2]<-"Room"
# metwithphy_nd<-merge(metwithphy_n,data_new[c(1,2,4:7)],by=c("Date","Room"))
# 
# write.csv(data_new,"output/DA_data_final.csv")

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

#--- read self shading factor (used mean sigma for each room)
sigma_s<-read.csv("Parameters/leaf_area_vs_self_shading.csv")

metwithphy_new<-merge(photom_15min,sigma_s[c(1,2,4,11)],by="Room")
metwithphy_new$DateTime_hr <- as.POSIXct(metwithphy_new$DateTime_hr,format="%Y-%m-%d %T",tz="UTC")
metwithphy_new$Date <- date(metwithphy_new$DateTime_hr)


# add seedling mass data

metwithphy_new<-merge(metwithphy_new,seedling_mass[c(1:2,5,7,9)],by=c("Date","Room"))

# names(metwithphy_new)[c(41,42,43)]<-c("Leafmass","Stemmass","Rootmass")

metwithphy_gpp <- metwithphy_new[order(metwithphy_new$DateTime_hr), ]


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

write.csv(metwithphy_gpp,"Parameters/data_for_attribution_analysis_V2.csv")


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# 
# room_1<-subset(metwithphy_gpp,metwithphy_gpp$Room==1 & metwithphy_gpp$Date==as.Date("2016-01-08"))
# 
# sum(room_1$R_leaf)
# 
# #-----
# gpp_att<-photom_15min