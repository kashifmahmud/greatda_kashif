
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------


## this script process met data and parameters for the attribution analysis.

#' add met data, seedling mass data and physiological parameters to a single file
#' Nore: all respiration data are in gC gC-1 s-1 

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
rday_aci$Rd.mean <- with(rday_aci, Rd.mean*10^-6*12.0107) # unit conversion to gC m-2 s-1

srl_mass<-read.csv("Parameters/modelled_data.csv")
srl_mass$SLA<-with(srl_mass,Leafarea/Leafmass)

sla_mean<-summaryBy(SLA~Room,FUN=mean,data=srl_mass,keep.names=T)

rday_aci<-merge(rday_aci,sla_mean,by="Room")

rday_aci$Rd.mean <- with(rday_aci, Rd.mean*SLA) # unit conversion to gC gC-1 s-1

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

# srl_mass<-read.csv("Parameters/modelled_data.csv")
# srl_mass$SLA<-with(srl_mass,Leafarea/Leafmass)


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

# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------


#- calculate leaf respiration rates (Rm Leaf) (mass basis)
#- assume Rday=0.7*Rdark 
#- units in gC gC-1 s-1
#- use Rdark measured at 25C with Q10=2 to calculate Rdark at growth temperatures in 15 min intervals
# Rday not included in GPP calculations
q10<-2.13
q10_d<-0.7
metwithphy$Rdark_leaf<-ifelse(metwithphy$period=="Day",with(metwithphy,Rd.mean*q10_d^((Tair-25)/10))*0.7, 
                         with(metwithphy,R_leaf*q10^((Tair-25)/10))) # units gC gC-1 s-1

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# calculate stem and root respiration

metwithphy$Rdark_stem<-with(metwithphy, R_stem*2^((Tair-25)/10)) # units: gC gC-1 s-1
metwithphy$Rdark_root<- with(metwithphy, R_root*2^((Tair-25)/10))# units: gC gC-1 s-1


names(metwithphy)


metwithphy<-metwithphy[c(1:3,5,7:8,17,20:44,46:52)]
names(metwithphy)[c(37:39)]<-c("R_leaf","R_stem","R_root")
write.csv(metwithphy,"Parameters/data_for_attribution_analysis.csv")


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# set up photosynthesis model 

#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------


# get gross GPP, Rday set to zero

photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha,
                                   theta = theta, Jmax =Jmax25, Vcmax =Vcmax25,delsC=delsV*10^3,delsJ=delsJ*10^3,EaV=EaV*10^3,EaJ=EaJ*10^3,
                                   Rd=0))
