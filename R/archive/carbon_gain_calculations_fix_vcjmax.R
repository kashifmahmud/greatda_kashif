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
source("R/loadLibraries.R")
source("R/functions_for_analysis.R")
source("R/generic_functions.R")
source("R/GREAT_functions.R")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#--- read met data and fitted parameters

#--- read GREAT climate data 15 min averages

met<-read.csv("Parameters/s39_climate_15_min_averages.csv")
met$Date<-as.Date(met$Date)

#-- calculate Age of seedlings
met$Age_of_plant<-with(met,Date-met$Date[[1]]) # dates starts from the begining of temperature treatments
met$Age_of_plant<-as.numeric(met$Age_of_plant)




#--- read g1

#short_g1<-read.csv("Parameters/great_g1_short_term.csv") # short-term tempearture response

short_g1<-read.csv("Parameters/great_g1_short_term.csv") # short-term tempearture response



#--- read alpha

alph<-read.csv("Parameters/great_alpha_long_term_nls_method.csv") # alpha from nls method


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# make dataframe with g1 and alpha

phy.data<-merge(short_g1[c(1:3)],alph,by="Room") 



# #add fixed biochemical parameters (Ea and deltaS of Vcmax and Jmax)
# 
# bioch_fixed<-read.csv("Parameters/great_aci_t_response_fixed.csv")
# # bioch_fixed$Date<-as.Date(bioch_fixed$Date)
# 
# phy.data<- merge(phy.data,bioch_fixed,by="Room",all=T) 
# 
# 
# # write.csv(phy.data,file = "Parameters/fixed_physiology_parameters.csv",row.names=FALSE,sep=",")

#--- add physiological parameters to med data

metwithphy<-merge(met,phy.data,by="Room")



#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#--- read vcmax, jmax and temperature response parameters (fitted to GREAT ACi data)

bioch<-read.csv("yplantupscale/data/great_aci_t_response.csv")[c(1:4,9:11)]
# bioch$Date<-as.Date(bioch$Date)

metwithphy<-merge(metwithphy,bioch,by=c("Room"),all=T)

# # set Vcmax and Jmax from 7/1/2016 to 2/2/2016 (similar to values on 3/2/2016)
# metwithphy$Vcmax25<-ifelse(metwithphy$Date< as.Date("2016-02-03"),subset(bioch,bioch$Date=="2016-02-03")$Vcmax[[1]],metwithphy$Vcmax25)
# metwithphy$Jmax25<-ifelse(metwithphy$Date< as.Date("2016-02-03"),subset(bioch,bioch$Date=="2016-02-03")$Jmax[[1]],metwithphy$Jmax25)


# # set Vcmax and Jmax from 17/2/2016 onwards (similar to values at growth temperatures on 16/2/2016)
# for(i in 1:length(unique(metwithphy$Room))){
#   
#   metwithphy$Vcmax25[which(metwithphy$Date>as.Date("2016-02-16") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-16" & bioch$Room==i)$Vcmax25
#   metwithphy$Jmax25[which(metwithphy$Date>as.Date("2016-02-16") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-16" & bioch$Room==i)$Jmax25
#   
# }
# 
# # set Vcmax and Jmax from the begining to 02/2/2016  (similar to values at growth temperatures on 02/2/2016)
# for(i in 1:length(unique(metwithphy$Room))){
#   
#   metwithphy$Vcmax25[which(metwithphy$Date<as.Date("2016-02-03") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-03" & bioch$Room==i)$Vcmax25
#   metwithphy$Jmax25[which(metwithphy$Date<as.Date("2016-02-03") & metwithphy$Room==i)]<-subset(bioch,bioch$Date=="2016-02-03" & bioch$Room==i)$Jmax25
#   
# }

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

# -- plot estimated Vcmax25 and Jmax25
# 
# COL <- rev(brewer.pal(6,"Spectral"))
# windows(80,100);
# par(cex.lab=1.5,mar=c(0,5,0.5,0.5),oma=c(4,4,0,0),cex.axis=1.5,las=1,mfrow=c(2,1))
# 
# plotBy(Vcmax25~Date|Room,type="l",data=metwithphy,legend=F,lwd=3,ylab="",xlab="")
# # with(bioch,points(Date,Vcmax25,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
# # axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
# #legend("topright",letters[1],bty="n",cex=1.2)
# 
# title(ylab=expression(V[cmax25]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=3.5)
# abline(v=as.Date("2016-02-16"),lty=3)
# abline(v=as.Date("2016-02-03"),lty=3)
# legend("bottomleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
# 
# 
# plotBy(Jmax25~Date|Room,type="l",data=metwithphy,legend=F,lwd=3,ylab="",xlab="")
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
# # with(bioch,points(Date,Jmax25,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
# 
# #legend("topright",letters[2],bty="n",cex=1.2)
# 
# # legend("topright",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
# 
# 
# title(ylab=expression(J[max25]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=3.5)
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3)
# abline(v=as.Date("2016-02-16"),lty=3)
# abline(v=as.Date("2016-02-03"),lty=3)
# 
# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#set up photosynthesis model 

#Day respiration: Basal rate (at 22C) and Q10 were assumed to be simila to dark respiration rate
#Values from Drake et al 2017 GCB; averaged across provenances (table 2)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------



#make negatine PAR values (recorded in night) to zero.

metwithphy$PAR<-ifelse(metwithphy$PAR<0,0,metwithphy$PAR)
metwithphy<-subset(metwithphy,!is.na(PAR))  

metwithphy <- metwithphy[order(metwithphy$DateTime_hr) , ]
# run room 1 


photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha,
                                   theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25
                                   , Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                   delsC = delsV*10^3,delsJ = 0.6312849*10^3,EdVC = 2e+05,EdVJ = 2e+05,EaV = EaV*10^3,
                                   EaJ = 42.74635*10^3,Rdayfrac=.7))


# photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,
#                                    theta = 0.85,  Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                   ))
# 
photom <- photo_gr[,c(1:6, 8)]
photom_15min <- cbind(photom,metwithphy) # add met data 


# convert GPP to g (Carbon) m-1s-1

photom_15min$Date <- as.Date(photom_15min$Date)
photom_15min$gCarbon <- with(photom_15min, ALEAF*15*60*10^-6*12.0107) #g Carbon m-1 s-1

# dCarbon <- dplyr::summarize(group_by(photom_15min,Date,Room),
#                             gCarbon = sum(gCarbon),
#                             ALEAF=mean(ALEAF),
#                             maxALEAF= max(ALEAF),
#                             Rd=mean(Rd),
#                             Tair=mean(Tair),
#                             VPD=mean(VPD),
#                             PAR_mean=mean(PAR),
#                             PAR_max=max(PAR),
#                             RH=mean(RH)
#                             )


# get some summary stats and daily carbon gain
# Aleaf_summary <- doBy::summaryBy(ALEAF+Rd ~ Date+Room, data=photom_15min, FUN=c(mean,min,max))



#-- get daily total carbon gain for each temperature treatment

dCarbon <- doBy::summaryBy(gCarbon ~ Date+Room, data=photom_15min, FUN=sum, keep.names=TRUE ) #g Carbon m-1 day-1
dCarbon <-subset(dCarbon,Date>as.Date("2016-01-07") & Date <= as.Date("2016-03-02"))

# dCarbon<-merge(Aleaf_summary,carb_sum,by=c("Date","Room"))

#names(dCarbon)[3] <- "carbon_day"
#get daily averages of met data


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#--- add daily met data to daily carbongain data

met_daily<-read.csv("Parameters/s39_climate_daily_averages.csv")

carb_with_met<-merge(dCarbon,met_daily,by=c("Date","Room"))

carb_with_met$PARFac<-cut(carb_with_met$PARsum_mol,breaks=c(0,50,100,150),labels=F)


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#-- add daily total stem/root/leaf mass data (units g(C))



srl_mass<-read.csv("Parameters/great_leaf_area_predicted.csv")
names(srl_mass)[3]<-"LA"
srl_mass$Date<-as.Date(srl_mass$Date,format="%Y-%m-%d")
carb_with_met_la<-merge(carb_with_met,srl_mass,by=c("Date","Room"))


# la_daily<-read.csv("Parameters/leaf_area_with_self_shading.csv")
# la_daily$Date<-as.Date(as.character(la_daily$Date),format="%d/%m/%Y")
# 
# 
# carb_with_met_la<-merge(carb_with_met,la_daily,by=c("Date","Room"))

#--- get daily total GPP (gC plant-1 day-1 )
carb_with_met_la$cday<-with(carb_with_met_la,(gCarbon*LA)) # total GPP per day not corrected for self shading


#--- multiply GPP by self shading factor

#--- read self shading factor (used mean sigma for each room)

sigma_s<-read.csv("Parameters/leaf_area_vs_self_shading.csv")
carb_with_met_la<-merge(carb_with_met_la,sigma_s[c(1,11)],by="Room")

carb_with_met_la$GPP<-with(carb_with_met_la,cday*self_s.mean)

#carb_with_met_la$GPP<-with(carb_with_met_la,gCarbon_LA*sigma)  # total GPP per day with self shading



#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#-- add leaf, stem and root respiration

Rdat_mean<-read.csv("Parameters/great_Resp_leaf_shoot_root_fh.csv")

carb_with_met_la<-merge(carb_with_met_la,Rdat_mean,by="Room")

# 
# carb_with_met_la<-carb_with_met_la[c(1:2,12,16,17,20:28)]
# names(carb_with_met_la)[c(5,7)]<-c("LA","dailyCarbon")


#------
#------
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#-- Gap filling in GPP data
#-- 2016-1-20 and 2016-1-21 are the rotation dates (rooms changed). these two days, met data cannot be used  
#-- assumed 2016-1-19 GPP values for these two days

gpp_na<-subset(carb_with_met_la,carb_with_met_la$Date %in% c(as.Date("2016-1-18"),as.Date("2016-1-19")))
gpp_na$Date<-ifelse(gpp_na$Date=="2016-1-18","2016-1-20","2016-1-21")  


carb_with_met_la<-rbind(carb_with_met_la,gpp_na) 



#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

write.csv(carb_with_met_la,file = "Parameters/great_daily_carbon_gain_LA.csv",row.names=FALSE,sep=",")

#-------------------------------------------------
#-------------------------------------------------


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

#-- calculate daily maintenance respiration of leaf/stem and root 
#-- respiration numbers are in gC gC-1(biomass)
#-- R (gC gC(biomass))*biomass


carb_with_met_la$Rm_leaf<-with(carb_with_met_la,R_leaf*Leafmass) # leaf mass (g) converted to g(C) 
carb_with_met_la$Rm_stem<-with(carb_with_met_la,R_stem*Stemmass)# stem mass (g) converted to g(C)
carb_with_met_la$Rm_root<-with(carb_with_met_la,R_root*Rootmass)# root mass (g) converted to g(C)


#-- get daily total Rm (gC)

carb_with_met_la$Rm<-with(carb_with_met_la,Rm_leaf+Rm_stem+Rm_root)


#-- get total GPP (sum over total duration) and Rm (sum over total duration) for each temperature treatment
#-- units in g(C)
dat5<-summaryBy(GPP+Rm~Room,data=carb_with_met_la,FUN =sum)
dat5$Tair.mean<-c(18,21.5,25,28.5,32.5,35.5)


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

#-- get final mass data and add to GPP dataframe

f_mass<-subset(srl_mass,srl_mass$Date==as.Date("2016-02-29"))

f_mass<-transform(f_mass, sum=rowSums(f_mass[c(4:6)])) # get total plant weight in each treatments (g(C))
names(f_mass)[7]<-"totdm_f"


dat5<-merge(dat5,f_mass[c(2:7)],by="Room")

#-- add final mass given in John's paper
dat6<-read.csv("Parameters/great_final_mass.csv") # units (g)

dat5<-merge(dat5,dat6[c(1:2)],by="Room")

# add initial mass to the dataframe (mass data on 8/1/2016)

dat5$totdm_i<-with(srl_mass[1,],sum(Leafmass,Stemmass,Rootmass)) #units in g(C) and assume similar initial mass for all treatments



#-- get biomass change during experimant period
dat5$fCarbon_calc<-((dat5$totdm_f-dat5$totdm_i)*1.3)+dat5$Rm.sum  #sum of GPP= sum of Rm+(1.3*Biomass change)

dat5$storage<-with(dat5,GPP.sum-fCarbon_calc)

write.csv(dat5,file = "Parameters/great_final_carbon_estimates.csv",row.names=FALSE,sep=",")


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#-- plot carbon balance

windows(40,40);
par(cex.lab=1.5,mar=c(5,5,1,0.5),oma=c(0,0,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))# 

#-- modelled GPP (C)
with(dat5,plot(GPP.sum~Tair.mean,pch=17,cex=2,ylim=c(0,25),ylab="",xlab="",col="red"))
# 
#-- add final mass (C) corrected for Rg (final mass*1.3)
with(dat5,points(fCarbon_calc~Tair.mean,pch=16,cex=2,ylab="",xlab=""))
# adderrorbars(x=dat5$Tair.mean,y=(dat5$totdm.mean*0.47*1.3)+dat5$Rm.sum,SE=dat3$totdm.standard.error*.47,direction="updown")

title(ylab=expression(C~(g)),cex.lab=2,line=2)
title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=2)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#--- plot storage C
windows(40,40);
par(cex.lab=1.5,mar=c(5,5,1,0.5),oma=c(0,0,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))# 

#-- modelled GPP (C)
with(dat5,plot(storage~Tair.mean,pch=17,cex=2,ylim=c(-2,10),ylab="",xlab="",col="red"))
# 
#-- add final mass (C) corrected for Rg (final mass*1.3)
with(dat5,points(fCarbon_calc~Tair.mean,pch=16,cex=2,ylab="",xlab=""))
# adderrorbars(x=dat5$Tair.mean,y=(dat5$totdm.mean*0.47*1.3)+dat5$Rm.sum,SE=dat3$totdm.standard.error*.47,direction="updown")

title(ylab=expression(C~(g)),cex.lab=2,line=2)
title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=2)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


# 
# 
# 
# 
# 