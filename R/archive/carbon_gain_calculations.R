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
library(plantecophys)
library(doBy)
library(plotBy)
library(RColorBrewer)
library(magicaxis)

source("R/functions_for_analysis.R")
source("R/generic_functions.R")
source("R/GREAT_functions.R")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#--- read met data and fitted parameters

#--- read GREAT climate data 15 min averages

met<-read.csv("Data/MET_DATA/s39_climate_15_min_averages.csv")
met$Date<-as.Date(met$Date)

# to set Room temperatures for the mean daily T

dat_summary<-summaryBy(Tair~Date+Room,FUN=c(mean,min,max),data=met)

tair_mean<-summaryBy(Tair~Room,data=met,FUN=mean,keep.names=T)

met$Tair<-NULL
met<-merge(met,tair_mean,by="Room")
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# #-- calculate Age of seedlings
# met$Age_of_plant<-with(met,Date-met$Date[[1]]) # dates starts from the begining of temperature treatments
# met$Age_of_plant<-as.numeric(met$Age_of_plant)

 


# #--- read g1
# 
#  #short_g1<-read.csv("Parameters/great_g1_long_term.csv") # short-term tempearture response
# 
#  short_g1<-read.csv("Parameters/great_g1_long_term.csv") # short-term tempearture response
# 
# 
# 
# #--- read alpha
# 
# alph<-read.csv("Parameters/great_alpha_long_term_nls_method.csv") # alpha from nls method

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# read Rday (corrected for rday frac =0.7)

rday<-read.csv("Parameters/great_Resp_leaf_area_fh.csv")[c(2:6)]

met<-merge(met,rday,by="Room")



#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# read dataframe with g1, alpha

# phy.data<-merge(short_g1[c(1:3)],alph,by="Room") 

phy.data<-read.csv("Parameters/great_alpha_g1_ea_dels.csv")
phy.data$Date<-NULL


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------




#--- add physiological parameters to med data

metwithphy<-merge(met,phy.data,by="Room")


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# 
# 
# #add fixed biochemical parameters (Ea and deltaS of Vcmax and Jmax)
# 
# bioch_fixed<-read.csv("Parameters/great_aci_t_response_fixed.csv")
# # bioch_fixed$Date<-as.Date(bioch_fixed$Date)

# phy.data<- merge(phy.data,bioch_fixed,by="Room",all=T) 


# write.csv(phy.data,file = "Parameters/fixed_physiology_parameters.csv",row.names=FALSE,sep=",")




#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#--- read vcmax and jmax  

bioch<-read.csv("Parameters/great_vcmax_jmax.csv")
bioch$Date<-as.Date(bioch$Date)

metwithphy<-merge(metwithphy,bioch,by=c("Date","Room"),all=T)

# # set Vcmax and Jmax from 7/1/2016 to 2/2/2016 (similar to values on 3/2/2016)
# metwithphy$Vcmax25<-ifelse(metwithphy$Date< as.Date("2016-02-03"),subset(bioch,bioch$Date=="2016-02-03")$Vcmax[[1]],metwithphy$Vcmax25)
# metwithphy$Jmax25<-ifelse(metwithphy$Date< as.Date("2016-02-03"),subset(bioch,bioch$Date=="2016-02-03")$Jmax[[1]],metwithphy$Jmax25)


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

# -- plot estimated Vcmax and Jmax

COL <- rev(brewer.pal(6,"Spectral"))
windows(80,100);
par(cex.lab=1.5,mar=c(0,5,0.5,0.5),oma=c(4,4,0,0),cex.axis=1.5,las=1,mfrow=c(2,1))

plotBy(Vcmax25~Date|Room,type="l",col=COL,data=metwithphy,legend=F,lwd=3,ylab="",xlab="",ylim=c(0,300))
# with(bioch,points(Date,Vcmax25,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
#legend("topright",letters[1],bty="n",cex=1.2)

title(ylab=expression(V[cmax]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=3.5)
abline(v=as.Date("2016-02-26"),lty=3)
abline(v=as.Date("2016-02-03"),lty=3)
# legend("bottomleft",legend=unique(metwithphy$Room),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


plotBy(Jmax25~Date|Room,type="l",col=COL,data=metwithphy,legend=F,lwd=3,ylab="",xlab="",ylim=c(0,300))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
# with(bioch,points(Date,Jmax25,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))

#legend("topright",letters[2],bty="n",cex=1.2)

legend("bottomleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


title(ylab=expression(J[max]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=3.5)
title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3)
abline(v=as.Date("2016-02-26"),lty=3)
abline(v=as.Date("2016-02-03"),lty=3)

# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------



#make negatine PAR values (recorded in night) to zero. 

metwithphy$PAR<-ifelse(metwithphy$PAR<0,0,metwithphy$PAR) 
metwithphy<-subset(metwithphy,!is.na(PAR))  
metwithphy <- metwithphy[order(metwithphy$DateTime_hr), ]
# run room 1 

# Note the Vcmax and Jmax numbers are at growth temperatures.
# no need of temperature response parameters

 # photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha,
 #                                                                          theta =0.2357 , Jmax =Jmax25, Vcmax =Vcmax25
 #                                                                        , Rd0 = 0.68, Q10 = 2.1,TrefR = 22, Rdayfrac=.7,Tcorrect=F))


startdate<-as.Date("2016-01-07")
enddate<-as.Date("2016-02-29")



#set up photosynthesis model 
# get gross GPP, Rday set to zero

photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha=alpha,
                                      theta = 0.28, Jmax =Jmax25, Vcmax =Vcmax25,Tcorrect=F,Rd=0))
  
# photo_gr<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha,
#                                    theta = 0.28, Jmax =Jmax25, Vcmax =Vcmax25,Tcorrect=F,Rd=ifelse(PAR>0,Rd*0.7,Rd)))

   
photom <- photo_gr[,c(1:6, 8)]
photom_15min <- cbind(photom,metwithphy) # add met data 


# convert GPP to g (Carbon) m-1s-1

photom_15min$Date <- as.Date(photom_15min$Date)
photom_15min$gCarbon <- with(photom_15min, ALEAF*15*60*10^-6*12.0107) #g Carbon m-1 s-1
# photom_15min$gCarbon_R <- with(photom_15min, Rd*15*60*10^-6*12.0107) #g Carbon m-1 s-1



#-- get daily total carbon gain for each temperature treatment

dCarbon <- doBy::summaryBy(gCarbon ~ Date+Room, data=photom_15min, FUN=sum, keep.names=TRUE) #g Carbon m-1 day-1  Dark Respiration not included

#dCarbon <- doBy::summaryBy(gCarbon+gCarbon_R~Date+Room, data=photom_15min, FUN=sum, keep.names=TRUE ) #g Carbon m-1 day-1

dCarbon <-subset(dCarbon,Date >startdate & Date <= enddate) 


#---------------------------------------------------------------------------------------------------------------
#--- plot daily carbon gain per m2

windows(40,40);
par(cex.lab=1.5,mar=c(0,7,2,0.5),oma=c(3,2,0,3),cex.axis=1,las=1,mfrow=c(1,1))
plotBy(gCarbon~Date|Room,type="l",col=COL,data=dCarbon,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,12))


title(ylab=expression(C[Day]~(gC~m^-2~d^-1),cex.lab=1.5))
legend("topleft",legend=unique(dCarbon$Room),fill=COL,cex=1.3,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#--- add daily met data to daily carbon gain data

met_daily<-read.csv("Data/MET_DATA/s39_climate_daily_averages.csv")
met_daily$Date<-as.Date(met_daily$Date)
carb_with_met<-merge(dCarbon,met_daily,by=c("Date","Room"))

carb_with_met$PARFac<-cut(carb_with_met$PARsum_mol,breaks=c(0,50,100,150),labels=F)


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#-- add daily total stem/root/leaf mass data (units g(C))

srl_mass<-read.csv("Parameters/great_leaf_area_predicted.csv")
names(srl_mass)[3]<-"LA"
srl_mass$Date<-as.Date(srl_mass$Date,format="%Y-%m-%d")
carb_with_met_la<-merge(carb_with_met,srl_mass,by=c("Date","Room"))

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#--- get daily total GPP (gC plant-1 day-1 )

carb_with_met_la$cday<-with(carb_with_met_la,(gCarbon*LA)) # total GPP per day not corrected for self shading
                                                           # (gross ALEAF), Rday=0
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------


#--- multiply GPP by self shading factor

#--- read self shading factor (used mean sigma for each room)

sigma_s<-read.csv("Parameters/leaf_area_vs_self_shading.csv")

carb_with_met_la<-merge(carb_with_met_la,sigma_s,by="Room")
carb_with_met_la$GPP<-with(carb_with_met_la,cday*(Intercept+Slope*LA))  # total GPP per day with self shading/sigma as a function of LA
# carb_with_met_la$GPP<-with(carb_with_met_la,cday*0.9)  # total GPP per day with self shading/sigma as a function of LA

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

#-------------------------------------------------
#-------------------------------------------------


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

#-- calculate daily maintenance respiration of leaf/stem and root 
#-- respiration numbers are in gC gC-1(biomass)
#-- R (gC gC(biomass))*biomass


carb_with_met_la$Rm_leaf<-with(carb_with_met_la,R_leaf*Leafmass) # leaf mass (gC) 
# carb_with_met_la$Rm_leaf<-with(carb_with_met_la,(gCarbon_R*LA)) # leaf area in m-2
carb_with_met_la$Rm_stem<-with(carb_with_met_la,R_stem*Stemmass)# stem mass (gC) 
carb_with_met_la$Rm_root<-with(carb_with_met_la,R_root*Rootmass)# root mass (gC) 


#-- correction for day respiration

carb_with_met_la$GPP<-with(carb_with_met_la,GPP+(Rm_leaf*0.15)) # assume 12 hr daylight and 30% suppression of Rdark 

write.csv(carb_with_met_la[c(1,2,11,27:36)],file = "Parameters/great_daily_carbon_gain_LA.csv",row.names=FALSE,sep=",")




#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------

#-- plot time vs leaf respiration

windows(100,100);
par(cex.lab=1.5,mar=c(0,7,2,0.5),oma=c(3,2,0,3),cex.axis=1,las=1,mfrow=c(2,2))
plotBy(Rm_leaf~Date|Room,type="l",col=COL,data=carb_with_met_la,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,0.2))
title(ylab=expression(Leaf[R]~(gC~d^-1),cex.lab=1.5,line=4))
legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)

# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

#-- plot time vs stem respiration

plotBy(Rm_stem~Date|Room,type="l",col=COL,data=carb_with_met_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Stem[R]~(gC~d^-1),cex.lab=1.5,line=4))
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

#-- plot time vs stem respiration

plotBy(Rm_root~Date|Room,type="l",col=COL,data=carb_with_met_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Root[R]~(gC~d^-1),cex.lab=1.5,line=4))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------


#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------


#-- get daily total Rm (gC)

# carb_with_met_la$Rm<-with(carb_with_met_la,Rm_stem+Rm_root)

carb_with_met_la$Rm<-with(carb_with_met_la,Rm_leaf+Rm_stem+Rm_root)


#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------


#-- plot time vs total maintenance respiration

plotBy(Rm~Date|Room,type="l",col=COL,data=carb_with_met_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Total[Rm]~(gC~d^-1),cex.lab=1.5,line=4))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)


#-- get total GPP (sum over total duration) and Rm (sum over total duration) for each temperature treatment
#-- units in g(C)
# dat5<-summaryBy(GPP+Rm+Rm_leaf~Room,data=carb_with_met_la,FUN =sum)

dat5<-summaryBy(GPP+Rm~Room,data=carb_with_met_la,FUN =sum)

dat5$Tair.mean<-c(18,21.5,25,28.5,32.5,35.5)


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

#-- get final mass data and add to GPP dataframe

f_mass<-subset(srl_mass,srl_mass$Date==enddate)
# f_mass<-subset(srl_mass,srl_mass$Date==as.Date("2016-02-24"))
f_mass<-transform(f_mass, sum=rowSums(f_mass[c(4:6)])) # get total plant weight in each treatments (g(C))
names(f_mass)[7]<-"totdm_f"


dat5<-merge(dat5,f_mass[c(2:7)],by="Room")

# #-- add final mass given in John's paper
# dat6<-read.csv("Parameters/great_final_mass.csv") # units (g)
# 
# dat5<-merge(dat5,dat6[c(1:2)],by="Room")
# 
# dat5$totdm_f_test<-with(dat5,totdm.mean*.48)


# add initial mass to the dataframe (mass data on 8/1/2016)

dat5$totdm_i<-with(srl_mass[1,],sum(Leafmass,Stemmass,Rootmass)) #units in g(C) and assume similar initial mass for all treatments



#-- get total C use 

dat5$fCarbon_calc<-((dat5$totdm_f-dat5$totdm_i)*1.3)+dat5$Rm.sum  #sum of GPP= sum of Rm (gC)+(1.3*Biomass change)gC

# dat5$fCarbon_calc<-((dat5$totdm_f_test-dat5$totdm_i)*1.3)+dat5$Rm.sum  #test with John's final biomass estimates


dat5$Rg<-((dat5$totdm_f-dat5$totdm_i)*.3)  # calculate growth respiration rate

dat5$storage<-with(dat5,GPP.sum-fCarbon_calc) # calculate balance C (storage??)
dat5$st_perc<-dat5$storage/dat5$GPP


write.csv(dat5,file = "Parameters/great_final_carbon_estimates.csv",row.names=FALSE,sep=",")


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#-- plot carbon balance

windows(100,60);
par(cex.lab=1.5,mar=c(5,5,1,0.5),oma=c(0,0,0,0),cex.axis=1.5,las=1,mfrow=c(1,2))# 

#-- modelled GPP (C)
with(dat5,plot(GPP.sum~Tair.mean,pch=17,cex=2,ylim=c(0,25),ylab="",xlab="",col="red"))
# 
#-- add final mass (C) corrected for Rg (final mass*1.3)
with(dat5,points(fCarbon_calc~Tair.mean,pch=16,cex=2,ylab="",xlab=""))
# adderrorbars(x=dat5$Tair.mean,y=(dat5$totdm.mean*0.47*1.3)+dat5$Rm.sum,SE=dat3$totdm.standard.error*.47,direction="updown")

title(ylab=expression(C~(g)),cex.lab=2,line=1.5)
title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.5)
legend("topright",c("GPP","Used C"),col=c("red","black"),pch=c(17,16),bty="n",ncol=2,cex=2)
legend("bottomleft",legend=(enddate),bty="n")
#  
dat5$Rm_perc<-with(dat5,Rm.sum/GPP.sum*100)
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# plot carbon partitioning

dat_part<-dat5[c("Room","totdm_i","GPP.sum","Rm.sum","Rg","Leafmass","Stemmass","Rootmass","storage")]

dat_part$Rm.sum<-with(dat_part,Rm.sum/GPP.sum)
dat_part$Rg<-with(dat_part, Rg/GPP.sum)
dat_part$Leafmass<-with(dat_part, (Leafmass-srl_mass[1,][[4]])/GPP.sum)
dat_part$Stemmass<-with(dat_part, (Stemmass-srl_mass[1,][[5]])/GPP.sum)
dat_part$Rootmass<-with(dat_part, (Rootmass-srl_mass[1,][[6]])/GPP.sum)
dat_part$storage<-with(dat_part, storage/GPP.sum)


# windows()
# barplot(t(dat_part[c(4:9)]), beside=FALSE, col=rainbow(6), ylim=c(0,2),ylab="",xlab="")
# box()

toplot<-as.matrix(dat_part[c(4:9)])
rownames(toplot)<-c("18","21.5","25","28.5","32.5","35.5")

# windows()
barplot(t(toplot), ylim=c(0, 1.5), ylab = "C Partitioning (%GPP)", xlab = expression(T[growth]),  
        col = rainbow(6),legend = c(expression(R[m]),expression(R[g]),expression(C[Leaf]),expression(C[Stem]),expression(C[Root]),expression(C[Storage])), 
        args.legend = list(x = "top", bty = "n",ncol=3))
# 
# #--- plot storage C
# windows(40,40);
# par(cex.lab=1.5,mar=c(5,5,1,0.5),oma=c(0,0,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))# 
# 
# #-- modelled GPP (C)
# with(dat5,plot(storage~Tair.mean,pch=17,cex=2,ylim=c(-2,10),ylab="",xlab="",col="red"))
# # 
# #-- add final mass (C) corrected for Rg (final mass*1.3)
# with(dat5,points(fCarbon_calc~Tair.mean,pch=16,cex=2,ylab="",xlab=""))
# # adderrorbars(x=dat5$Tair.mean,y=(dat5$totdm.mean*0.47*1.3)+dat5$Rm.sum,SE=dat3$totdm.standard.error*.47,direction="updown")
# 
# title(ylab=expression(C~(g)),cex.lab=2,line=2)
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=2)

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


# #-- plot sigma vs time
# 
# COL <- rev(brewer.pal(6,"Spectral"))
#  windows(40,40);
#  par(cex.lab=1.5,mar=c(0,5,0.5,0.5),oma=c(4,4,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))
# # 
# plotBy(sigma~Date|Room,type="l",col=COL,data=carb_with_met_la,legend=F,lwd=3,ylab="",xlab="",ylim=c(.85,.95))
# # # with(bioch,points(Date,Vcmax25,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
#  axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
# 
#  legend("topright",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
#  
#  
#   title(ylab=expression(sigma),cex.lab=2,line=3.5)
#   title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=2)
#  
# 
# #--------------
# #--------------
# # #-- get the total carbon gain per each treatment
# # 
# carb_total_sum <- doBy::summaryBy(gCarbon_LA ~ Room, data=carb_with_met_la, FUN=sum, keep.names=TRUE )
# mean_T <- doBy::summaryBy(Tair ~ Room, data=met_daily, FUN=mean, keep.names=TRUE )
# carb_total_sum<-merge(carb_total_sum,mean_T,by="Room")
# 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(gCarbon_LA~Tair,data=carb_total_sum,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg="red",xlim=c(15,40),ylim=c(0,20))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# 
# 
# title(xlab=expression(T[growth]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(Total~Carbon~(g)),cex.lab=1.3,line=2)
# 
# 
# #------------------------------------------------------------------------------------------------------------
# #------------------------------------------------------------------------------------------------------------
# 
# # get some plots
# 
# #------------------------------------------------------------------------------------------------------------
# #------------------------------------------------------------------------------------------------------------
# 
#--- plot temperature response of daily carbon gain
# 
# 
# palette(rev(brewer.pal(8,"Set1")))
# COL.1=palette()[c(2,1,3,4)]
# 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(GPP~Tair,data=carb_with_met_la,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.1[PARFac],xlim=c(15,40),ylim=c(0,1.5))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# 
# 
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(Daily~Carbon~(g~day^-1)),cex.lab=1.3,line=2)
# 
# legend("topright",c("< 50","50-100","> 100"),fill=COL.1,cex=1.2,title=expression(Daily~Total~PAR~(mol)),bty="n",horiz=T)

# 
# #--------------------------------------------------------------------------
# 
# # windows();
# # #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# # par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# # plot(gCarbon_LA~Tair,data=carb_with_met_la,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.1[PARFac],xlim=c(15,40),ylim=c(0,2))
# # magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# # #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# # 
# # 
# # title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3,line=2)
# # title(ylab=expression(Daily~Carbon~(g~plant^-1~day^-1)),cex.lab=1.3,line=2)
# # 
# # legend("topright",c("< 50","50-100","> 100"),fill=COL.1,cex=1.2,title=expression(Daily~Total~PAR~(mol)),bty="n",horiz=T)
# 
# 
# #-------------------------------------------------------------------------
#--- plot time vs daily carbon gain
# 
# COL.1=palette()[c(1:6)]
# 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(GPP~Date,data=carb_with_met_la,pch=21,cex=1.5,ylab="",xlab="",col="black",bg=COL.1[Room],ylim=c(0,1.5))
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-01"),to=max(carb_with_met$Date),by="week"),labels=T)
# 
# legend("top",legend=c("18","21.5","25","28.5","32.5","35.5"),fill=COL.1,cex=1.2,title=expression(Daily~T[air]~~(degree*C)),bty="n",horiz=T)
# title(ylab=expression(Daily~Carbon~(g~day^-1)),cex.lab=1.3,line=2)
# # 
# #-------------------------------------------------------------------------
# #-------------------------------------------------------------------------
windows();
plotBy(GPP~Date|Room,type="l",col=COL,data=carb_with_met_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Daily~Carbon~(g~day^-1)),cex.lab=1.3,line=2)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)
legend("top",legend=c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=1.2,title=expression(Daily~T[air]~~(degree*C)),bty="n",ncol=2)

# # Plot T-response of carbon gain in a sunny day
# sunday <- subset(carb_with_met_la, PARsum_mol >= 100)
# unique(sunday$Date)
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(gCarbon_LA~Tair,data=subset(carb_with_met_la,carb_with_met_la$Date=="2016-02-02"),pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg="red",xlim=c(15,40),ylim=c(0,.4))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# 
# 
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(Daily~Carbon~(g~plant^-1~day^-1)),cex.lab=1.3,line=2)
# 
# 
# #-------------------------------------------------------------------------
# #-------------------------------------------------------------------------
# 
# # Plot T-response of carbon gain in a cloudy day
# sunday <- subset(carb_with_met_la, PARsum_mol < 50)
# unique(sunday$Date)
# 
# #windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# #par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# points(gCarbon_LA~Tair,data=subset(carb_with_met_la,carb_with_met_la$Date=="2016-02-04"),pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg="black",xlim=c(15,40),ylim=c(0,.4))
# #magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# 
# 
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(Daily~Carbon~(g~plant^-1~day^-1)),cex.lab=1.3,line=2)

# #-------------------------------------------------------------------------
# #-------------------------------------------------------------------------
# 
# 
# #--- plot temperature response of daily maximum photosynthesis
# 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(ALEAF.max~Tair,data=carb_with_met,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.1[PARFac],xlim=c(15,40),ylim=c(0,40))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# 
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(Daily~Maximum~Photosynthesis[modelled](mu*mol~m^2~s^-2)),cex.lab=1.3,line=2)
# 
# legend("topright",c("< 50","50-100","> 100"),fill=COL.1,cex=1.2,title=expression(Daily~Total~PAR~(mol)),bty="n",horiz=T)
# 
# 
# #---
# 
# #--- 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot( Rd.max~Tair,data=carb_with_met,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.1[2],xlim=c(15,40),ylim=c(0,5))
# points( Rd.mean~Tair,data=carb_with_met,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.1[1],xlim=c(15,40),ylim=c(0,5))
# 
# 
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# #adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
# 
# title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(Daily~Mean~R[leaf](mu*mol~m^2~s^-2)),cex.lab=1.3,line=2)
# 
# 
# #---
# 
# #--- 
# 
# 
# # #--- plot time vs daily carbon gain
# # 
# COL.1=palette()[c(1:6)]
# 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(gCarbon~Date,data=carb_with_met,pch=21,cex=1.5,ylab="",xlab="",col="black",bg=COL.1[Room],ylim=c(0,15))
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-01"),to=max(carb_with_met$Date),by="week"),labels=T)
# 
# legend("top",legend=c("18","21.5","25","28.5","32.5","35.5"),fill=COL.1,cex=1.2,title=expression(Daily~T[air]~~(mu*mol~m^2~s^-1)),bty="n",horiz=T)
# title(ylab=expression(Daily~Carbon~(g~m^-2~day^-1)),cex.lab=1.3,line=2)
# 
# # 
# # 
# # #----
# # 
# # #--Plot Met data
# # 
#  windows();
# # #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
#  par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
#  plotBy(VPD~Date|Room,type="l",col=COL,data=carb_with_met,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,5))
#  plotBy(PARmax~Date|Room,type="l",col=COL,data=carb_with_met,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,2000))
#  
 
 #
# # # 
# # 
# # #test with different dates
# # camp1<-subset(carb_with_met,carb_with_met$Date>"2016-02-01" & carb_with_met$Date<"2016-02-04")
# # with(camp1,plot(Tair,ALEAF.max,col=Room,main="camp_1"))
# # 
# # 
# # camp2<-subset(carb_with_met,carb_with_met$Date>"2016-02-23" & carb_with_met$Date<"2016-02-26")
# # with(camp2,plot(Tair,ALEAF.max,col=Room,main="camp_2"))
