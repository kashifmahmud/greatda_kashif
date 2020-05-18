#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- parameter estimate for daily gpp estimations (GREAT Experiment). 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#' first source functions
source("R/functions_for_analysis.R")

#' stomatal conductance parameters (g1, residual conductance, g0 assumed as 0)
#' used Medlyn et al 2011  stomatal conductance model
#' 
#' data: ACi curves (first observation of each curve ([CO2]~400 ppm))
#-------------------------------------------------------
# read GREAT ACi data

path<-getwd()
gr_aci<-read.csv(paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L0.csv"))

#- get ambient CO2 levels (first observation)

gr_aci.a<-gr_aci[firstobs(~Curve,data=gr_aci),]
gr_aci.a<-subset(gr_aci.a,gr_aci.a$CO2S>300) #remove one datapoint with CO2S~280

#- fit stomatal conductance model for each growth temperatures seperately
#- within each growth temperature, variation of leaf temperatures was ignored
#- so the parameters have effects of both long-term and short-term effects of temperature

tofit.2 <- split(gr_aci.a,gr_aci.a$Room)

g1fits<- get_topts(lapply(tofit.2,FUN=fitStom))
g1fits$Room<-c(1,4,6)
g1fits$Tgrowth<-c(18,28.5,35.5)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#estimate alpha (quantum yield of electron transport)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#alpha=4/PAR*(Ci+2*GammaStar)/(Ci-GammaStar)*(Photo+Rday)

#need some function for temperature response of Rday
#Assume temperature response of Rday is similar to Temperature response of dark respiration (????)

#basal rate for Rday is 0.68 mumolm-2s-1 (GREAT provenance average (short term temperature response of Rdark))
#Q10 for Rday is 2.1 (GREAT provenance average (short term temperature response Rdark))

#--------------------------------------------------

#- get the data (i.e., "long-term" photosynthesis dataset measured at 4 different PAR levels)

 aq <- getAQ()
# 
# #- estimate Gammastar (Source for kinetic parameters: Bernacchi et al 2001; default in Plantecophys)
# 
# aq$GamStar<-TGammaStar(aq$Tleaf,Patm=aq$Press) #estimate Gammastar 
# 
# 
# #- estimate Rday 
# 
# aq$Rday<-TRdayQ10(Tleaf=aq$Tleaf)
# 
# 
# #estimate alpha from low PAR gas exchange measurements
# 
# aq_low<-subset(aq,aq$PARi<150 & aq$W_treatment=="w") # only considered well watered seedlings
# 
# #- calculate alpha
# aq_low$alpha<-with(aq_low,4/PARi*((Ci-2*GamStar)/(Ci-GamStar))*(Photo+Rday))
# 
# #- get mean alpha for each growth treatment
# 
# alpha_mean<-summaryBy(alpha~Room,data=aq_low,FUN=c(mean,std.error))
# alpha_mean$Tgrowth<-c(18,21.5,25,28.5,32,35.5)


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#estimate Vcmax, Jmax and their temperature response parameters

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#fit ACi curves
path<-getwd()
pathtodata<-paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data")


fg<-makecurves(pathtodata,"/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L1.csv")
dg <- makedata.1(pathtodata,"/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L1.csv",fg)

dg<-subset(dg,!dg$Curve %in% c(98,99,103)) #remove three outlier
dg.1<-split(dg,paste(dg$Room)) # ignore provenance differences (Drake et al 2017 GCB)
#dg.1[3]<-NULL  # no data

#- fit temperature response of Vcmax
gr.vcmax<-get_topts(lapply(dg.1,function(x)fitvcmax_mm(x,random="Replicate",return="Peak")))
gr.vcmax$DataSet<-"GRATE"
gr.vcmax$Room<-names(dg.1)


#fit temperature response of Jmax

dg.a<-subset(dg,!is.na(Jmax) & dg$Jmax<400) #function do not work with NAs, so remove few NAs and 2 very high Jmax values
dg.2<-split(dg.a,paste(dg.a$Room))

gr.jmax<-get_topts(lapply(dg.2,function(x)fitjmax_mm(x,random="Replicate",return="Peak")))
gr.jmax$DataSet<-"GRATE"
gr.jmax$Room<-names(dg.2)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#- set up photosynthesis model model
#- first try with ACi data for testing purpose

#- get all parameters to a single dataframe

phy.data<-merge(gr.vcmax[c(1:3,14)],gr.jmax[c(1:3,14)],by="Room") #merge vcmax and jmax to one dataframe

short_g1<-read.csv("Parameters/great_g1_short_term.csv") # short-term tempearture response


phy.data<-merge(short_g1,phy.data,by="Room") #add g1

alph<-read.csv("Parameters/great_alpha_long_term_nls_method.csv") # alpha from nls method

phy.data<-merge(alph,phy.data,by="Room") #add alpha

#- add physiological data to the dataframe with VPD, Tleaf ect... (Ambient datapoints in ACi curves)
aciwithphy<-merge(gr_aci.a,phy.data,by="Room")

#write.csv(phy.data,paste0(path,"/parameters.csv"),row.names=FALSE,sep=",")

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#- model photosynthesis for room 1, 4 and 6 using specified parameters for each room

#room 1
photo_gr_1<-with(subset(aciwithphy,aciwithphy$Room==1),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05
))

photo_gr_1$Room<-1
#----------------------------
#----------------------------
#room 4
photo_gr_4<-with(subset(aciwithphy,aciwithphy$Room==4),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05
))
photo_gr_4$Room<-4
#----------------------------
#----------------------------
#room 6
photo_gr_6<-with(subset(aciwithphy,aciwithphy$Room==6),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05
))


with(photo_gr_1,plot(Tleaf,(Ac-Rd),pch=16,ylab="Ac or AJ or Anet"))
with(photo_gr_1,points(Tleaf,(Aj-Rd),col="red",pch=16))
with(photo_gr_1,points(Tleaf,ALEAF,col="green",pch=16))

legend("bottomleft",c("Ac","Aj","Anet"),col=c("black","red","green"),pch=16,bty="n")

with(photo_gr_4,plot(Tleaf,(Ac-Rd),pch=16,ylab="Ac or AJ or Anet"))
with(photo_gr_4,points(Tleaf,(Aj-Rd),col="red",pch=16))
with(photo_gr_4,points(Tleaf,ALEAF,col="green",pch=16))
legend("bottomleft",c("Ac","Aj","Anet"),col=c("black","red","green"),pch=16,bty="n")


with(photo_gr_6,plot(Tleaf,(Ac-Rd),pch=16,ylab="Ac or AJ or Anet"))
with(photo_gr_6,points(Tleaf,(Aj-Rd),col="red",pch=16))
with(photo_gr_6,points(Tleaf,ALEAF,col="green",pch=16))
legend("bottomleft",c("Ac","Aj","Anet"),col=c("black","red","green"),pch=16,bty="n")




# photo_gr_6$Room<-6
# 
# photo_model<-rbind(photo_gr_1,photo_gr_4,photo_gr_6)
# 
# photo_model_all<-cbind(aciwithphy[c("Tleaf","Photo","Room")],photo_model[c("ALEAF","Room")])
# 
# with(photo_model_all,plot(Photo,ALEAF,xlim=c(0,30),ylim=c(0,30),xlab="Measured A",ylab="Modelled A"))
# fit<-summary(lm(ALEAF~Photo,data=photo_model_all))
# abline(a=coef(fit)[1],b= coef(fit)[2],col="red",lwd=3)
# abline(a=0,b=1)
# legend("topleft",c("1:1 line","Regression line"),col=c("black","red"),lty=1,lwd=3)
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#- plots

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#- plot measured and modelled photosynthesis

windows(100,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,3))


with(subset(aciwithphy,aciwithphy$Room==1),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_1,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 18C",cex=2,bty="n")
title(ylab=expression(A[net]~(mu*mol~m^-2~s^-1)),cex.lab=2,outer=T)
title(xlab=expression(T[leaf]~(degree*C)),cex.lab=2,outer=T)


with(subset(aciwithphy,aciwithphy$Room==4),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_4,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 28.5C",cex=2,bty="n")

with(subset(aciwithphy,aciwithphy$Room==6),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_6,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 35.5C",cex=2,bty="n")

#--------------------------------------------------------------------------------------------------------------

#- plot measured and modelled Ci

windows(100,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,3))


with(subset(aciwithphy,aciwithphy$Room==1),plot(Tleaf,Ci,col="black",pch=16,ylim=c(200,400),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_1,points(Tleaf,Ci,col="red",pch=16,ylim=c(200,400),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 18C",cex=2,bty="n")
title(ylab=expression(C[i]~(ppm)),cex.lab=2,outer=T)
title(xlab=expression(T[growth]~(degree*C)),cex.lab=2,outer=T)


with(subset(aciwithphy,aciwithphy$Room==4),plot(Tleaf,Ci,col="black",pch=16,ylim=c(200,400),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_4,points(Tleaf,Ci,col="red",pch=16,ylim=c(200,400),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 28.5C",cex=2,bty="n")


with(subset(aciwithphy,aciwithphy$Room==6),plot(Tleaf,Ci,col="black",pch=16,ylim=c(200,400),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_6,points(Tleaf,Ci,col="red",pch=16,ylim=c(200,400),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 35.5C",cex=2,bty="n")

#--------------------------------------------------------------------------------------------------------------

#- plot measured and modelled Conductance

windows(100,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,3))


with(subset(aciwithphy,aciwithphy$Room==1),plot(Tleaf,Cond,col="black",pch=16,ylim=c(0,2),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_1,points(Tleaf,GS,col="red",pch=16,ylim=c(0,2),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 18C",cex=2,bty="n")
title(ylab=expression(g[s]),cex.lab=2,outer=T)
title(xlab=expression(T[growth]~(degree*C)),cex.lab=2,outer=T)


with(subset(aciwithphy,aciwithphy$Room==4),plot(Tleaf,Cond,col="black",pch=16,ylim=c(0,2),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_4,points(Tleaf,GS,col="red",pch=16,ylim=c(0,2),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 18C",cex=2,bty="n")



with(subset(aciwithphy,aciwithphy$Room==6),plot(Tleaf,Cond,col="black",pch=16,ylim=c(0,2),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_6,points(Tleaf,GS,col="red",pch=16,ylim=c(0,2),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","Growth temperature = 18C",cex=2,bty="n")



#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#- set up photosynthesis model model with shot-term photosynthesis data

#- read short-term photosynthesis data


aq <- getAQ()
avt <- subset(aq,aq$campaign==1 & aq$W_treatment=="w" )


#- add physiological data to the dataframe with VPD, Tleaf ect... (Ambient datapoints in ACi curves)
phy.data<-read.csv("Parameters/great_alpha_vcmax_jmax_g1.csv")


avtwithphy<-merge(avt,phy.data,by="Room")
avtwithphy$Age_of_plant<-26

#- model photosynthesis for room 1, 4 and 6 using specified parameters for each room

#room 1
photo_gr_avt1<-with(subset(avtwithphy,avtwithphy$Room==1),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,Rdayfrac=.7
))

photo_gr_avt1$Room<-1
#----------------------------
#----------------------------
#room 4
photo_gr_avt4<-with(subset(avtwithphy,avtwithphy$Room==4),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =(Jmax25-(Jmax25*(-0.0302*(40-Age_of_plant)))), Vcmax =(Vcmax25-(Vcmax25*(-0.0302*(40-Age_of_plant)))), Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,Rdayfrac=.7
))
photo_gr_avt4$Room<-4
#----------------------------
#----------------------------
#room 6
photo_gr_avt6<-with(subset(avtwithphy,avtwithphy$Room==6),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =(Jmax25-(Jmax25*(-0.0334*(40-Age_of_plant)))), Vcmax =(Vcmax25-(Vcmax25*(-0.0334*(40-Age_of_plant)))), Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,Rdayfrac=.7
))
photo_gr_avt6$Room<-6
# 
(112.70-(112.70*(0.009387171*(40-26))))
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#- plot measured and modelled photosynthesis

windows(100,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,3))


with(subset(avtwithphy,avtwithphy$Room==1),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_avt1,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","18C",cex=1,bty="n")
title(ylab=expression(A[net]~(mu*mol~m^-2~s^-1)),cex.lab=2,outer=T)
title(xlab=expression(T[leaf]~(degree*C)),cex.lab=2,outer=T)


with(subset(avtwithphy,avtwithphy$Room==4),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_avt4,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","28.5C",cex=1,bty="n")

with(subset(avtwithphy,avtwithphy$Room==6),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_avt6,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","35.5C",cex=1,bty="n")

#--------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#- set up photosynthesis model model with long-term photosynthesis data (campaign 2)

#- read short-term photosynthesis data


avt_last <- subset(aq,aq$campaign==2 & aq$W_treatment=="w")

#- add physiological data to the dataframe with VPD, Tleaf ect... (Ambient datapoints in ACi curves)
avtwithphy_last<-merge(avt_last,phy.data,by="Room")
avtwithphy_last$Age_of_plant<-49

#- model photosynthesis for room 1, 4 and 6 using specified parameters for each room

#room 1
photo_gr_avt1<-with(subset(avtwithphy_last,avtwithphy_last$Room==1),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                   gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                   theta = 0.85, Jmax =(Jmax25-(Jmax25*(0.009387171*(40-Age_of_plant)))), Vcmax =(Vcmax25-(Vcmax25*(0.009387171*(40-Age_of_plant)))), Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                   delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,Rdayfrac=.7
))

photo_gr_avt1$Room<-1
#----------------------------
#----------------------------
#room 4
photo_gr_avt4<-with(subset(avtwithphy_last,avtwithphy_last$Room==4),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                   gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                   theta = 0.85, Jmax =(Jmax25-(Jmax25*(-0.0302*(40-Age_of_plant)))), Vcmax =(Vcmax25-(Vcmax25*(-0.0302*(40-Age_of_plant)))), Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                   delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,Rdayfrac=.7
))
photo_gr_avt4$Room<-4
#----------------------------
#----------------------------
#room 6
photo_gr_avt6<-with(subset(avtwithphy_last,avtwithphy_last$Room==6),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
                                                                   gsmodel = "BBOpti", g1=g1,alpha = alpha,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                   theta = 0.85, Jmax =(Jmax25-(Jmax25*(-0.0334*(40-Age_of_plant)))), Vcmax =(Vcmax25-(Vcmax25*(-0.0334*(40-Age_of_plant)))), Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                   delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,Rdayfrac=.7
))
photo_gr_avt6$Room<-6


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#- plot measured and modelled photosynthesis

windows(100,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,3))


with(subset(avtwithphy_last,avtwithphy_last$Room==1),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_avt1,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","18C",cex=1,bty="n")
title(ylab=expression(A[net]~(mu*mol~m^-2~s^-1)),cex.lab=2,outer=T)
title(xlab=expression(T[leaf]~(degree*C)),cex.lab=2,outer=T)


with(subset(avtwithphy_last,avtwithphy_last$Room==4),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_avt4,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","28.5C",cex=1,bty="n")

with(subset(avtwithphy_last,avtwithphy_last$Room==6),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr_avt6,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,las=1,cex.axis=1.1)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=2,bty="n")
legend("bottomright","35.5C",cex=1,bty="n")

#--------------------------------------------------------------------------------------------------------------
