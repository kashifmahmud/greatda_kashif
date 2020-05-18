#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- Testing daily GPP estimates of Leaf scale model with measured GPP data (GREAT Experiment). 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


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

# 
# #--- read g1
# 
# # short_g1<-read.csv("Parameters/great_g1_long_term.csv") # short-term tempearture response
# 
# short_g1<-read.csv("Parameters/great_g1_short_term.csv") # short-term tempearture response
# 
# 
# 
# #--- read alpha
# 
# alph<-read.csv("Parameters/great_alpha_long_term_nls_method.csv") # alpha from nls method
# 
# 
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# 
# # make dataframe with g1 and alpha
# 
# phy.data<-merge(short_g1[c(1:3)],alph,by="Room") 
# 
# 
# 
# #add fixed biochemical parameters (Ea and deltaS of Vcmax and Jmax)
# 
# bioch_fixed<-read.csv("Parameters/great_aci_t_response_fixed.csv")
# # bioch_fixed$Date<-as.Date(bioch_fixed$Date)
# 
# phy.data<- merge(phy.data,bioch_fixed,by="Room",all=T) 
# 
# 
# 
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# #--- read vcmax, jmax and temperature response parameters (fitted to GREAT ACi data)
# 
# bioch<-read.csv("Parameters/great_aci_t_response.csv")[c(1,2,9,17)]
# bioch$Date<-as.Date(bioch$Date)
# 
# phy.data<- merge(phy.data,subset(bioch,bioch$Date==as.Date("2016-02-03")),by="Room",all=T) 
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#--Get data and add physiological parameters

aq <- getAQ()
phy.data<-read.csv("Parameters/great_alpha_vcmax_jmax_25.csv")
aq_phys<-merge(aq,phy.data,by="Room")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# read Rday (corrected for rday frac =0.7)

rday<-read.csv("Parameters/great_Resp_leaf_shoot_root_25C.csv")

aq_phys<-merge(aq_phys,rday,by="Room")

aq.list <- subset(aq_phys, aq_phys$campaign==1 & aq_phys$W_treatment=="w") # campaign 1 well watered trees low and high light

# aq_phys<-merge(aq.list,phy.data,by="Room")


#-- set up photosynthesis model 

photo_gr<-with(aq.list,Photosyn(VPD=VpdL, Ca=CO2S,Tleaf=Tleaf,PPFD=PARi,gsmodel = "BBOpti", g1=g1, alpha = alpha,
                                theta = theta, Jmax =Jmax25, Vcmax =Vcmax25,delsC=delsV*10^3,delsJ=delsJ*10^3,EaV=EaV*10^3,EaJ=EaJ*10^3,
                                   Rd0=R_leaf_a,Rdayfrac=.7,Q10=2.1,TrefR=25))
photo_model_all<-cbind(aq.list,photo_gr)

photo_model_all$resid<-with(photo_model_all,(ALEAF-Photo)/Photo)


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
palette(rev(brewer.pal(6,"Spectral")))
COL=palette()[c(1:6)]

library(scales)
# windows(100,40);


png(file="output/model_testing.png",width=300,height=500)
par(cex.lab=1.5,mar=c(5,5,0,0.5),oma=c(2,0,2,0),cex.axis=1.5,las=1,mfrow=c(3,1))


with(photo_model_all,plot(Photo,ALEAF,col=alpha("black",0.2),pch=16,ylim=c(0,40),xlim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
abline(0,1)
title(ylab=expression(A[Leaf~model]),cex.lab=2)
title(xlab=expression(A[Measured]),cex.lab=2)
abline(summary(lm(ALEAF~Photo,data=photo_model_all))$coefficients[1],summary(lm(ALEAF~Photo,data=photo_model_all))$coefficients[2],lty=3,col="red",lwd=3)


with(aq.list,plot(Tleaf,Photo,col=alpha("black",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
with(photo_gr,points(Tleaf,ALEAF,col=alpha("red",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
title(xlab=expression(T[leaf]),cex.lab=2)
legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=1,bty="n",ncol=2)
title(ylab=expression(A[Leaf~model]),cex.lab=2)
# lm1<-summary(lm(ALEAF~Photo,data=photo_model_all))
# plot(resid(lm1),col=Room)


with(photo_model_all,plot(Tleaf,resid,col=COL[LightFac],pch=16,ylim=c(-1,3),xlim=c(15,45),cex=1,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
title(xlab=expression(T[Growth]),cex.lab=2)
title(ylab=expression(Residuals),cex.lab=2)
legend("topleft",legend=(c("100","500","1000","1500")),col=COL[unique(photo_gr_summary$LightFac)],pch=19,title="PAR",ncol=2)

dev.off()
#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#plot mean photosynthesis and modelled photosynthesis

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#plot mean photosynthesis and modelled photosynthesis
# windows(100,40);
# par(cex.lab=1.5,mar=c(.5,.5,.5,5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,2))

photo_gr_summary<-summaryBy(ALEAF+Photo+resid+Tleaf~Room+LightFac,data=photo_model_all,FUN=c(mean,std.error))
# photo_gr_summary$resid<-with(photo_gr_summary,(ALEAF.mean-Photo.mean)/Photo.mean)

with(photo_gr_summary,plot(Photo.mean,ALEAF.mean,col=COL[Room],pch=16,ylim=c(0,30),xlim=c(0,30),cex=2.5,ylab="",xlab="",axes=F))
adderrorbars(x=photo_gr_summary$Photo.mean,y=photo_gr_summary$ALEAF.mean,SE=photo_gr_summary$ALEAF.std.error*1.96,direction="updown")
adderrorbars(x=photo_gr_summary$Photo.mean,y=photo_gr_summary$ALEAF.mean,SE=photo_gr_summary$Photo.std.error*1.96,direction="leftright")
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
abline(0,1)
abline(summary(lm(ALEAF~Photo,data=photo_model_all))$coefficients[1],summary(lm(ALEAF~Photo,data=photo_model_all))$coefficients[2],lty=3,col="red")
legend("topleft",legend=unique(photo_gr_summary$Room),col=COL[unique(photo_gr_summary$Room)],pch=19,title="Room",ncol=2,cex=1.3)
title(ylab=expression(A[Leaf~model]),cex.lab=2)
title(xlab=expression(A[Measured]),cex.lab=2,line=2)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
with(photo_gr_summary,plot(resid.mean~Tleaf.mean,col=COL[LightFac],pch=16,ylim=c(-.5,.5),xlim=c(15,40),cex=1.5,ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
legend("topleft",legend=(c("100","500","1000","1500")),col=COL[unique(photo_gr_summary$LightFac)],pch=19,title="PAR",ncol=2)
abline(h=0,lty=3)
title(xlab=expression(T[Growth]),cex.lab=2)
title(ylab=expression(Residuals),cex.lab=2)

dev.off()

# #-------------------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------------------

# # test for campaign 2 AQ data
# 
# aq.list_2 <- subset(aq_phys, aq_phys$campaign==2 & aq_phys$W_treatment=="w") # campaign 2 well watered trees low and high light
# 
# # aq_phys<-merge(aq.list,phy.data,by="Room")
# 
# 
# #-- set up photosynthesis model 
# 
# photo_gr_c2<-with(aq.list_2,Photosyn(VPD=VpdL, Ca=CO2S,Tleaf=Tleaf,PPFD=PARi,gsmodel = "BBOpti", g1=g1, alpha = alpha,
#                                 theta = theta, Jmax =Jmax25, Vcmax =Vcmax25
#                                 , Rd0 = R_leaf_a, Q10 = 2.1,TrefR = 25,
#                                 delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,EaV = EaV*10^3,
#                                 EaJ = EaJ*10^3,Rdayfrac=.7))
# photo_model_all_c2<-cbind(aq.list_2,photo_gr_c2)
# 
# 
# #-------------------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------------------
# 
# windows(100,40);
# par(cex.lab=1.5,mar=c(.5,.5,.5,5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,2))
# 
# 
# with(photo_model_all_c2,plot(Photo,ALEAF,col=alpha("black",0.2),pch=16,ylim=c(0,40),xlim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# abline(0,1)
# title(ylab=expression(A[Leaf~model]),cex.lab=2,outer=T)
# title(xlab=expression(A[Measured]),cex.lab=2,outer=T,adj=.2)
# 
# 
# 
# with(aq.list_2,plot(Tleaf,Photo,col=alpha("black",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# with(photo_gr_c2,points(Tleaf,ALEAF,col=alpha("red",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# title(xlab=expression(T[leaf]),cex.lab=2,outer=T,adj=.72)
# legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=1,bty="n",ncol=2)
# 
# 
# 
# 
# 
# 
# 
# 
# #-------------------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------------------
# 
# #-- set up photosynthesis model with adjusted EaJ and delsJ
# 
# photo_gr<-with(aq_phys,Photosyn(VPD=VpdL, Ca=CO2S,Tleaf=Tleaf,PPFD=PARi,gsmodel = "BBOpti", g1=g1, alpha = alpha,
#                                 theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25
#                                 , Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                 delsC = delsV*10^3,EdVC = 2e+05,EdVJ = 2e+05,EaV = EaV*10^3,
#                                 Rdayfrac=.7,EaJ=42.74635*10^3,delsJ=0.6312849*10^3))
# 
# 
# photo_model_all<-cbind(aq_phys,photo_gr)
# 
# 
# windows(100,40);
# par(cex.lab=1.5,mar=c(.5,.5,.5,5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,2))
# 
# with(photo_model_all,plot(Photo,ALEAF,col=alpha("black",0.2),pch=16,ylim=c(0,40),xlim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# abline(0,1)
# title(ylab=expression(A[Leaf~model]),cex.lab=2,outer=T)
# title(xlab=expression(A[Measured]),cex.lab=2,outer=T,adj=.2)
# 
# 
# 
# with(aq.list,plot(Tleaf,Photo,col=alpha("black",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# with(photo_gr,points(Tleaf,ALEAF,col=alpha("red",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# title(xlab=expression(T[leaf]),cex.lab=2,outer=T,adj=.72)
# legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=1,bty="n",ncol=2)
# summary(lm(ALEAF~Photo,data=photo_model_all))
# 
# #-------------------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------------------
# 
# #plot mean photosynthesis and modelled photosynthesis
# windows(100,40);
# par(cex.lab=1.5,mar=c(.5,.5,.5,5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,2))
# 
# photo_gr_summary<-summaryBy(ALEAF+Photo+Tleaf~Room+LightFac,data=photo_model_all,FUN=c(mean,std.error))
# photo_gr_summary$resid<-with(photo_gr_summary,(ALEAF.mean-Photo.mean)/Photo.mean)
# with(photo_gr_summary,plot(Photo.mean,ALEAF.mean,col=COL[Room],pch=16,ylim=c(0,40),xlim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# adderrorbars(x=photo_gr_summary$Photo.mean,y=photo_gr_summary$ALEAF.mean,SE=photo_gr_summary$ALEAF.std.error*1.96,direction="updown")
# adderrorbars(x=photo_gr_summary$Photo.mean,y=photo_gr_summary$ALEAF.mean,SE=photo_gr_summary$Photo.std.error*1.96,direction="leftright")
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# abline(0,1)
# abline(summary(lm(ALEAF~Photo,data=photo_model_all))$coefficients[1],summary(lm(ALEAF~Photo,data=photo_model_all))$coefficients[2],lty=3,col="red")
# legend("topleft",legend=unique(photo_gr_summary$Room),col=COL[unique(photo_gr_summary$Room)],pch=19,title="Room")
# title(ylab=expression(A[Leaf~model]),cex.lab=2,outer=T)
# title(xlab=expression(A[Measured]),cex.lab=2,outer=T,adj=.2)
# 
# #-------------------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------------------
# with(photo_gr_summary,plot(resid~Tleaf.mean,col=COL[LightFac],pch=16,ylim=c(-.5,.5),xlim=c(15,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# legend("topleft",legend=unique(photo_gr_summary$LightFac),col=COL[unique(photo_gr_summary$LightFac)],pch=19,,title="PAR")
# abline(h=0,lty=3)
# title(ylab=expression(Residuals),cex.lab=1.8,outer=T,line=-22)
# title(xlab=expression(T[Growth]),cex.lab=2,outer=T,adj=.72)
# 
# 
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------
# 
# #--test the model with short term photosynthesis measurements
# 
# avt <- getAvT()
# 
# avt_phys<-merge(avt,phy.data,by="Room")
# 
# avt.list <- subset(avt_phys, avt_phys$W_treatment=="w") # campaign 1 well watered trees low and high light
# 
# 
# #-- set up photosynthesis model 
# 
# photo_gr_c3<-with(avt.list,Photosyn(VPD=VpdL, Ca=CO2S,Tleaf=Tleaf,PPFD=PARi,gsmodel = "BBOpti", g1=g1, alpha = alpha,
#                                      theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25
#                                      , Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                      delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05,EaV = EaV*10^3,
#                                      EaJ = EaJ*10^3,Rdayfrac=.7))
# photo_model_all_c3<-cbind(avt.list,photo_gr_c3)
# 
# 
# #-------------------------------------------------------------------------------------------------------------------
# #-------------------------------------------------------------------------------------------------------------------
# 
# windows(100,40);
# par(cex.lab=1.5,mar=c(.5,.5,.5,5),oma=c(7,7,1,7),cex.axis=1.5,las=1,mfrow=c(1,2))
# 
# 
# with(photo_model_all_c3,plot(Photo,ALEAF,col=alpha("black",0.2),pch=16,ylim=c(0,40),xlim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# abline(0,1)
# title(ylab=expression(A[Leaf~model]),cex.lab=2,outer=T)
# title(xlab=expression(A[Measured]),cex.lab=2,outer=T,adj=.2)
# 
# 
# 
# with(avt.list,plot(Tleaf,Photo,col=alpha("black",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# with(photo_model_all_c3,points(Tleaf,ALEAF,col=alpha("red",0.5),pch=16,ylim=c(0,40),cex=1.5,ylab="",xlab="",axes=F))
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1)
# title(xlab=expression(T[leaf]),cex.lab=2,outer=T,adj=.72)
# legend("top",c("Measured","Modelled"),pch=16,col=c("black","red"),cex=1,bty="n",ncol=2)
# 
# 
