#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#- This plots temperature response curves forlight-saturated photosynthetic rates
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

library(nlshelper)


#-----------------------------------------------------------------------------------------




#--- Asat
#- get the AQ data (i.e., "long-term" photosynthesis dataset)
aq <- getAQ()

Asatdata <- subset(aq,campaign==1 & W_treatment=="w" & !is.na(Photo))

#- plot temperature response of Asat (long-term)
#- get TREATMENT MEANS for photosynthesis (average across replicate seedlings)
#- one point in figure represents one provenance, but provenances not showed differently as they are similar (Drake et al 2017, GCB)

dat3 <- summaryBy(PARi+Photo+Tair+Tleaf~Room+location+LightFac,FUN=c(mean,std.error),data=Asatdata,na.rm=T)


#get a colour palette

#palette(rev(brewer.pal(11,"BrBG")))

palette(rev(brewer.pal(9,"YlOrRd")))
COL=palette()[c(1,2,3,8)]


windows();
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
plot(Photo.mean~Tleaf.mean,data=dat3,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL[LightFac],ylim=c(0,30))
adderrorbars(x=dat3$Tleaf.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")
#adderrorbars(x=dat3$Tleaf.mean,y=dat3$Photo.mean,SE=dat3$Tleaf.standard.error,direction="leftright")


magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.3,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)

legend(28,30,c("100","500","1000","1500"),fill=COL,cex=0.9,title=expression(PAR~(mu*mol~m^-2~s^-1)),bty="n",horiz=T)


#- estimate temperature optimum for photosynthesis
#- AQ data-long-term (campaign 1)

#tofit_long<-split(dat3,paste(dat3$LightFac))

tofit_long<-split(Asatdata,paste(Asatdata$LightFac))
topt_long<-get_topts(lapply(tofit_long,function(x)fit.nlme(x,yvar="Photo",xvar="Tleaf",random="Prov")))
topt_long$LightFac<-c(1,2,3,4)


#- add Topt for photosynthesis to same plot 

points(aopt~topt,data=topt_long,pch=21,cex=2,ylab="",xlab="",col="black",bg=COL[LightFac],ylim=c(0,30),lwd=1)

adderrorbars(x=topt_long$topt,y=topt_long$aopt,SE=1.96*(topt_long$topt.se),direction="leftright")


# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(topt~LightFac,data=topt_long,pch=21,cex=2,ylab="",xlab="",axes=F,col="black",bg="lightgrey",ylim=c(20,30))
# adderrorbars(x=topt_long$LightFac,y=topt_long$topt,SE=1.96*(topt_long$topt.se),direction="updown")
# 
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# title(xlab=expression(Measurement~T[opt]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(PAR[i]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)




#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
palette(rev(brewer.pal(8,"Set1")))
COL.1=palette()[c(1:6)]


windows();
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
plot(Photo.mean~PARi.mean,data=dat3,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.1[Room],xlim=c(0,1500),ylim=c(0,30))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
#adderrorbars(x=dat3$PARi.mean,y=dat3$Photo.mean,SE=dat3$Photo.standard.error,direction="updown")


title(xlab=expression(Measurement~PAR~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)

legend(0,30,c("18","21.5","25","28.5","32","35.5"),fill=COL.1,cex=0.9,title=expression(Measurement~T[leaf]~(degree*C)),bty="n",horiz=T)



#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#function to fit non-rectangular hypobola and get fitted parameters
# --- Nonlinear least squares regression (non-rectangular hyperbola).  
# --- 4 parameter model: Amax (max gross photosytnthetic rate), Rd (dark respiration), AQY (apparent quantum yield), Theta (curvature parameter, dimensionless) ---

# fit_nrh<-function(data){
#   
#   curve.nlslrc = nls(Photo ~ (1/(2*theta))*(AQY*PARi+Amax-sqrt((AQY*PARi+Amax)^2-4*AQY*theta*Amax*PARi))-Rd,
#                      start=list(Amax=(max(data$Photo)),AQY=0.05,Rd=-min(data$Photo),theta=1),data=data) 
#   
#   summary(curve.nlslrc) #summary of model fit
#   
# }
# 
# curve.nlslrc = nls(Photo ~ (1/(2*theta))*(AQY*PARi+Amax-sqrt((AQY*PARi+Amax)^2-4*AQY*theta*Amax*PARi))-Rd,
#                    start=list(Amax=25,AQY=0.05,Rd=-.5,theta=0.8),data=Asatdata) 
# plot_nls(curve.nlslrc)

#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#- plot temperature response of photosynthesis in campaign 2 data
#- only provenance B (Warm-edge) measured

Asatdata.c2 <- subset(aq,campaign==2 & W_treatment=="w")

dat4 <- summaryBy(PARi+Photo+Tair+Tleaf~Room+LightFac,FUN=c(mean,standard.error),data=Asatdata.c2,na.rm=T)

windows();
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
plot(Photo.mean~Tleaf.mean,data=dat4,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL[LightFac],ylim=c(0,30))
adderrorbars(x=dat4$Tleaf.mean,y=dat4$Photo.mean,SE=dat4$Photo.standard.error,direction="updown")
#adderrorbars(x=dat4$Tleaf.mean,y=dat4$Photo.mean,SE=dat4$Tleaf.standard.error,direction="leftright")


magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.3,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)

legend(28,30,c("100","500","1000","1500"),fill=COL,cex=0.9,title=expression(PAR~(mu*mol~m^-2~s^-1)),bty="n",horiz=T)


#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#- plot temperature response of photosynthesis in ACi data

# read GREAT ACi data

path<-getwd()
gr_aci<-read.csv(paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L0.csv"))

#- get ambient CO2 levels (first observation)

gr_aci.a<-gr_aci[firstobs(~Curve,data=gr_aci),]
gr_aci.a<-subset(gr_aci.a,gr_aci.a$CO2S>300) #remove one datapoint with CO2S~280

#- assign the temperature levels
gr_aci.a$TleafFac <- cut(gr_aci.a$Tleaf,breaks=c(15,22,27,32,38,45),labels=1:5)



dat5 <- summaryBy(Photo+Tair+Tleaf~Room+Provenance+TleafFac,FUN=c(mean,std.error),data=gr_aci.a,na.rm=T)
dat5<-subset(dat5,dat5$Provenance %in% c("A","B"))

palette(rev(brewer.pal(8,"Set1")))
COL.2=palette()[c(1:6)]

windows();
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
plot(Photo.mean~Tleaf.mean,data=dat5,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.2[Room],ylim=c(0,30))
adderrorbars(x=dat5$Tleaf.mean,y=dat5$Photo.mean,SE=dat5$Photo.standard.error,direction="updown")
#adderrorbars(x=dat5$Tleaf.mean,y=dat5$Photo.mean,SE=dat5$Tleaf.standard.error,direction="leftright")


magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.3,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)

legend("bottomleft",c("18","28.5","35.5"),fill=COL.2[c(1,4,6)],cex=0.9,title=expression(Growth~Temperature~(degree*C)),bty="n",horiz=T)

arrows(x0=20,y0=26,x1=20,y1=24,lwd=3,col=COL.2[1])
arrows(x0=28,y0=19,x1=30,y1=17,lwd=3,col=COL.2[4])
arrows(x0=35,y0=20,x1=35,y1=18,lwd=3,col=COL.2[6])

#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#- fit.nlme(dat5,yvar="Photo.mean",xvar="Tleaf.mean",random="Provenance")
#- estimate temperature optimum for photosynthesis
#- ACI data-Short term

tofit_aci<-split(gr_aci.a,gr_aci.a$Room)

aci.topts<-get_topts(lapply(tofit_aci[c(1,3)],function(x)fit.nlme(x,yvar="Photo",xvar="Tleaf",random="Provenance")))
aci.topts$Room<-names(tofit_aci)[c(1,3)]

aci.topts.1<-get_topts(lapply(tofit_aci[2],function(x)fitquad(x)))
aci.topts.1$Room<-names(tofit_aci[2])

aci.topts<-rbind(aci.topts,aci.topts.1)


#- add Topt for photosynthesis to same plot 

points(aopt~topt,data=aci.topts,pch=21,cex=2,ylab="",xlab="",col="black",bg=COL.2[as.numeric(Room)],ylim=c(0,30),lwd=3)

adderrorbars(x=aci.topts$topt,y=aci.topts$aopt,SE=1.96*(aci.topts$topt.se),direction="leftright")



#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

# #- estimate temperature optimum for photosynthesis
# #- AQ data-long-term (campaign 1)
# 
# #tofit_long<-split(dat3,paste(dat3$LightFac))
# 
# tofit_long<-split(Asatdata,paste(Asatdata$LightFac))
# topt_long<-get_topts(lapply(tofit_long,function(x)fit.nlme(x,yvar="Photo",xvar="Tleaf",random="Prov")))
# topt_long$LightFac<-c(100,500,1000,1500)
# 
# 
# windows();
# #pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
# par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
# plot(topt~LightFac,data=topt_long,pch=21,cex=2,ylab="",xlab="",axes=F,col="black",bg="lightgrey",ylim=c(20,30))
# adderrorbars(x=topt_long$LightFac,y=topt_long$topt,SE=1.96*(topt_long$topt.se),direction="updown")
# 
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
# title(xlab=expression(Measurement~T[opt]~(degree*C)),cex.lab=1.3,line=2)
# title(ylab=expression(PAR[i]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)


#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#- estimate temperature optimum for photosynthesis
#- AQ data-long-term (campaign 2): Warm-Edge prov only


# tofit_long_c2<-split(Asatdata.c2,paste(Asatdata.c2$LightFac))
# topt_long<-get_topts(lapply(tofit_long_c2,function(x)fitquad(x)))

#cannot fit models. no peaks in photosynthesis.

#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#- plot short-term photosynthesis data

a_short<-getAvT()

#- assign light levels to a factor variable
a_short$LightFac <- NA
a_short$LightFac[which(a_short$PARi<150)] <- 1
a_short$LightFac[which(a_short$PARi>400) ] <- 4 
a_short$LightFac <- as.factor(a_short$LightFac)

#- plot temperature response of Asat (short-term)
#- get TREATMENT MEANS for photosynthesis (average across replicate seedlings)
#- one point in figure represents one provenance, but provenances not showed differently as they are similar (Drake et al 2017, GCB)

a_short_means <- summaryBy(PARi+Photo+Tair+Tleaf~Room+Prov+LightFac,FUN=c(mean,standard.error),data=a_short,na.rm=T)

palette(rev(brewer.pal(11,"BrBG")))
COL.3=palette()[c(11,8)]


windows();
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
plot(Photo.mean~Tleaf.mean,data=a_short_means,pch=21,cex=1.5,ylab="",xlab="",axes=F,col="black",bg=COL.3[LightFac],ylim=c(0,35))
adderrorbars(x=a_short_means$Tleaf.mean,y=a_short_means$Photo.mean,SE=a_short_means$Photo.standard.error,direction="updown")
#adderrorbars(x=a_short_means$Tleaf.mean,y=a_short_means$Photo.mean,SE=a_short_means$Tleaf.standard.error,direction="leftright")


magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.3,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)

legend("topleft",c("100","1500"),fill=COL.3,cex=0.9,title=expression(PAR~(mu*mol~m^-2~s^-1)),bty="n",horiz=T)


#- estimate temperature optimum for photosynthesis
#- AQ data-long-term (campaign 1)

#tofit_long<-split(a_short_means,paste(a_short_means$LightFac))

tofit_short<-split(a_short,paste(a_short$LightFac))
topt_short<-get_topts(lapply(tofit_short,function(x)fit.nlme(x,yvar="Photo",xvar="Tleaf",random="Prov")))
topt_short$LightFac<-factor(c(1,4))


#- add Topt for photosynthesis to same plot 

points(aopt~topt,data=topt_short,pch=21,cex=2,ylab="",xlab="",col="black",bg=COL.3[LightFac],ylim=c(0,30),lwd=1)

adderrorbars(x=topt_short$topt,y=topt_short$aopt,SE=1.96*(topt_short$topt.se),direction="leftright")

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

# get temperature response of warm-edge provenance for comparison

#- long-term data
prov_B_l <- subset(dat3,location=="Warm-edge")


palette(rev(brewer.pal(9,"YlOrRd")))
COL=palette()[c(1,2,3,8)]

windows(100,40);
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0),mfrow=c(1,3))
plot(Photo.mean~Tleaf.mean,data=prov_B_l,pch=21,cex=2,ylab="",xlab="",axes=F,col="black",bg=COL[LightFac],ylim=c(0,30))
adderrorbars(x=prov_B_l$Tleaf.mean,y=prov_B_l$Photo.mean,SE=prov_B_l$Photo.standard.error,direction="updown")
#adderrorbars(x=prov_B_l$Tleaf.mean,y=prov_B_l$Photo.mean,SE=prov_B_l$Tleaf.standard.error,direction="leftright")


magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.5,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.5,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=2)

legend(22,30,c("100","500","1000","1500"),fill=COL,cex=1.3,title=expression(PAR~(mu*mol~m^-2~s^-1)),bty="n",horiz=T)
legend("bottomright","in-situ (03/02/2016)",cex=1.3,bty="n")

arrows(x0=20.5,y0=21,x1=20.5,y1=19,lwd=3,col=COL.2[1])
arrows(x0=35,y0=25,x1=33,y1=25,lwd=3,col=COL.2[4])
arrows(x0=39.5,y0=19,x1=39.5,y1=17,lwd=3,col=COL.2[6])


#- plot campaign 2 in the same plot

plot(Photo.mean~Tleaf.mean,data=dat4,pch=21,cex=2,ylab="",xlab="",axes=F,col="black",bg=COL[LightFac],ylim=c(0,30))
adderrorbars(x=dat4$Tleaf.mean,y=dat4$Photo.mean,SE=dat4$Photo.standard.error,direction="updown")
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.5,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.5,line=2)
#title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)
legend(22,30,c("100","500","1000","1500"),fill=COL,cex=1.3,title=expression(PAR~(mu*mol~m^-2~s^-1)),bty="n",horiz=T)
legend("bottomright","in-situ (26/02/2016)",cex=1.3,bty="n")


arrows(x0=20.5,y0=27.5,x1=20.5,y1=25.5,lwd=3,col=COL.2[1])
arrows(x0=32.8,y0=17,x1=32.8,y1=15,lwd=3,col=COL.2[4])
arrows(x0=40.3,y0=14,x1=40.3,y1=12,lwd=3,col=COL.2[6])


#- plot warm edge prov from ACi data

asat_B_aci<-subset(dat5,dat5$Provenance %in% c("B"))


palette(rev(brewer.pal(8,"Set1")))
COL.2=palette()[c(1:6)]

plot(Photo.mean~Tleaf.mean,data=asat_B_aci,pch=21,cex=2,ylab="",xlab="",axes=F,col="black",bg=COL.2[Room],ylim=c(0,30))
adderrorbars(x=asat_B_aci$Tleaf.mean,y=asat_B_aci$Photo.mean,SE=asat_B_aci$Photo.standard.error,direction="updown")
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.5,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.5,line=2)
#title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)

legend("topright",c("18","28.5","35.5"),fill=COL.2[c(1,4,6)],cex=1.3,title=expression(Growth~Temperature~(degree*C)),bty="n",horiz=T)

arrows(x0=20,y0=26,x1=20,y1=24,lwd=3,col=COL.2[1])
arrows(x0=31,y0=10,x1=31,y1=12,lwd=3,col=COL.2[4])
arrows(x0=35.5,y0=11,x1=35.5,y1=13,lwd=3,col=COL.2[6])

legend("bottomright","ACi data (16-25/02/2016)",cex=1.3,bty="n")



