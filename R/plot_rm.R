#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------

palette(rev(brewer.pal(6,"Spectral")))
COL=palette()[c(1:6)]
#-- plot time vs leaf respiration

# windows(100,100);
png(file="output/leaf_respiration.png",width=800,height=350)
par(cex.lab=1.5,mar=c(2,5,0,0.5),oma=c(2,0,2,0),cex.axis=1,las=1,mfrow=c(1,3))


plotBy(R_leaf~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,0.15))
title(ylab=expression(Leaf[R]~(gC~d^-1),cex.lab=1,line=4))
legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=1.2,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)
legend("bottomright",legend="Total leaf respiration")


plotBy(R_leaf_n~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,0.15))
title(ylab=expression(Leaf[R]~(gC~d^-1),cex.lab=1,line=4))
# legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=1.2,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)
legend("bottomright",legend="Nitht time")



plotBy(R_leaf_d~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,0.15))
title(ylab=expression(Leaf[R]~(gC~d^-1),cex.lab=1,line=4))
# legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=1.2,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
legend("bottomright",legend="Day time")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

dev.off()
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)



png(file="output/stem_root_respiration.png",width=500,height=350)
par(cex.lab=1.5,mar=c(2,5,0,0.5),oma=c(2,0,2,0),cex.axis=1,las=1,mfrow=c(1,2))

#-- plot time vs stem respiration

plotBy(R_stem~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Stem[R]~(gC~d^-1),cex.lab=1,line=4))
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

#-- plot time vs stem respiration

plotBy(R_root~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Root[R]~(gC~d^-1),cex.lab=1.5,line=4))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

# dev.off()
#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# read leaf, stem and root Rdark values for each room  (values at 25C measured at final harvest)

rday<-read.csv("Parameters/great_Resp_leaf_shoot_root_25C.csv")

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
# add Rday estimated form ACi curves

rday_aci<-read.csv("Parameters/rday_aci_fits.csv")

rday<-merge(rday,rday_aci,by="Room")

rday$Tair<-c(18,21.5,25,28.5,32.5,35.5)

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

png(file="output/leaf_r_day_vs_night.png",width=500,height=500)
par(cex.lab=1.5,mar=c(4,5,0.5,0.5),oma=c(0,0,0,0),cex.axis=1,las=1,mfrow=c(1,1))


plot(R_leaf_a~Tair,type="l",col=COL,data=rday,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,2))
points(R_leaf_a~Tair,pch=21,cex=1.5,bg=COL,data=rday,legend=F,lwd=3,ylab="",xlab="",ylim=c(0,2))
adderrorbars(rday$Tair,rday$R_leaf_a,rday$R_leaf_a_se,"updown")


# add day respiration
par(new=T)

plot(Rd.mean~Tair,type="l",col=COL[6],data=rday,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,2))
points(Rd.mean~Tair,pch=21,cex=1.5,bg=COL,data=rday,legend=F,lwd=3,ylab="",xlab="",ylim=c(0,2))
adderrorbars(rday$Tair,rday$Rd.mean,rday$Rd.std.error,"updown")


legend(18,0.3,legend=c("Night","Day"),lty=1,col=c(COL[1],COL[6]),title=expression(R[leaf]~(mu*mol~m^-2~s^-1)),bty="n",lwd=3)

title(ylab=expression(Leaf[R]~(mu*mol~m^-2~s^-1),cex.lab=1.5,line=2))
title(xlab=expression(T[growth]~(degree*C),cex.lab=1.5,line=2))

dev.off()


#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#-- plot time vs total maintenance respiration

plotBy(Rm~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Total[Rm]~(gC~d^-1),cex.lab=1,line=4))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

dev.off()
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

img <- png::readPNG("output/leaf_stem_root_respiration.png")
grid::grid.raster(img)
