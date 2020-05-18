#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------

palette(rev(brewer.pal(6,"Spectral")))
COL=palette()[c(1:6)]
#-- plot time vs leaf respiration

# windows(100,100);
png(file="output/leaf_stem_root_respiration.png",width=500,height=500)
par(cex.lab=1.5,mar=c(2,5,0,0.5),oma=c(2,0,2,0),cex.axis=1,las=1,mfrow=c(2,2))
plotBy(R_leaf~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,0.15))
title(ylab=expression(Leaf[R]~(gC~d^-1),cex.lab=1,line=4))
legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=1.2,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)

# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

#-- plot time vs stem respiration

plotBy(R_stem~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Stem[R]~(gC~d^-1),cex.lab=1,line=4))
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

#-- plot time vs stem respiration

plotBy(R_root~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Root[R]~(gC~d^-1),cex.lab=1.5,line=4))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

#-=---------------------------------------------------------------------------------------------------------------
#-=---------------------------------------------------------------------------------------------------------------


#-- plot time vs total maintenance respiration

plotBy(Rm~Date|Room,type="l",col=COL,data=dCarbon_with_la,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(Total[Rm]~(gC~d^-1),cex.lab=1,line=4))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T,cex=.8)

dev.off()
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------

img <- png::readPNG("output/leaf_stem_root_respiration.png")
grid::grid.raster(img)
