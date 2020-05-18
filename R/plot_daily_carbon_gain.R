

COL <- rev(brewer.pal(6,"Spectral"))
# #---------------------------------------------------------------------------------------------------------------
# #--- plot daily carbon gain per m2
# 
# windows(60,40);
png(file="output/daily_carbon.png",width=700,height=500)

par(mfrow=c(2,1),cex.lab=1.5,mar=c(3,5,0,0),oma=c(2,2,2,2),cex.axis=1,las=1)
plotBy(gCarbon~Date|Room,type="l",col=COL,data=dCarbon,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,15))


title(ylab=expression(C[Day]~(gC~m^-2~d^-1),cex.lab=.8))
# legend("topleft",legend=unique(dCarbon$Room),fill=COL,cex=1.3,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
 
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------

# #--- plot daily carbon gain per plant
# 
# windows(60,40);
# par(cex.lab=1.5,mar=c(0,7,2,0.5),oma=c(3,2,0,3),cex.axis=1,las=1,mfrow=c(1,1))
plotBy(gCarbon_la~Date|Room,type="l",col=COL,data=dCarbon,legend=F,lwd=3,ylab="",xlab="",cex=1,ylim=c(0,1))


title(ylab=expression(C[Day]~(gC~plant^-1),cex.lab=.8))
legend("topleft",legend=unique(dCarbon$Room),fill=COL,cex=1.2,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

 dev.off()
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------


img <- png::readPNG("output/daily_carbon.png")
grid::grid.raster(img)

