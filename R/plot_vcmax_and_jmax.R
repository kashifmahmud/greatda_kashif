# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------

# # # -- plot time vs estimated Vcmax and Jmax
# # 
COL <- rev(brewer.pal(6,"Spectral"))
# windows(80,100);

png(file="output/vcmax_jmax.png",width=600,height=500)

par(cex.lab=1.5,mar=c(2,6,0.5,0.5),oma=c(1,1,0,0),cex.axis=1.5,las=1,mfrow=c(2,1))

plotBy(Vcmax25~Date|Room,type="l",col=COL,data=metwithphy,legend=F,lwd=3,ylab="",xlab="",ylim=c(0,150))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

title(ylab=expression(V[cmax(25)]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=3.5)
abline(v=as.Date("2016-02-26"),lty=3)
abline(v=as.Date("2016-02-03"),lty=3)

plotBy(Jmax25~Date|Room,type="l",col=COL,data=metwithphy,legend=F,lwd=3,ylab="",xlab="",ylim=c(0,250))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

legend("bottomleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


title(ylab=expression(J[max(25)]~(mu*mol~m^-2~s^-1)),cex.lab=1.5,line=3.5)
title(xlab=expression(Daily~mean~T[air]~(degree*C)),cex.lab=1.3)
abline(v=as.Date("2016-02-26"),lty=3)
abline(v=as.Date("2016-02-03"),lty=3)

dev.off()
# # #---------------------------------------------------------------------------------------------------------------------------------
# # #---------------------------------------------------------------------------------------------------------------------------------

img <- png::readPNG("output/vcmax_jmax.png")
grid::grid.raster(img)
