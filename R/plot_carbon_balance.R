#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

#-- plot carbon balance

# windows(100,60);

png(file="output/carbon_balance.png",width=800,height=500)
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
barplot(t(toplot), ylim=c(0, 1.5), ylab = "C Partitioning (%GPP)", xlab = expression(Daily~mean~T[air]),  
        col = rainbow(6),legend = c(expression(R[m]),expression(R[g]),expression(C[Leaf]),expression(C[Stem]),expression(C[Root]),expression(C[Storage])), 
        args.legend = list(x = "top", bty = "n",ncol=3))

dev.off()

img <- png::readPNG("output/carbon_balance.png")
grid::grid.raster(img)


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------