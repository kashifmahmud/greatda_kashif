#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

# Interpolate leaf area, leaf mass, stem mass and root mass

# harvest_dat<-read.csv(paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/harvest_data.csv"))


harvest_dat<-read.csv("Parameters/modelled_data.csv") # leaf area, leaf mass, stem mass and root mass from allometric models
harvest_dat$Date<-as.Date(harvest_dat$Date,format="%Y-%m-%d")


with(harvest_dat,plot(Date,Leafarea,col=Room,pch=16))



datesout <- seq(as.Date("2016/1/8"), as.Date("2016/2/29"), by="1 day") # get a sequance of dated 


# function to interpolate variables

interp_var <- function(d, datesout,var){
  
  sp <- spline(x=d$Date, y=d[[var]], xout=datesout) #smooth curve
  # sp <- approx(x=d$Date, y=d[[var]], xout=datesout) # linear interpolation
  
  out<-data.frame(Date=datesout, Room=unique(d$Room), var=sp$y)
  return(out)
}


#- leaf area
leaf_area <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Leafarea"))
names(leaf_area)[3]<-"leafarea_p"

#- leaf mass
Leafmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Leafmass"))
names(Leafmass)[3]<-"Leafmass"

#- stem mass

Stemmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Stemmass"))
names(Stemmass)[3]<-"Stemmass"

#- root mass
Rootmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Rootmass"))
names(Rootmass)[3]<-"Rootmass"


seedling_mass<-cbind(leaf_area,Leafmass[3],Stemmass[3],Rootmass[3])


write.csv(seedling_mass,file = "Parameters/great_leaf_area_predicted.csv",row.names=FALSE,sep=",")



#--Plot interpolated leaf area

COL <- rev(brewer.pal(6,"Spectral"))
windows(100,100);

par(cex.lab=1.5,mar=c(0,5,2,0.5),oma=c(2,0,0,0),cex.axis=1,las=1,mfrow=c(2,2))
plotBy(leafarea_p~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Leafarea,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
# with(harvest_dat_h,points(Date,Leafarea,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Leafarea,SE=harvest_dat_h$Leafarea_SE,direction="updown")

# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
#legend("topright",letters[1],bty="n",cex=1.2)

title(ylab=expression(Leaf[area]~(m^-2)),cex.lab=1.5,line=2.5)

legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


#--Plot leaf mass

plotBy(Leafmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Leafmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Leaf[mass]~(gC),cex.lab=1.5,line=1))

# with(harvest_dat_h,points(Date,Leafmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Leafmass,SE=harvest_dat_h$Leafmass_SE,direction="updown")


#--Plot stem mass

plotBy(Stemmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Stemmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Stem[mass]~(gC),cex.lab=1.5,line=2.5))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

# with(harvest_dat_h,points(Date,Stemmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Stemmass,SE=harvest_dat_h$Stemmass_SE,direction="updown")

#--Plot Root mass

plotBy(Rootmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Rootmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Root[mass]~(gC),cex.lab=1.5,line=3))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

# with(harvest_dat_h,points(Date,Rootmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Rootmass,SE=harvest_dat_h$Rootmass_SE,direction="updown")

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

