


# function to plot DA model parameters

plot_da_parameters<-function(result){
  
  # plot fitte parameters K, Y, al,as and ar
  
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[2]])
  }
  
  params = do.call("rbind", paramslist)
  
  toplot<-split(params,params$variable)
  
  palette(rev(brewer.pal(6,"Spectral")))
  COL=palette()[c(1:6)]
  
  #- plot Asat
  
pdf(file="output/DA_Figures/parameters.pdf",width=6,height=8)
  
  par(mfrow=c(3,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)
  
  ylabs<-c(expression(k~(gC~g^-1~C~d^-1)),expression(Y),expression(a[f]),expression(a[w]),expression(a[r]),expression(s[f]~(gC~g^-1~C~d^-1)))
  
  
  for(i in 1:length(toplot)){
    
    
    toplot.a<-toplot[[i]]
    
    plotBy(Parameter~Date|treatment,data=toplot.a,legend=F,type="l",las=1,ylim=c(0,1),lwd=3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",col=COL,
           xlab="")
    
    magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
    axis.Date(side=1,at=seq.Date(from=min(toplot.a$Date),to=max(toplot.a$Date),by="week"),labels=T)
    legend("topright",paste("(",letters[i],")",sep=""),bty="n",cex=1.5,text.font=2)
    
    
    title(ylab=ylabs[i],
          outer=F,cex.lab=1.5,line=2)
    
    
    # add standard errors
    
    toplot.b<-split(toplot.a,toplot.a$treatment)
    
    for(i in 1:length(toplot.b)){
      
      toplot.c<-toplot.b[[i]]
      
      polygon(x = c(toplot.c$Date, rev(toplot.c$Date)), y = c((toplot.c$Parameter-toplot.c$Parameter_SD), 
                                                              rev(toplot.c$Parameter+toplot.c$Parameter_SD)),  col = alpha(COL[i],0.5), border = NA)
    }
    
    
  }
  
  #--------------------------------------------------------------------
  # 
  # # add storage to the same plot
  # 
  # nsclist <- vector(mode = "list", length = nlevels(treat.group))
  # for (i in 1:nlevels(treat.group)) {
  #   nsclist[[i]] <- data.frame(result[[i]][[7]])
  # }
  # 
  # nsc = do.call("rbind", nsclist)
  # 
  # # toplot.nsc<-split(params,params$variable)
  # 
  # plotBy(Cstorage.modelled~Date|treatment,data=nsc,legend=F,type="l",las=1,ylim=c(0,1),lwd=3,cex.lab=2,xaxt="n",yaxt="n",
  #        ylab="",col=COL,
  #        xlab="")
  # 
  # magaxis(side=c(2,4),labels=c(1,0),las=1)
  # axis.Date(side=1,at=seq.Date(from=min(toplot.a$Date),to=max(toplot.a$Date),by="week"),labels=T)
  # legend("topright",paste("(",letters[i],")",sep=""),bty="n",cex=1.5,text.font=2)
  # title(ylab=ylabs[i],
  #       outer=F,cex.lab=1.5,line=2)
  # 
  # # add legend
  
  par(mar=c(0,0,0,0))
  plot(1,4,legend=F,type="l",las=1,ylim=c(0,1),lwd=3,cex.lab=2,xaxt="n",yaxt="n",
       ylab="",col="white",
       xlab="",axes=F)
  legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
         pch=16,cex=1.5)
  #--------------------------------------------------------------------
  #--------------------------------------------------------------------
  dev.off()
}

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

# function to plot modelled data

plot_da_data<-function(result){
  
  # plot modelled mass data
  
  
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[4]])
  }
  
  meas_dat <- do.call("rbind", paramslist)
  meas_dat<-subset(meas_dat,!is.na(value))
  
  
  #modelled data
  toplot<-split(meas_dat,meas_dat$variable)
  toplot[4]<-NULL
  
  #----------------------------------------------------------
  
  # # measured data
  # 
  # paramslist <- vector(mode = "list", length = nlevels(treat.group))
  # for (i in 1:nlevels(treat.group)) {
  #   paramslist[[i]] <- data.frame(result[[i]][[3]])
  # }
  # 
  # meas_dat <- do.call("rbind", paramslist)
  # meas_dat_d<-subset(meas_dat,!is.na(value))
  # 
  # toplot_dat<-split(meas_dat_d,meas_dat_d$variable)
  
  #----------------------------------------------------------
  
  # standard errors
  
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[5]])
  }
  
  meas_dat_se <- do.call("rbind", paramslist)
  meas_dat_se<-subset(meas_dat_se,!is.na(value))
  
  toplot_dat_se<-split(meas_dat_se,meas_dat_se$variable)
  toplot_dat_se[4]<-NULL
  
  
  
  # get colours
  
  palette(rev(brewer.pal(6,"Spectral")))
  COL=palette()[c(1:6)]
  

  
  pdf(file="output/DA_Figures/data.pdf",width=6,height=6)
  
  par(mfrow=c(2,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)
  
ylabs<-c(expression(Leaf~mass~(gC~plant^-1)),expression(Stem~mass~(gC~plant^-1)),expression(Root~mass~(gC~plant^-1)),expression(NSC~(gC~plant^-1)))
  
  data<-c("Leafmass","Stemmass","Rootmass")
  serrors<-c("Leafmass_SE","Stemmass_SE","Rootmass_SE")
  
  for(i in 1:length(toplot)){
    
    # plot data
    
    toplot.a<-toplot[[i]]
    
    if(i==1){YLIM<-c(0,3)}
    if(i==2){YLIM<-c(0,2)}
    if(i==3){YLIM<-c(0,1.5)}
    
    
    plotBy(value~Date|treatment,data=toplot.a,legend=F,type="l",lwd=3,las=1,ylim=YLIM,cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",bg=COL,
           xlab="")
    
    magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
    axis.Date(side=1,at=seq.Date(from=min(toplot.a$Date),to=max(toplot.a$Date),by="week"),labels=T)
    legend("topleft",paste("(",letters[i],")",sep=""),bty="n",cex=1.5,text.font=2)
    
    
    title(ylab=ylabs[i],
          outer=F,cex.lab=1.5,line=2)
    
    #--------------------------------------------------------------------


    # add measured data
    
    toplot_data_se<-toplot_dat_se[[i]]
       
   
     
    par(new=T)
    
    
    plotBy(parameter~Date|treatment,data=toplot_data_se,legend=F,pch=16,las=1,ylim=YLIM,cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",col="white",
           xlab="")
    
    adderrorbars(x=toplot_data_se$Date,y=toplot_data_se$parameter,SE=toplot_data_se$value,direction="updown",
                 barlen = 0.01)
    
    par(new=T)
    
    plotBy(parameter~Date|treatment,data=toplot_data_se,legend=F,pch=16,las=1,ylim=YLIM,cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",col=COL[unique(toplot_data_se$treatment)],
           xlab="")
    }
  
  
  # 
  # # add storage to the same plot
  # 
  nsclist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    nsclist[[i]] <- data.frame(result[[i]][[7]])
  }

  nsc = do.call("rbind", nsclist)

  # toplot.nsc<-split(params,params$variable)

  plotBy(Cstorage.modelled~Date|treatment,data=nsc,legend=F,type="l",las=1,ylim=c(0,1.5),lwd=3,cex.lab=2,xaxt="n",yaxt="n",
         ylab="",col=COL,
         xlab="")

  magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
  axis.Date(side=1,at=seq.Date(from=min(toplot.a$Date),to=max(toplot.a$Date),by="week"),labels=T)
  legend("topright",paste("(",letters[4],")",sep=""),bty="n",cex=1.5,text.font=2)
  
  title(ylab=ylabs[4],
        outer=F,cex.lab=1.5,line=2)

  # # add legend
  
    legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
         pch=16,cex=1.2)
  #--------------------------------------------------------------------
  #--------------------------------------------------------------------
  dev.off()
}

# plot_da_data(result=result)

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

# function to plot modelled NSC data



plot_da_nsc<-function(result){
  
  # plot modelled mass data
  
  # mass
  
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[4]])
  }
  
  meas_dat <- do.call("rbind", paramslist)
  
  leaf_mass<-subset(meas_dat,meas_dat$variable=="Mleaf.modelled")[c(1,3,4)]
  names(leaf_mass)[2]<-"leafmass"
  
  stem_mass<-subset(meas_dat,meas_dat$variable=="Mwood.modelled")[c(1,3,4)]
  names(stem_mass)[2]<-"woodmass"
  
  root_mass<-subset(meas_dat,meas_dat$variable=="Mroot.modelled")[c(1,3,4)]
  names(root_mass)[2]<-"rootmass"
  
  mass_dat<-merge(leaf_mass,stem_mass,by=c("Date","treatment"))
  mass_dat2<-merge(mass_dat,root_mass,by=c("Date","treatment"))
  
  
  # NSC
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[7]])
  }
  
  meas_nsc <- do.call("rbind", paramslist)
  
  mass_dat3<-merge(mass_dat2,meas_nsc,by=c("Date","treatment"))
  
  
  
  
  # get colours
  
  palette(rev(brewer.pal(6,"Spectral")))
  COL=palette()[c(1:6)]
  
  
  
  pdf(file="output/DA_Figures/NSC_data.pdf",width=8,height=3.5)
  
  par(mfrow=c(1,3),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)
  
  # leaf
  
    plotBy(Sleaf.modelled/leafmass~Date|treatment,data=mass_dat3,legend=F,type="l",lwd=3,las=1,ylim=c(0,.5),cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",bg=COL[unique(mass_dat3$treatment)],
           xlab="")
    
    magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
    axis.Date(side=1,at=seq.Date(from=min(mass_dat3$Date),to=max(mass_dat3$Date),by="week"),labels=T)
    legend("topright",paste("(",letters[1],")",sep=""),bty="n",cex=1.5,text.font=2)
    
    
    title(ylab=expression(NSC[leaf]~(gC~gC^-1)),
          outer=F,cex.lab=1.5,line=2)
    
    #--------------------------------------------------------------------
    
    # wood
    
    plotBy(Swood.modelled/woodmass~Date|treatment,data=mass_dat3,legend=F,type="l",lwd=3,las=1,ylim=c(0,.5),cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",bg=COL,
           xlab="")
    
    magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
    axis.Date(side=1,at=seq.Date(from=min(mass_dat3$Date),to=max(mass_dat3$Date),by="week"),labels=T)
    legend("topright",paste("(",letters[2],")",sep=""),bty="n",cex=1.5,text.font=2)
    
    
    title(ylab=expression(NSC[wood]~(gC~gC^-1)),
          outer=F,cex.lab=1.5,line=2)
    
    #--------------------------------------------------------------------
    
    # root
    
    plotBy(Sroot.modelled/rootmass~Date|treatment,data=mass_dat3,legend=F,type="l",lwd=3,las=1,ylim=c(0,.5),cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",bg=COL,
           xlab="")
    
    magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
    axis.Date(side=1,at=seq.Date(from=min(mass_dat3$Date),to=max(mass_dat3$Date),by="week"),labels=T)
    legend("topright",paste("(",letters[3],")",sep=""),bty="n",cex=1.5,text.font=2)
    
    
    title(ylab=expression(NSC[root]~(gC~gC^-1)),
          outer=F,cex.lab=1.5,line=2)
    
    
  legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
         pch=16,cex=1.2)
  #--------------------------------------------------------------------
  #--------------------------------------------------------------------
  dev.off()
}

# plot_da_nsc(result=result)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

# function to plot modelled data (log scale)

plot_da_data_log<-function(result){
  
  # plot modelled mass data
  
  
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[4]])
  }
  
  meas_dat <- do.call("rbind", paramslist)
  meas_dat<-subset(meas_dat,!is.na(value))
  
  
  #modelled data
  toplot<-split(meas_dat,meas_dat$variable)
  toplot[4]<-NULL
  
  #----------------------------------------------------------
  
  # # measured data
  # 
  # paramslist <- vector(mode = "list", length = nlevels(treat.group))
  # for (i in 1:nlevels(treat.group)) {
  #   paramslist[[i]] <- data.frame(result[[i]][[3]])
  # }
  # 
  # meas_dat <- do.call("rbind", paramslist)
  # meas_dat_d<-subset(meas_dat,!is.na(value))
  # 
  # toplot_dat<-split(meas_dat_d,meas_dat_d$variable)
  
  #----------------------------------------------------------
  
  # standard errors
  
  paramslist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    paramslist[[i]] <- data.frame(result[[i]][[5]])
  }
  
  meas_dat_se <- do.call("rbind", paramslist)
  meas_dat_se<-subset(meas_dat_se,!is.na(value))
  
  toplot_dat_se<-split(meas_dat_se,meas_dat_se$variable)
  toplot_dat_se[4]<-NULL
  
  
  
  # get colours
  
  palette(rev(brewer.pal(6,"Spectral")))
  COL=palette()[c(1:6)]
  
  
  
  pdf(file="output/DA_Figures/log10_data.pdf",width=6,height=6)
  
  par(mfrow=c(2,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)
  
  ylabs<-c(expression(Leaf~mass~(gC~plant^-1)),expression(Stem~mass~(gC~plant^-1)),expression(Root~mass~(gC~plant^-1)),expression(NSC~(gC~plant^-1)))
  
  data<-c("Leafmass","Stemmass","Rootmass")
  serrors<-c("Leafmass_SE","Stemmass_SE","Rootmass_SE")
  
  for(i in 1:length(toplot)){
    
    # plot data
    
    toplot.a<-toplot[[i]]
    
    if(i==1){YLIM<-c(0,3)}
    if(i==2){YLIM<-c(0,2)}
    if(i==3){YLIM<-c(0,1.5)}
    
    
    plotBy(log10(value)~Date|treatment,data=toplot.a,legend=F,type="l",lwd=3,las=1,cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",bg=COL,
           xlab="")
    
    magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
    axis.Date(side=1,at=seq.Date(from=min(toplot.a$Date),to=max(toplot.a$Date),by="week"),labels=T)
    legend("topleft",paste("(",letters[i],")",sep=""),bty="n",cex=1.5,text.font=2)
    
    
    title(ylab=ylabs[i],
          outer=F,cex.lab=1.5,line=2)
    
    #--------------------------------------------------------------------
    
    
    # add measured data
    
    toplot_data_se<-toplot_dat_se[[i]]
    
    
    
    par(new=T)
    
    
    plotBy(log10(parameter)~Date|treatment,data=toplot_data_se,legend=F,pch=16,las=1,cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",col="white",
           xlab="")
    
    adderrorbars(x=toplot_data_se$Date,y=log10(toplot_data_se$parameter),SE=toplot_data_se$value,direction="updown",
                 barlen = 0.01)
    
    par(new=T)
    
    plotBy(log10(parameter)~Date|treatment,data=toplot_data_se,legend=F,pch=16,las=1,cex=1.3,cex.lab=2,xaxt="n",yaxt="n",
           ylab="",col=COL[unique(toplot_data_se$treatment)],
           xlab="")
  }
  
  
  # 
  # # add storage to the same plot
  # 
  nsclist <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    nsclist[[i]] <- data.frame(result[[i]][[7]])
  }
  
  nsc = do.call("rbind", nsclist)
  
  # toplot.nsc<-split(params,params$variable)
  
  plotBy(log10(Cstorage.modelled)~Date|treatment,data=nsc,legend=F,type="l",las=1,lwd=3,cex.lab=2,xaxt="n",yaxt="n",
         ylab="",col=COL,
         xlab="")
  
  magaxis(side=c(2,4),labels=c(1,0),las=1,ratio=0.4,tcl=0.2)
  axis.Date(side=1,at=seq.Date(from=min(toplot.a$Date),to=max(toplot.a$Date),by="week"),labels=T)
  legend("topright",paste("(",letters[4],")",sep=""),bty="n",cex=1.5,text.font=2)
  
  title(ylab=ylabs[4],
        outer=F,cex.lab=1.5,line=2)
  
  # # add legend
  
  legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
         pch=16,cex=1.2)
  #--------------------------------------------------------------------
  #--------------------------------------------------------------------
  dev.off()
}

# plot_da_data(result=result)

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#' Function for smoothplots of GAMs. 
fitgam <- function(X,Y,dfr, k=-1, R=NULL){
  dfr$Y <- dfr[,Y]
  dfr$X <- dfr[,X]
  if(!is.null(R)){
    dfr$R <- dfr[,R]
    model <- 2
  } else model <- 1
  dfr <- droplevels(dfr)
  
  
  if(model ==1){
    g <- gam(Y ~ s(X, k=k), data=dfr)
  }
  if(model ==2){
    g <- gamm(Y ~ s(X, k=k), random = list(R=~1), data=dfr)
  }
  
  return(g)
}


#' Plot a generalized additive model
#' @param x Variable for X axis (unquoted)
#' @param y Variable for Y axis (unquoted)
#' @param data Dataframe containing x and y
#' @param kgam the \code{k} parameter for smooth terms in gam.
#' @param random An optional random effect (quoted)
#' @param log Whether to add log axes for x or y (but no transformations are done).
#' @param fitoneline Whether to fit only 
smoothplot <- function(x,y,g=NULL,data,
                       fittype=c("gam","lm"),
                       kgam=4,
                       random=NULL,
                       randommethod=c("lmer","aggregate"),
                       log="",
                       fitoneline=FALSE,
                       pointcols=NULL,
                       linecols=NULL, 
                       xlab=NULL, ylab=NULL,
                       polycolor=alpha("lightgrey",0.7),
                       axes=TRUE,
                       ...){
  
  fittype <- match.arg(fittype)
  randommethod <- match.arg(randommethod)
  
  if(!is.null(substitute(g))){
    data$G <- as.factor(eval(substitute(g),data))
  } else {
    fitoneline <- TRUE
    data$G <- 1
  }
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  data <- droplevels(data)
  
  data <- data[!is.na(data$X) & !is.na(data$Y) & !is.na(data$G),]
  
  if(is.null(pointcols))pointcols <- palette()
  if(is.null(linecols))linecols <- palette()
  
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)
  
  # If randommethod = aggregate, average by group and fit simple gam.
  if(!is.null(random) && randommethod == "aggregate"){
    data$R <- data[,random]
    
    data <- summaryBy(. ~ R, FUN=mean, na.rm=TRUE, keep.names=TRUE, data=data,
                      id=~G)
    R <- NULL
  }
  
  
  if(!fitoneline){
    
    d <- split(data, data$G)
    
    if(fittype == "gam"){
      fits <- lapply(d, function(x)try(fitgam("X","Y",x, k=kgam, R=random)))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- lapply(d, function(x)lm(Y ~ X, data=x))
    }
    hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  } else {
    if(fittype == "gam"){
      fits <- list(fitgam("X","Y",data, k=kgam, R=random))
      if(!is.null(random))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- list(lm(Y ~ X, data=data))
    }
    hran <- list(range(data$X, na.rm=TRUE))
    
  }
  
  with(data, plot(X, Y, xaxt="n",yaxt="n", pch=16, col=pointcols[G],
                  xlab=xlab, ylab=ylab, ...))
  
  if(axes){
    if(log=="xy")magaxis(side=1:2, unlog=1:2)
    if(log=="x"){
      magaxis(side=1, unlog=1)
      axis(2)
      box()
    }
    if(log=="y"){
      magaxis(side=2, unlog=2)
      axis(1)
      box()
    }
    if(log==""){
      axis(1)
      axis(2)
      box()
    }
  }
  
  for(i in 1:length(fits)){
    
    if(fittype == "gam"){
      nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
      if(!inherits(fits[[i]], "try-error")){
        p <- predict(fits[[i]],nd,se.fit=TRUE)
        addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=polycolor[i])
        lines(nd$X, p$fit, col=linecols[i], lwd=2)
      }
    }
    if(fittype == "lm"){
      pval <- summary(fits[[i]])$coefficients[2,4]
      LTY <- if(pval < 0.05)1 else 5
      predline(fits[[i]], col=linecols[i], lwd=2, lty=LTY)
    }
  }
  
  return(invisible(fits))
}


addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.7),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}

predline <- function(fit, from=NULL, to=NULL, col=alpha("lightgrey",0.7), ...){
  
  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)
  
  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]
  
  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)
  
  addpoly(newdat[[1]], pred$lwr, pred$upr, col=col)
  
  #ablinepiece(fit, from=from, to=to, ...)
  lines(pred$fit~newdat[,1])
}

#'@title Add a line to a plot
#'@description As \code{abline}, but with \code{from} and \code{to} arguments. 
#'If a fitted linear regression model is used as asn argument, it uses the min and max values of the data used to fit the model.
#'@param a Intercept (optional)
#'@param b Slope (optional)
#'@param reg A fitted linear regression model (output of \code{\link{lm}}).
#'@param from Draw from this X value
#'@param to Draw to this x value
#'@param \dots Further parameters passed to \code{\link{segments}}
#'@export
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
  
}

#----------------------------------------------------------------------------------------------------------------
