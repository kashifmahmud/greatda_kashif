# test run 1

# seed 18
source("R/DA_model_functions/DA_Head.R")
result_r1 <- clusterMap(cluster, mcmc.great, 
                     treat.group=treat.group,
                     MoreArgs=list(chainLength=5000,with.storage=T, model.comparison=F, model.optimization=F, 
                                   no.param.per.var=2))



# seed 5
# source("R/DA_model_functions/DA_Head.R")

result_r2 <- clusterMap(cluster, mcmc.great, 
                        treat.group=treat.group,
                        MoreArgs=list(chainLength=5000,with.storage=T, model.comparison=F, model.optimization=F, 
                                      no.param.per.var=2))

# seed 50
# source("R/DA_model_functions/DA_Head.R")

result_r3 <- clusterMap(cluster, mcmc.great, 
                        treat.group=treat.group,
                        MoreArgs=list(chainLength=5000,with.storage=T, model.comparison=F, model.optimization=F, 
                                      no.param.per.var=2))



source("R/plot_DA_results.R")
# seed 18 bic 16

plot_da_nsc(result=result_r1)

# seed 5 bic 25
plot_da_nsc(result=result_r2)


# seed 50 bic 22
plot_da_nsc(result=result_r3)
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

# plot fitte parameters K

# Run 1
paramslist <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  paramslist[[i]] <- data.frame(result_r1[[i]][[2]])
}

params_r1 = do.call("rbind", paramslist)

toplot_r1<-split(params_r1,params_r1$variable)

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#Run 2

paramslist <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  paramslist[[i]] <- data.frame(result_r2[[i]][[2]])
}

params_r2 = do.call("rbind", paramslist)

toplot_r2<-split(params_r2,params_r2$variable)

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#Run 3

paramslist <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  paramslist[[i]] <- data.frame(result_r3[[i]][[2]])
}

params_r3 = do.call("rbind", paramslist)

toplot_r3<-split(params_r3,params_r3$variable)

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------


# plot K

palette(rev(brewer.pal(6,"Spectral")))
COL=palette()[c(1:6)]

plot_parameters<-function(data,YLIM=c(0,1),yno=1,letter=1){
  
  ylabs<-c(expression(k~(gC~g^-1~C~d^-1)),expression(Y),expression(a[f]),expression(a[w]),expression(a[r]),expression(s[f]~(gC~g^-1~C~d^-1)))

  plotBy(Parameter~Date|treatment,data=data,legend=F,type="l",las=1,ylim=YLIM,lwd=3,cex.lab=2,xaxt="n",yaxt="n",
         ylab="",col=COL,
         xlab="")
  
  magaxis(side=c(2,4),labels=c(1,0),las=1)
  axis.Date(side=1,at=seq.Date(from=min(data$Date),to=max(data$Date),by="week"),labels=T)
  legend("topright",letter,bty="n",cex=1.5,text.font=2)
  
  
  
  toplot.b<-split(data,data$treatment)
  
  for(i in 1:length(toplot.b)){
    
    toplot.c<-toplot.b[[i]]
    
    polygon(x = c(toplot.c$Date, rev(toplot.c$Date)), y = c((toplot.c$Parameter-toplot.c$Parameter_SD), 
                                                            rev(toplot.c$Parameter+toplot.c$Parameter_SD)),  col = alpha(COL[i],0.5), border = NA)
  }
  
  title(ylab=ylabs[yno],
        outer=F,cex.lab=1.5,line=2)
  
  
  
}


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# plot K
pdf(file="output/DA_Figures/parameters_K_test_run.pdf",width=8,height=3.5)
par(mfrow=c(1,3),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

plot_parameters(toplot_r1[[1]],yno=1,letter="test-1")
plot_parameters(toplot_r2[[1]],yno=1,letter="test-2")
plot_parameters(toplot_r3[[1]],yno=1,letter="test-3")

dev.off()

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# plot Y

pdf(file="output/DA_Figures/parameters_Y_test_run.pdf",width=8,height=3.5)
par(mfrow=c(1,3),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

plot_parameters(toplot_r1[[2]],yno=2,letter="test-1",YLIM=c(0,0.5))
plot_parameters(toplot_r2[[2]],yno=2,letter="test-2",YLIM=c(0,0.5))
plot_parameters(toplot_r3[[2]],yno=2,letter="test-3",YLIM=c(0,0.5))
legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2)
dev.off()

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# plot al

pdf(file="output/DA_Figures/parameters_af_test_run.pdf",width=8,height=3.5)
par(mfrow=c(1,3),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

plot_parameters(toplot_r1[[3]],yno=3,letter="test-1",YLIM=c(0,1))
plot_parameters(toplot_r2[[3]],yno=3,letter="test-2",YLIM=c(0,1))
plot_parameters(toplot_r3[[3]],yno=3,letter="test-3",YLIM=c(0,1))
legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2)
dev.off()


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# plot aw

pdf(file="output/DA_Figures/parameters_aw_test_run.pdf",width=8,height=3.5)
par(mfrow=c(1,3),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

plot_parameters(toplot_r1[[4]],yno=4,letter="test-1",YLIM=c(0,1))
plot_parameters(toplot_r2[[4]],yno=4,letter="test-2",YLIM=c(0,1))
plot_parameters(toplot_r3[[4]],yno=4,letter="test-3",YLIM=c(0,1))
legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2)
dev.off()

#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# plot ar

pdf(file="output/DA_Figures/parameters_ar_test_run.pdf",width=8,height=3.5)
par(mfrow=c(1,3),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

plot_parameters(toplot_r1[[5]],yno=5,letter="test-1",YLIM=c(0,1))
plot_parameters(toplot_r2[[5]],yno=5,letter="test-2",YLIM=c(0,1))
plot_parameters(toplot_r3[[5]],yno=5,letter="test-3",YLIM=c(0,1))
legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2)

dev.off()


#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

treat.group = as.factor(c("1","4","6")) # Assign few treatments to check the results



result_r4 <- clusterMap(cluster, mcmc.great, 
                        treat.group=treat.group,
                        MoreArgs=list(chainLength=15000,with.storage=T, model.comparison=F, model.optimization=F, 
                                      no.param.per.var=2))


# Run 1
paramslist <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  paramslist[[i]] <- data.frame(result_r4[[i]][[2]])
}

params_r4 = do.call("rbind", paramslist)

toplot_r4<-split(params_r4,params_r4$variable)

# result_r5 <- clusterMap(cluster, mcmc.great, 
#                         treat.group=treat.group,
#                         MoreArgs=list(chainLength=15000,with.storage=T, model.comparison=F, model.optimization=F, 
#                                       no.param.per.var=2))
# 
# 
# result_r6 <- clusterMap(cluster, mcmc.great, 
#                         treat.group=treat.group,
#                         MoreArgs=list(chainLength=15000,with.storage=T, model.comparison=F, model.optimization=F, 
#                                       no.param.per.var=2))


#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

# plot K
pdf(file="output/DA_Figures/parameters_K_test_run_room.pdf",width=8,height=3.5)
par(mfrow=c(1,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

COL=palette()[c(1:6)]
plot_parameters(data=toplot_r1[[1]],yno=1,letter="Trt-6")

COL=palette()[c(1,4,6)]
plot_parameters(data=toplot_r4[[1]],yno=1,letter="Trt-3")


COL=palette()[c(1:6)]
legend("bottomleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2,ncol=3)

dev.off()

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# plot Y
pdf(file="output/DA_Figures/parameters_Y_test_run_room.pdf",width=8,height=3.5)
par(mfrow=c(1,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

COL=palette()[c(1:6)]
plot_parameters(data=toplot_r1[[2]],yno=2,letter="Trt-6")

COL=palette()[c(1,4,6)]
plot_parameters(data=toplot_r4[[2]],yno=2,letter="Trt-3")


COL=palette()[c(1:6)]
legend("topleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2,ncol=3)

dev.off()

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# plot al

pdf(file="output/DA_Figures/parameters_af_test_run_room.pdf",width=8,height=3.5)
par(mfrow=c(1,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

COL=palette()[c(1:6)]
plot_parameters(data=toplot_r1[[3]],yno=3,letter="Trt-6")

COL=palette()[c(1,4,6)]
plot_parameters(data=toplot_r4[[3]],yno=3,letter="Trt-3")


COL=palette()[c(1:6)]
legend("bottom",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1.2,ncol=3)

dev.off()

