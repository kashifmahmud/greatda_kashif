# Matching C balance of the entire experiment considering C inputs and outputs
treat.group = unique(as.factor(data.all$Room)) # Assign all treatments
C.balance = data.frame(matrix(ncol = 7, nrow = length(treat.group)))
names(C.balance) = c("Room","GPP","Rm","biomass","growth.resp","C.output","storage")
C.balance$Room = treat.group
  
for (v in 1:length(treat.group)) {
  data.set = subset(data.all,(Room %in% treat.group[v]))
  # data.set[nrow(data.set),c(10:17)] = data.set[nrow(data.set)-1,c(10:17)]
  data.set[,c(10:17)] = na.spline(data.set[,c(10:17)])
  # plot(data.set$Date, data.set$LM)
  data.set = subset(data.set, Date <= as.Date("2016-02-24"))
  
  C.balance$GPP[v] = sum(data.set$GPP) 
  
  C.balance$Rm[v] = sum(data.set$R_leaf * data.set$LM + data.set$R_wood * data.set$WM + data.set$R_root * data.set$RM)
  C.balance$biomass[v] = (data.set$LM[nrow(data.set)] - data.set$LM[1]) + (data.set$WM[nrow(data.set)] - data.set$WM[1]) + (data.set$RM[nrow(data.set)] - data.set$RM[1])
  C.balance$growth.resp[v] = 0.3 * C.balance$biomass[v]
  C.balance$C.output[v] = C.balance$Rm[v] + C.balance$biomass[v] + C.balance$growth.resp[v]
  C.balance$storage[v] = C.balance$GPP[v] - C.balance$C.output[v]
}

C.balance.fraction = C.balance[, c(3:5,7)]
C.balance.fraction[,] = C.balance.fraction[,] / C.balance[,2] * 100
row.names(C.balance.fraction) <- c("1","2","3","4","5","6")
C.balance.fraction = abs(C.balance.fraction)
  
C.balance = C.balance[,-c(6)]
colnames(C.balance) <- c("Room", "GPP (g C)", "Rm (g C)", "Cs (g C)", "Rg (g C)", "Cn (g C)")
# C.balance = C.balance[,c(10,1,2,3,4,7,5,6,8,9)]
write.csv(C.balance, file = "output/C_partitioning_great.csv", row.names = FALSE)


cbPalette = c("gray", "orange", "skyblue", "green3", "#009E73", "yellow3", "#0072B2", "#D55E00")
png("output/Figure_1a_C_balance_great.png", units="px", width=1200, height=1000, res=200)
par(mfrow = c(1, 1), mar=c(5, 4, 2, 6))
# bb = barplot(as.matrix(t(Ct.fraction.group)), ylim=c(0, 107), ylab = "C Partitioning (%)", xlab = "Treatments (Container size)",  
#         col = rainbow(20),legend = colnames(Ct.fraction.group), 
#         args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))
bb = barplot(as.matrix(t(C.balance.fraction)), ylim=c(0, max(rowSums(C.balance.fraction))+10), ylab = "C Partitioning (%)", xlab = "Container size (L))",  
             col = cbPalette[1:4],legend = c(expression(R["m,tot"]),expression(C[s]),expression(R[g]),expression(C[n])), 
             args.legend = list(x = "topright", bty = "n", inset=c(-0.17, 0)))
# text( bb, Ct.fraction.group[,1]+Ct.fraction.group[,2]+Ct.fraction.group[,3]+Ct.fraction.group[,4]+Ct.fraction.group[,5]+Ct.fraction.group[,6]+Ct.fraction.group[,7]-1, labels = round(Ct.group[,9],1), cex=.9)
text( bb, rowSums(C.balance.fraction)+0.5, labels = round(C.balance[,2],1), pos = 3, cex=1, col="red")

dev.off()


#-------------------------------------------------------------------------------------
# Plot time course of GPP and all components of C utilised (LM, WM, RM, Rm) 
for (v in 1:length(treat.group)) {
  data.set = subset(data.all,(Room %in% treat.group[v]))
  data.set[,c(10:17)] = na.spline(data.set[,c(10:17)])
  data.set$Rm = data.set$R_leaf*data.set$LM + data.set$R_wood*data.set$WM + data.set$R_root*data.set$RM
  data.set[,c("GPP","Rm")] = cumsum(data.set[,c("GPP","Rm")])
  data.set$LM = data.set$LM - data.set$LM[1]
  data.set$WM = data.set$WM - data.set$WM[1]
  data.set$RM = data.set$RM - data.set$RM[1]
  # data.set$LM[1] = data.set$WM[1] = data.set$RM[1] = 0.001
  
  if (v == 1) {
    data.set.final = data.set
  }
  if (v > 1) {
    data.set.final = rbind(data.set.final,data.set)
  }
}
keeps = c("Date","Room","GPP","LM","WM","RM","Rm")
data.set.final = data.set.final[ , keeps, drop = FALSE]
data.set.final$LM_WM = data.set.final$LM + data.set.final$WM
data.set.final$LM_WM_RM = data.set.final$LM_WM + data.set.final$RM
data.set.final$LM_WM_RM_Rm = data.set.final$LM_WM_RM + data.set.final$Rm
data.set.final[data.set.final <= 0] <- 0

data.set.final.melt = melt(data.set.final[,c("Date","Room","GPP","LM","LM_WM","LM_WM_RM","LM_WM_RM_Rm")], id.vars=c("Room","Date"))

plot = list() 
pd <- position_dodge(0) # move the overlapped errorbars horizontally
cbPalette = c("firebrick", "orange", "green3", "yellow3", "#0072B2", "#D55E00", "gray")
title = as.character(c("Room#1","Room#2","Room#3","Room#4","Room#5","Room#6"))
for (v in 1:length(treat.group)) {
  data.set.final.melt.set = subset(data.set.final.melt,(Room %in% treat.group[v]))
  plot[[v]] = ggplot() + 
    geom_line(position=pd,data = data.set.final.melt.set, aes(x = Date, y = value, group = variable, colour=variable)) +
    ylab(paste("Cumulative C (g C)")) + xlab("Date") +
    labs(colour="C components") + 
    # scale_y_continuous(trans = 'log10') +
    scale_colour_manual(breaks=c("GPP","LM","LM_WM","LM_WM_RM","LM_WM_RM_Rm"), labels=c("GPP","LM","LM + WM","LM + WM + RM",expression("LM + WM + RM + "~R[m])),
                        values=cbPalette) +
    # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
    theme_bw() +
    annotate("text", x = max(data.set.final.melt.set$Date)-3, y = min(data.set.final.melt.set$value), size = font.size-9, label = paste(title[v])) +
    # theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(legend.title = element_text(colour="black", size=font.size-4)) +
    theme(legend.text = element_text(colour="black", size = font.size-5)) +
    theme(legend.key.height=unit(0.6,"line")) +
    theme(legend.position = c(0.25,0.7),legend.direction = "vertical", legend.text.align = 0) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    # theme(plot.title = element_text(hjust = 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  # plot[[v]] = ggplot() + 
  #   geom_line(position=pd,data = data.set.final.melt.set, aes(x = Date, y = value, group = variable, colour=variable)) +
  #   ylab(paste("Cumulative C (g C)")) + xlab("Date") +
  #   labs(colour="C components") + 
  #   scale_y_continuous(trans = 'log10') +
  #   scale_colour_manual(breaks=c("GPP","LM","LM_WM","LM_WM_RM","LM_WM_RM_Rm"), labels=c("GPP","LM","LM + WM","LM + WM + RM",expression("LM + WM + RM + "~R[m])),
  #                       values=cbPalette) +
  #   # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
  #   theme_bw() +
  #   annotate("text", x = min(data.set.final.melt.set$Date)+5, y = max(data.set.final.melt.set$value)-2, size = font.size-7, label = paste(title[v])) +
  #   # theme(plot.title = element_text(size = 20, face = "bold")) +
  #   theme(legend.title = element_text(colour="black", size=font.size-2)) +
  #   theme(legend.text = element_text(colour="black", size = font.size-3)) +
  #   theme(legend.key.height=unit(0.6,"line")) +
  #   theme(legend.position = c(0.75,0.28),legend.direction = "vertical", legend.text.align = 0) +
  #   theme(legend.key = element_blank()) +
  #   theme(text = element_text(size=font.size)) +
  #   theme(axis.title.x = element_blank()) +
  #   theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  #   # theme(plot.title = element_text(hjust = 0)) +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
   
  if (v>1) {
    plot[[v]] = plot[[v]] + guides(colour=FALSE)
  }
}

png("output/Figure_C_timecourse.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()



#-------------------------------------------------------------------------------------

