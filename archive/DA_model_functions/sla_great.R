# Calculate SLA from data
sla = data.frame(matrix(ncol = 2, nrow = length(treat.group)))
names(sla) = c("sla","sla.harvest")
sla$Room = treat.group

keeps = c("Date","Room","LA","LA_SE","LM","LM_SE")
data = data.all[ , keeps, drop = FALSE]
data[,c("LM","LM_SE")] = data[,c("LM","LM_SE")] / c1
data[,c("LA","LA_SE")] = data[,c("LA","LA_SE")] * 10000

for (v in 1:length(treat.group)) {
  data.set = subset(data,(Room %in% treat.group[v]))
  # data.set[,c(3:4)] = na.spline(data.set[,c(3:4)])
  data.set = data.set[complete.cases(data.set[, "LM"]), ]
  data.set$sla = data.set$LA / data.set$LM
  data.set$sla_SE = (((data.set$LA_SE/data.set$LA)^2 + (data.set$LM_SE/data.set$LM)^2)^0.5) * data.set$sla
  if (v == 1) {
    data.set.final = data.set
  }
  if (v > 1) {
    data.set.final = rbind(data.set.final,data.set)
  }
}

data.set.final.melt = melt(data.set.final[,c("Date","Room","sla")], id.vars=c("Room","Date"))
names(data.set.final.melt)[4] = c("sla")
sla.melted.error = melt(data.set.final[,c("Date","Room","sla_SE")], id.vars=c("Room","Date"))
names(sla.melted.error)[4] = c("sla_SE")

data.set.final.melt = merge(data.set.final.melt[,c("Date","Room","sla")],sla.melted.error[,c("Date","Room","sla_SE")], by=c("Room","Date"))

plot = list() 
pd <- position_dodge(1) # move the overlapped errorbars horizontally
cbPalette = c("firebrick", "orange", "green3", "yellow3", "#0072B2", "#D55E00", "gray")
font.size = 12
# title = as.character(c("Room#1","Room#2","Room#3","Room#4","Room#5","Room#6"))
  plot[[4]] = ggplot(data = data.set.final.melt, aes(x = Date, y = sla, group = as.factor(Room), colour=as.factor(Room))) + 
    geom_point(position=pd) +
    geom_line(position=pd,data = data.set.final.melt, aes(x = Date, y = sla, group = as.factor(Room), colour=as.factor(Room))) +
    geom_errorbar(position=pd,aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=4) +
    ylab(expression("SLA"~"("*cm^"2"*" "*g^"-1"*")")) + xlab("Date") +
    labs(colour="Room") + 
    # scale_y_continuous(trans = 'log10') +
    # scale_colour_manual(breaks=c("GPP","LM","LM_WM","LM_WM_RM","LM_WM_RM_Rm"), labels=c("GPP","LM","LM + WM","LM + WM + RM",expression("LM + WM + RM + "~R[m])),
    #                     values=cbPalette) +
    # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
    theme_bw() +
    annotate("text", x = mean(data.set.final.melt$Date), y = max(data.set.final$sla)+max(data.set.final$sla_SE), size = font.size-9, label = paste("SLA time course")) +
    # theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(legend.title = element_text(colour="black", size=font.size-2)) +
    theme(legend.text = element_text(colour="black", size = font.size-3)) +
    theme(legend.key.height=unit(0.5,"line")) +
    theme(legend.position = c(0.6,0.08),legend.direction = "horizontal", legend.text.align = 0) +
    theme(legend.key = element_blank()) +
    theme(text = element_text(size=font.size)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
    # theme(plot.title = element_text(hjust = 0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Import final harvest data
final.harvest = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160217-20160224_L3_with date.csv")
final.harvest$Date = as.Date(final.harvest$Date, format = "%d/%m/%Y")
final.harvest = final.harvest[final.harvest$W_treatment %in% as.factor("w"),]
# final.harvest$Room = NULL
final.harvest$Leafarea = rowSums(final.harvest[,c("Leafarea", "Leafarea_sub")], na.rm=T)
final.harvest$Leafmass = rowSums(final.harvest[,c("Leafmass", "Leafmass_sub")], na.rm=T)
# final.harvest$Stemmass = rowSums(final.harvest[,c("Stemmass", "Stemmass_sub")], na.rm=T) # Clarify from John about the experiment procedure
final.harvest$Rootmass = rowSums(final.harvest[,c("Rootmass", "Rootmass_sub")], na.rm=T)
final.harvest[c("Leafarea_sub","Leafmass_sub","Stemmass_sub","Rootmass_sub")] = NULL
final.harvest[ , c("Leafmass", "Stemmass", "Rootmass")] = final.harvest[ , c("Leafmass", "Stemmass", "Rootmass")]/1000 # Unit conversion from mg to g
# final.harvest$Leafarea = final.harvest$Leafarea / 10000  # Unit conversion from cm3 to m3
keeps = c("Room","Prov","Height","Leafarea","Leafmass")
final.harvest.set = final.harvest[ , keeps, drop = FALSE]
final.harvest.set$sla = final.harvest.set$Leafarea / final.harvest.set$Leafmass

final.harvest.sla <- summaryBy(sla ~ Room+Prov, data=final.harvest.set, FUN=c(mean,standard.error))
names(final.harvest.sla)[3:4] = c("sla", "sla_SE")

final.harvest.sla.mean <- summaryBy(sla ~ Room, data=final.harvest.set, FUN=c(mean,standard.error))
names(final.harvest.sla.mean)[2:3] = c("sla", "sla_SE")
final.harvest.sla.mean$Prov = as.factor("Mean")
final.harvest.sla = rbind(final.harvest.sla,final.harvest.sla.mean)

#-------------------------------------------------------------------------------------
# Import weekly height diameter data for 3 months (Diameter is in mm; Height is in cm)
height.dia <- read.csv("Data/GHS39_GREAT_MAIN_HEIGHTDIAMETER_20160108-20160229_L2.csv")
height.dia = height.dia[height.dia$W_treatment %in% as.factor("w"),]
height.dia$D = rowMeans(height.dia[,c("D1", "D2")], na.rm=TRUE)
height.dia = height.dia[complete.cases(height.dia), ]
height.dia$Date = as.Date(height.dia$Date, format = "%d/%m/%Y")

#-------------------------------------------------------------------------------------
# Import initial harvest data
initial.harvest = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160107_L2.csv")
# initial.harvest[ , c("Leafmass", "Stemmass", "Rootmass")] = initial.harvest[ , c("Leafmass", "Stemmass", "Rootmass")]
initial.harvest$Height = initial.harvest$Height/10 # Unit conversion from mm to cm
initial.harvest$Date = as.Date("2016-01-07")
initial.harvest$Room = 0
keeps = c("Room","Prov","Height","Leafarea","Leafmass")
initial.harvest.set = initial.harvest[ , keeps, drop = FALSE]
initial.harvest.set$sla = initial.harvest.set$Leafarea / initial.harvest.set$Leafmass

#-------------------------------------------------------------------------------------
# Import intermediate harvest 1 data
int.harvest.1 = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160129_L2.csv")
int.harvest.1 = int.harvest.1[int.harvest.1$W_treatment %in% as.factor("w"),]
# int.harvest.1[ , c("Leafmass", "Stemmass", "Rootmass")] = int.harvest.1[ , c("Leafmass", "Stemmass", "Rootmass")]
int.harvest.1$Date = as.Date("2016-01-29")
int.harvest.1 = unique(merge(int.harvest.1, height.dia[,c("Room","Pot")]))
keeps = c("Room","Prov","Height","Leafarea","Leafmass")
int.harvest.1.set = int.harvest.1[ , keeps, drop = FALSE]
int.harvest.1.set$sla = int.harvest.1.set$Leafarea / int.harvest.1.set$Leafmass

#-------------------------------------------------------------------------------------
# Import intermediate harvest 2 data
int.harvest.2 = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160210_L2.csv")
int.harvest.2 = int.harvest.2[int.harvest.2$W_treatment %in% as.factor("w"),]
# int.harvest.2[ , c("Leafmass", "Stemmass", "Rootmass")] = int.harvest.2[ , c("Leafmass", "Stemmass", "Rootmass")]
int.harvest.2$Date = as.Date("2016-02-10")
int.harvest.2 = unique(merge(int.harvest.2, height.dia[,c("Room","Pot")]))
keeps = c("Room","Prov","Height","Leafarea","Leafmass")
int.harvest.2.set = int.harvest.2[ , keeps, drop = FALSE]
int.harvest.2.set$sla = int.harvest.2.set$Leafarea / int.harvest.2.set$Leafmass

int.harvest = rbind(int.harvest.1, int.harvest.2)
int.harvest$sla = int.harvest$Leafarea / int.harvest$Leafmass
keeps = c("Room","Prov","sla")
int.harvest.set = int.harvest[ , keeps, drop = FALSE]
int.harvest.sla <- summaryBy(sla ~ Room+Prov, data=int.harvest.set, FUN=c(mean,standard.error))
names(int.harvest.sla)[3:4] = c("sla", "sla_SE")

int.harvest.sla.mean <- summaryBy(sla ~ Room, data=int.harvest.set, FUN=c(mean,standard.error))
names(int.harvest.sla.mean)[2:3] = c("sla", "sla_SE")
int.harvest.sla.mean$Prov = as.factor("Mean")
int.harvest.sla = rbind(int.harvest.sla,int.harvest.sla.mean)

#-------------------------------------------------------------------------------------
# Import the leaf punch datasets
leaf.punch.1 <-read.csv("Data/GHS39_GREAT_MAIN_LEAFAREA-PUNCHES_20160129_L2.csv")
leaf.punch.1$Date <- as.Date("2016-01-29")

leaf.punch.2 <-read.csv("Data/GHS39_GREAT_MAIN_LEAFAREA-PUNCHES_20160209_L2.csv")
leaf.punch.2$Date <- as.Date("2016-02-09")

leaf.punch <- rbind(leaf.punch.1,leaf.punch.2)
leaf.punch = leaf.punch[complete.cases(leaf.punch), ]

leaf.punch$sla <- with(leaf.punch,Puncharea/(Punchmass/1000))         # in cm2 g-1
# leaf.punch$LMA <- with(leaf.punch,(Punchmass/1000)/(Puncharea/10000)) # in g m-2

#dat$prov2 <- dat$prov
#dat$prov <- as.factor(substr(dat$pot,start=1,stop=1)) # overwrite "prov" to be A, B, or C. No Bw or Bd allowed.
leaf.punch$Room <- as.factor(leaf.punch$Room)
# leaf.punch$prov_trt <- as.factor(paste(leaf.punch$Prov,leaf.punch$Room,sep="-"))

#- assign drought treatments
leaf.punch$W_treatment <- factor(leaf.punch$W_treatment,levels=c("w","d"))
leaf.punch = leaf.punch[leaf.punch$W_treatment %in% as.factor("w"),]

leaf.punch.sla <- summaryBy(sla ~ Room+Prov, data=leaf.punch, FUN=c(mean,standard.error))
names(leaf.punch.sla)[3:4] = c("sla", "sla_SE")
#-------------------------------------------------------------------------------------

pd <- position_dodge(0.1) # move the overlapped errorbars horizontally
cbPalette = c("blue", "firebrick", "green3","yellow3")
# title = as.character(c("Room#1","Room#2","Room#3","Room#4","Room#5","Room#6"))
plot[[3]] = ggplot(data = final.harvest.sla, aes(x = Room, y = sla, group = as.factor(Prov), colour=as.factor(Prov))) + 
  geom_point(position=pd) +
  geom_line(position=pd,data = final.harvest.sla, aes(x = Room, y = sla, group = as.factor(Prov), colour=as.factor(Prov))) +
  geom_errorbar(position=pd,aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=0.5) +
  ylab(expression("SLA"~"("*cm^"2"*" "*g^"-1"*")")) + xlab("Room") +
  labs(colour="Provenance") + 
  # scale_y_continuous(trans = 'log10') +
  scale_colour_manual(breaks=c("A","C","B"), labels=c("Cool−origin","Central","Warm−origin"), values=cbPalette) +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
  theme_bw() +
  annotate("text", x = mean(as.numeric(final.harvest.sla$Room)), y = max(final.harvest.sla$sla)+40, size = font.size-9, label = paste("Final harvest data from 17 Feb to 24 Feb")) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=font.size-2)) +
  theme(legend.text = element_text(colour="black", size = font.size-3)) +
  theme(legend.key.height=unit(0.6,"line")) +
  theme(legend.position = c(0.8,0.15),legend.direction = "vertical", legend.text.align = 0) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  # theme(plot.title = element_text(hjust = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
plot[[3]] = plot[[3]] + guides(colour=FALSE)

plot[[2]] = ggplot(data = int.harvest.sla, aes(x = Room, y = sla, group = as.factor(Prov), colour=as.factor(Prov))) + 
  geom_point(position=pd) +
  geom_line(position=pd,data = int.harvest.sla, aes(x = Room, y = sla, group = as.factor(Prov), colour=as.factor(Prov))) +
  geom_errorbar(position=pd,aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=0.5) +
  ylab(expression("SLA"~"("*cm^"2"*" "*g^"-1"*")")) + xlab("Room") +
  labs(colour="Provenance") + 
  # scale_y_continuous(trans = 'log10') +
  scale_colour_manual(breaks=c("A","C","B","Mean"), labels=c("Cool−origin","Central","Warm−origin","Mean"), values=cbPalette) +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
  theme_bw() +
  annotate("text", x = mean(as.numeric(int.harvest.sla$Room)), y = max(int.harvest.sla$sla)+30, size = font.size-9, label = paste("Intermediate harvest data on 29 Jan and 10 Feb")) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=font.size-2)) +
  theme(legend.text = element_text(colour="black", size = font.size-3)) +
  theme(legend.key.height=unit(0.6,"line")) +
  theme(legend.position = c(0.4,0.18),legend.direction = "vertical", legend.text.align = 0) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  # theme(plot.title = element_text(hjust = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# plot[[2]] = plot[[2]] + guides(colour=FALSE)

plot[[1]] = ggplot(data = leaf.punch.sla, aes(x = Room, y = sla, group = as.factor(Prov), colour=as.factor(Prov))) + 
  geom_point(position=pd) +
  geom_line(position=pd,data = leaf.punch.sla, aes(x = Room, y = sla, group = as.factor(Prov), colour=as.factor(Prov))) +
  geom_errorbar(position=pd,aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=0.5) +
  xlab("Room") + ylab(expression("SLA"~"("*cm^"2"*" "*g^"-1"*")")) + 
  labs(colour="Provenance") + 
  # scale_y_continuous(trans = 'log10') +
  scale_colour_manual(breaks=c("A","C","B"), labels=c("Cool−origin","Central","Warm−origin"), values=cbPalette) +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g C plant-1)")) +
  theme_bw() +
  annotate("text", x = mean(as.numeric(leaf.punch.sla$Room)), y = max(leaf.punch.sla$sla)+30, size = font.size-9, label = paste("Leaf Punch data on 29 Jan and 09 Feb")) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=font.size-2)) +
  theme(legend.text = element_text(colour="black", size = font.size-3)) +
  theme(legend.key.height=unit(0.6,"line")) +
  theme(legend.position = c(0.8,0.15),legend.direction = "vertical", legend.text.align = 0) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  # theme(plot.title = element_text(hjust = 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sla.points = final.harvest.sla[final.harvest.sla$Prov %in% as.factor("Mean"),]
sla.points$Date = as.Date("2016-02-21")
pd = position_dodge(1)
plot[[4]] = plot[[4]] + 
  geom_point(position=pd, data = sla.points, aes(x = Date, y = sla, group = as.factor(Room), colour=as.factor(Room)), shape=17, size=2) +
  geom_errorbar(position=pd, data = sla.points, aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=4)
 
int.harvest.1$sla = int.harvest.1$Leafarea / int.harvest.1$Leafmass
int.harvest.1.mean <- summaryBy(sla ~ Room+Date, data=int.harvest.1, FUN=c(mean,standard.error))
names(int.harvest.1.mean)[3:4] = c("sla", "sla_SE")
int.harvest.1.mean$Prov = as.factor("Mean")
# int.harvest.1.mean$Date = as.Date("2016-01-29")
plot[[4]] = plot[[4]] + 
  geom_point(position=pd, data = int.harvest.1.mean, aes(x = Date, y = sla, group = as.factor(Room), colour=as.factor(Room)), shape=17, size=2) +
  geom_errorbar(position=pd, data = int.harvest.1.mean, aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=4)

int.harvest.2$sla = int.harvest.2$Leafarea / int.harvest.2$Leafmass
int.harvest.2.mean <- summaryBy(sla ~ Room+Date, data=int.harvest.2, FUN=c(mean,standard.error))
names(int.harvest.2.mean)[3:4] = c("sla", "sla_SE")
int.harvest.2.mean$Prov = as.factor("Mean")
# int.harvest.1.mean$Date = as.Date("2016-01-29")
plot[[4]] = plot[[4]] + 
  geom_point(position=pd, data = int.harvest.2.mean, aes(x = Date, y = sla, group = as.factor(Room), colour=as.factor(Room)), shape=17, size=2) +
  geom_errorbar(position=pd, data = int.harvest.2.mean, aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=4) +
  annotate("text", x = max(sla.points$Date), y = min(data.set.final.melt$sla)+20, size = font.size-10, label = paste("Triangles = Harvest data \n Circles = Estimated data"))
  

png("output/Figure_sla_great.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
harvest.sla.height = rbind(int.harvest.1.set, int.harvest.2.set,final.harvest.set)

keeps = c("Room","Prov","Height","Leafarea","Leafmass")
int.harvest.1.set = int.harvest.1[ , keeps, drop = FALSE]
int.harvest.1.set$sla = int.harvest.1.set$Leafarea / int.harvest.1.set$Leafmass
int.harvest.1.set$harvest = as.factor("First")
int.harvest.2.set = int.harvest.2[ , keeps, drop = FALSE]
int.harvest.2.set$sla = int.harvest.2.set$Leafarea / int.harvest.2.set$Leafmass
int.harvest.2.set$harvest = as.factor("Second")
final.harvest.set = final.harvest[ , keeps, drop = FALSE]
final.harvest.set$sla = final.harvest.set$Leafarea / final.harvest.set$Leafmass
final.harvest.set$harvest = as.factor("Final")
harvest.sla.height = rbind(int.harvest.1.set, int.harvest.2.set,final.harvest.set)


plot = list() 
plot[[1]] = ggplot() +
  stat_smooth(data = harvest.sla.height, aes(x = Height, y = sla, group = as.factor(Room), colour=as.factor(Room)), method=lm, se = FALSE) +
  geom_point(data = harvest.sla.height, aes(x = Height, y = sla, group = interaction(Room,harvest), colour=as.factor(Room), shape=as.factor(harvest))) +
  labs(colour="Rooms",shape="Harvest") + ylab(expression("SLA"~"("*cm^"2"*" "*g^"-1"*")")) + xlab("Height (cm)") + 
  ggtitle("SLA vs Height with temperature") +
  theme_bw() + theme(legend.position = c(0.9,0.82),legend.key.height=unit(0.5,"line")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/Figure_sla_height_great.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
leaf.punch.sla <- summaryBy(sla ~ Room+Date, data=leaf.punch, FUN=c(mean,standard.error))
names(leaf.punch.sla)[3:4] = c("sla", "sla_SE")
leaf.punch.sla$datatype = as.factor("leafpunch")

int.harvest.sla <- summaryBy(sla ~ Room+Date, data=int.harvest, FUN=c(mean,standard.error))
names(int.harvest.sla)[3:4] = c("sla", "sla_SE")
int.harvest.sla$datatype = as.factor("harvest")

leaf.sla.compare = rbind(leaf.punch.sla,int.harvest.sla)
  
plot = list() 
plot[[1]] = ggplot() +
  geom_point(data = leaf.sla.compare, aes(x = Date, y = sla, group = interaction(Room,datatype), colour=as.factor(Room), shape=as.factor(datatype))) + 
  geom_line(data = leaf.sla.compare, aes(x = Date, y = sla, group = interaction(Room,datatype), colour=as.factor(Room), linetype=as.factor(datatype))) +
  scale_x_date(breaks = date_breaks("11 day"), labels = date_format("%b %d")) +
  labs(colour="Rooms",shape="Data type",linetype="Data type") + ylab(expression("SLA"~"("*cm^"2"*" "*g^"-1"*")")) + xlab("Date") + 
  # ggtitle("SLA vs Height with temperature") +
  theme_bw() + theme(legend.position = c(0.7,0.6),legend.key.height=unit(0.6,"line"),legend.direction = "horizontal") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("output/Figure_sla_leafpunch.png", units="px", width=1600, height=1300, res=220)
print (do.call(grid.arrange,  plot))
dev.off()

