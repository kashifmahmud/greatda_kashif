# test SLA with LA and LM
data.attrib <- read.csv("parameters/data_for_attribution_analysis_v2.csv")
# data.attrib_v0 <- read.csv("parameters/data_for_attribution_analysis.csv")
data.attrib$Date = as.Date(data.attrib$Date, format="%d/%m/%Y")
# data.attrib$Date = as.Date(data.attrib$Date)
data.attrib$biomass = data.attrib$Leafmass + data.attrib$Stemmass + data.attrib$Rootmass
data.attrib$Intercept = mean(data.attrib$Intercept)
data.attrib$Slope = mean(data.attrib$Slope)
# data.attrib$SLA = mean(data.attrib$SLA)

##################------------------------------
# take only Room 1 and 4 for attribution analysis
data.attrib.set = subset(data.attrib, Room %in% c(1,4,6))

# average the daily data
data.attrib.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

##################------------------------------
# plot SLA, LA and LM
p1 = ggplot(data = data.attrib.daily, aes(x = Date, y = LA,  group = Room, colour=factor(Room))) +
  geom_point(size=0.01) +
  geom_line(data = data.attrib.daily, aes(x = Date, y = LA,  group = Room, colour=factor(Room))) +
  ylab(expression(Allometry-based ~ LA ~ (m^{2}))) +
  scale_colour_manual(name="", breaks=c("1","4","6"), labels=c("Temp 18.5 C","Temp 25.8 C","Tem 35.5 C"), values=cbPalette[c(3,1,2)]) +
  # annotate("text", x = mean(data.attrib.daily$Date), y = min(data.attrib.daily$SLA), size = font.size-4, label = paste("Triangles = Harvest")) +
  # annotate("text", x = min(data.attrib.daily$Date), y = max(data.attrib.daily$SLA)*0.98, size = font.size, label = paste(title[iteration])) +
  theme_bw() +
  theme(legend.position = c(0.25,0.85)) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0.25, 0.45), units="line")) +
  theme(text = element_text(size=font.size)) +
  theme(legend.key.height=unit(0.65,"line")) +
  # theme(legend.key.width=unit(2,"line")) +
  theme(axis.title = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(vjust=-2))

p2 = ggplot(data = data.attrib.daily, aes(x = Date, y = LM,  group = Room, colour=factor(Room))) +
  geom_point(size=0.01) +
  geom_line(data = data.attrib.daily, aes(x = Date, y = LM,  group = Room, colour=factor(Room))) +
  ylab(expression(DA-modelled ~ LM ~ (gC))) +
  scale_colour_manual(name="", breaks=c("1","4","6"), labels=c("Temp 18.5 C","Temp 25.8 C","Tem 35.5 C"), values=cbPalette[c(3,1,2)]) +
  # annotate("text", x = mean(data.attrib.daily$Date), y = min(data.attrib.daily$SLA), size = font.size-4, label = paste("Triangles = Harvest")) +
  # annotate("text", x = min(data.attrib.daily$Date), y = max(data.attrib.daily$SLA)*0.98, size = font.size, label = paste(title[iteration])) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85)) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0.25, 0.45), units="line")) +
  theme(text = element_text(size=font.size)) +
  theme(legend.key.height=unit(0.65,"line")) +
  # theme(legend.key.width=unit(2,"line")) +
  theme(axis.title = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(vjust=-2)) + guides(color=FALSE)

sla.harvest.all = read.csv("processed_data/sla.harvest.all.csv")
sla.harvest.all = subset(sla.harvest.all, Room %in% as.factor(c("1","4","6")))
sla.harvest.all$Date = as.Date(sla.harvest.all$Date)

p3 = ggplot(data = sla.harvest.all, aes(x = Date, y = sla,  group = Room, colour=factor(Room))) +
  geom_point(position=pd,data = sla.harvest.all, aes(x = Date, y = sla,  group = Room, colour=factor(Room)), size=3, pch=17) +
  geom_errorbar(position=pd,data = sla.harvest.all, aes(ymin=sla-sla_SE, ymax=sla+sla_SE), colour="grey", width=2) +
  geom_point(data = data.attrib.daily, aes(x = Date, y = LA/LM,  group = Room, colour=factor(Room)), size=0.01) +
  geom_line(data = data.attrib.daily, aes(x = Date, y = LA/LM,  group = Room, colour=factor(Room))) +
  ylab(expression(Attribution ~ SLA ~ (m^{2} ~ gC^{-1}))) +
  scale_colour_manual(name="", breaks=c("1","4","6"), labels=c("Temp 18.5 C","Temp 25.8 C","Tem 35.5 C"), values=cbPalette[c(3,1,2)]) +
  annotate("text", x = mean(data.attrib.daily$Date), y = min(data.attrib.daily$SLA), size = font.size-8, label = paste("Triangles = Harvest")) +
  # annotate("text", x = min(data.attrib.daily$Date), y = max(data.attrib.daily$SLA)*0.98, size = font.size, label = paste(title[iteration])) +
  theme_bw() +
  theme(legend.position = c(0.75,0.85)) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank(), plot.margin=unit(c(0.25, 0.25, 0.25, 0.45), units="line")) +
  theme(text = element_text(size=font.size)) +
  theme(legend.key.height=unit(0.65,"line")) +
  # theme(legend.key.width=unit(2,"line")) +
  theme(axis.title = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(vjust=-2)) + guides(color=FALSE)

png("output/1.SLA_great_rm1_4_6.png", units="px", width=1000, height=500, res=100)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



