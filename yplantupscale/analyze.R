
# sunny day, february
yruns1 <- readRDS("yruns.rds")

# complete diffuse conditions (same day otherwise)
yruns2 <- readRDS("yruns2.rds")

# Plot simulation results for one of the plants:
plot(yruns1[[1]])
plot(yruns2[[1]])

yres1 <- cbind(summary(yruns1), eucsumm)
yres2 <- cbind(summary(yruns2), eucsumm)


with(yres1, plot(nleavesp, totA / totA0, ylim=c(0.7,1)))
with(yres2, points(nleavesp, totA / totA0, pch=19))

