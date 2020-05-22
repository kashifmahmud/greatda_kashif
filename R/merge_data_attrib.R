# merge attribution data (from Dushan) and new DA results

listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[2]])
}
summary.param = do.call("rbind", listOfDataFrames)

listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[4]])
}
summary.output = do.call("rbind", listOfDataFrames)

listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[5]])
}
summary.error = do.call("rbind", listOfDataFrames)

if (with.storage==T) {
  listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
  for (i in 1:nlevels(treat.group)) {
    listOfDataFrames[[i]] <- data.frame(result[[i]][[7]])
  }
  summary.storage = do.call("rbind", listOfDataFrames)
}


