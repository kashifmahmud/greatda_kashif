# # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
# Try wide priors
#-------------------------------------------------------------------------------------
param.Y <- matrix(c(0.2,0.3,0.4) , nrow=1, ncol=3, byrow=T)
param.af <- matrix(c(0,0.5,1) , nrow=1, ncol=3, byrow=T) # Forcing af not to be negetive
param.as <- matrix(c(0,0.3,1) , nrow=1, ncol=3, byrow=T)
# param.sf <- matrix(c(0,0.0005,0.001) , nrow=1, ncol=3, byrow=T) # All Groups having same sf
# param.sr <- matrix(c(0,0.0005,0.001) , nrow=1, ncol=3, byrow=T) # All Groups having same sr

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
if (no.param > 1) {
  param.Y <- rbind(param.Y, c(-(param.Y[3]-param.Y[1])/2/nrow(data.set), 0, (param.Y[3]-param.Y[1])/2/nrow(data.set)))
  param.af <- rbind(param.af, c(-(param.af[3]-param.af[1])/2/nrow(data.set), 0, (param.af[3]-param.af[1])/2/nrow(data.set)))
  param.as <- rbind(param.as, c(-(param.as[3]-param.as[1])/2/nrow(data.set), 0, (param.as[3]-param.as[1])/2/nrow(data.set)))
  # param.sf <- rbind(param.sf, c(-(param.sf[3]-param.sf[1])/nrow(data.set), 0, (param.sf[3]-param.sf[1])/nrow(data.set)))
  # param.sr <- rbind(param.sr, c(-(param.sr[3]-param.sr[1])/nrow(data.set), 0, (param.sr[3]-param.sr[1])/nrow(data.set)))
}

if (no.param > 2) {
  param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data.set))/(2*nrow(data.set)^2), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data.set))/(2*nrow(data.set)^2)))
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data.set))/(2*nrow(data.set)^2), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data.set))/(2*nrow(data.set)^2)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data.set))/(2*nrow(data.set)^2), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data.set))/(2*nrow(data.set)^2)))
  # param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data.set))/(2*nrow(data.set)^2), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data.set))/(nrow(data.set)^2)))
  # param.sr <- rbind(param.sr, c((param.sr[1,1]-param.sr[1,3]-param.sr[2,3]*nrow(data.set))/(2*nrow(data.set)^2), 0, (param.sr[1,3]-param.sr[1,1]-param.sr[2,1]*nrow(data.set))/(nrow(data.set)^2)))
}
if (no.param > 3) {
  param.Y <- rbind(param.Y, c((param.Y[1,1]-param.Y[1,3]-param.Y[2,3]*nrow(data.set)-param.Y[3,3]*(nrow(data.set)^2))/(2*nrow(data.set)^3), 0, (param.Y[1,3]-param.Y[1,1]-param.Y[2,1]*nrow(data.set)-param.Y[3,1]*(nrow(data.set)^2))/(2*nrow(data.set)^3)))
  param.af <- rbind(param.af, c((param.af[1,1]-param.af[1,3]-param.af[2,3]*nrow(data.set)-param.af[3,3]*(nrow(data.set)^2))/(2*nrow(data.set)^3), 0, (param.af[1,3]-param.af[1,1]-param.af[2,1]*nrow(data.set)-param.af[3,1]*(nrow(data.set)^2))/(2*nrow(data.set)^3)))
  param.as <- rbind(param.as, c((param.as[1,1]-param.as[1,3]-param.as[2,3]*nrow(data.set)-param.as[3,3]*(nrow(data.set)^2))/(2*nrow(data.set)^3), 0, (param.as[1,3]-param.as[1,1]-param.as[2,1]*nrow(data.set)-param.as[3,1]*(nrow(data.set)^2))/(2*nrow(data.set)^3)))
  # param.sf <- rbind(param.sf, c((param.sf[1,1]-param.sf[1,3]-param.sf[2,3]*nrow(data.set)-param.sf[3,3]*(nrow(data.set)^2))/(2*nrow(data.set)^3), 0, (param.sf[1,3]-param.sf[1,1]-param.sf[2,1]*nrow(data.set)-param.sf[3,1]*(nrow(data.set)^2))/(2*nrow(data.set)^3)))
  # param.sr <- rbind(param.sr, c((param.sr[1,1]-param.sr[1,3]-param.sr[2,3]*nrow(data.set)-param.sr[3,3]*(nrow(data.set)^2))/(2*nrow(data.set)^3), 0, (param.sr[1,3]-param.sr[1,1]-param.sr[2,1]*nrow(data.set)-param.sr[3,1]*(nrow(data.set)^2))/(2*nrow(data.set)^3)))
}
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
if (with.storage==F) {
  param = data.frame(param.Y,param.af,param.as)
  names(param) <- c("Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max")
  pMinima <- param[ ,c("Y_min","af_min","as_min")]
  pMaxima <- param[ ,c("Y_max","af_max","as_max")]
  pValues <- param[ ,c("Y","af","as")] # Starting point of the chain
} else { # (with.storage==T) 
  param.k <- matrix(c(0.3,0.6,1) , nrow=1, ncol=3, byrow=T)
  if (no.param > 1) {
    param.k <- rbind(param.k, c(-(param.k[3]-param.k[1])/2/nrow(data.set), 0, (param.k[3]-param.k[1])/2/nrow(data.set)))
  } 
  if (no.param > 2) {
    param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data.set))/(2*nrow(data.set)^2), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data.set))/(2*nrow(data.set)^2)))
  }
  if (no.param > 3) {
    param.k <- rbind(param.k, c((param.k[1,1]-param.k[1,3]-param.k[2,3]*nrow(data.set)-param.k[3,3]*(nrow(data.set)^2))/(2*nrow(data.set)^3), 0, (param.k[1,3]-param.k[1,1]-param.k[2,1]*nrow(data.set)-param.k[3,1]*(nrow(data.set)^2))/(2*nrow(data.set)^3)))
  }
  param = data.frame(param.k,param.Y,param.af,param.as)
  names(param) <- c("k_min","k","k_max","Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max")
  pMinima <- param[ ,c("k_min","Y_min","af_min","as_min")]
  pMaxima <- param[ ,c("k_max","Y_max","af_max","as_max")]
  pValues <- param[ ,c("k","Y","af","as")] # Starting point of the chain
}
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain
#-------------------------------------------------------------------------------------

