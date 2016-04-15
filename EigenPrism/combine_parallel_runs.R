# This script combines the results from n.runs parallel runs of
# Brzyski_figure_4_with_EigenPrism.R, and produces some figures.

n.iter.per.run <- 20 # number of iterations per run, see the Ruby script
n.runs <- 200/n.iter.per.run # number of runs, see the Ruby script

FDR.combined      <- rep(0, 7)
FDR.sd.combined   <- rep(0, 7)
pow.combined      <- rep(0, 7)
pow.sd.combined   <- rep(0, 7)
sigma.combined    <- rep(0, 7)
sigma.sd.combined <- rep(0, 7)

for (run.num in 1:n.runs) {
  load(paste0("/home/agossman/github/grpSLOPE_examples/EigenPrism/RData/Brzyski_figure_4_with_EigenPrism_", run.num, ".RData"))
  FDR.combined      <- FDR.combined + n.iter.per.run * FDR
  FDR.sd.combined   <- FDR.sd.combined + (FDR.sd^2 * (n.iter.per.run-1))
  pow.combined      <- pow.combined + n.iter.per.run * pow
  pow.sd.combined   <- pow.sd.combined + (pow.sd^2 * (n.iter.per.run-1))
  sigma.combined    <- sigma.combined + n.iter.per.run * sigma 
  sigma.sd.combined <- sigma.sd.combined + (sigma.sd^2 * (n.iter.per.run-1))
}

FDR.combined      <- FDR.combined / (n.runs*n.iter.per.run) 
FDR.sd.combined   <- sqrt(FDR.sd.combined / (n.runs*n.iter.per.run - 1))
pow.combined      <- pow.combined / (n.runs*n.iter.per.run)
pow.sd.combined   <- sqrt(pow.sd.combined / (n.runs*n.iter.per.run - 1))
sigma.combined    <- sigma.combined / (n.runs*n.iter.per.run)
sigma.sd.combined <- sqrt(sigma.sd.combined / (n.runs*n.iter.per.run - 1))

# Save combined results, and remove everything else from the environment
save(n.runs, n.iter.per.run, n.relevant, fdr, FDR.combined, FDR.sd.combined, 
     pow.combined, pow.sd.combined, sigma.combined, sigma.sd.combined, 
     file="/home/agossman/github/grpSLOPE_examples/EigenPrism/RData/Brzyski_figure_4_with_EigenPrism_combined_results.RData")
rm(list=ls())
load("/home/agossman/github/grpSLOPE_examples/EigenPrism/RData/Brzyski_figure_4_with_EigenPrism_combined_results.RData")

#-------------------------------------------------------
# analog to Figure 4 (a) - q=0.1, Brzyski et. al. (2015)
#-------------------------------------------------------

postscript(file="/home/agossman/github/grpSLOPE_examples/EigenPrism/img/Brzyski_figure_4a_with_EigenPrism.eps", horizontal=FALSE, width=6, height=6)

# plot FDR -------------------------
plot(n.relevant, FDR.combined, ylim=c(0,0.25), type="b", lty=2,
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="Group SLOPE with EigenPrism")

# FDR nominal level
abline(fdr, 0)

legend(9, 0.24, c("gFDR, q=0.1", "Theoretical upper bound"), lty=c(2,1), pch=c(1,NA))

# FDR error bars
FDR.se.combined <- FDR.sd.combined/sqrt(n.runs*n.iter.per.run)
segments(n.relevant, FDR.combined-2*FDR.se.combined, 
         n.relevant, FDR.combined+2*FDR.se.combined, col="blue")
segments(n.relevant-1, FDR.combined-2*FDR.se.combined, 
         n.relevant+1, FDR.combined-2*FDR.se.combined, col="blue")
segments(n.relevant-1, FDR.combined+2*FDR.se.combined, 
         n.relevant+1, FDR.combined+2*FDR.se.combined, col="blue")

dev.off()

#-------------------------------------------------------
# analog to Figure 4 (b) - q=0.1, Brzyski et. al. (2015)
#-------------------------------------------------------

postscript(file="/home/agossman/github/grpSLOPE_examples/EigenPrism/img/Brzyski_figure_4b_with_EigenPrism.eps", horizontal=FALSE, width=6, height=6)

# plot power -------------------------
plot(n.relevant, pow.combined, ylim=c(0,1), type="b", col=1, pch=1,
     xlab="Number of relevant groups", ylab="Estimated power", 
     main="Group SLOPE with Eigen Prism")

# Power error bars
pow.se.combined <- pow.sd.combined/sqrt(n.runs*n.iter.per.run)
segments(n.relevant, pow.combined-2*pow.se.combined, n.relevant, 
         pow.combined+2*pow.se.combined, col="blue")
segments(n.relevant-1, pow.combined-2*pow.se.combined, n.relevant+1, 
         pow.combined-2*pow.se.combined, col="blue")
segments(n.relevant-1, pow.combined+2*pow.se.combined, n.relevant+1, 
         pow.combined+2*pow.se.combined, col="blue")

dev.off()
