# This script combines the results from n.runs parallel runs of
# Brzyski_figure_4_with_EigenPrism.R, and produces some figures.

n.iter <- 20 # number of iterations per run, see the Ruby script
n.runs <- 200/n.iter # number of runs, see the Ruby script

FDR.combined      <- rep(0, 7)
FDR.sd.combined   <- rep(0, 7)
pow.combined      <- rep(0, 7)
pow.sd.combined   <- rep(0, 7)
sigma.combined    <- rep(0, 7)
sigma.sd.combined <- rep(0, 7)

for (i in 1:n.runs) {
  load(paste0("Brzyski_figure_4_with_EigenPrism_", i, ".RData"))
  FDR.combined      <- FDR.combined + n.iter * FDR
  FDR.sd.combined   <- FDR.sd.combined + (FDR.sd^2 * (n.iter-1))
  pow.combined      <- pow.combined + n.iter * pow
  pow.sd.combined   <- pow.sd.combined + (pow.sd^2 * (n.iter-1))
  sigma.combined    <- sigma.combined + n.iter * sigma 
  sigma.sd.combined <- sigma.sd.combined + (sigma.sd^2 * (n.iter-1))
}

FDR.combined      <- FDR.combined / (n.runs*n.iter) 
FDR.sd.combined   <- sqrt(FDR.sd.combined / (n.runs*n.iter - 1))
pow.combined      <- pow.combined / (n.runs*n.iter)
pow.sd.combined   <- sqrt(pow.sd.combined / (n.runs*n.iter - 1))
sigma.combined    <- sigma.combined / (n.runs*n.iter)
sigma.sd.combined <- sqrt(sigma.sd.combined / (n.runs*n.iter - 1))

# Save combined results
save(list=ls(all.names=TRUE), file="/home/agossman/github/grpSLOPE_examples/EigenPrism/RData/Brzyski_figure_4_with_EigenPrism_combined_results.RData")

#-------------------------------------------------------
# analog to Figure 4 (a) - q=0.1, Brzyski et. al. (2015)
#-------------------------------------------------------

postscript(file="Brzyski_figure_4a_with_EigenPrism.eps", horizontal=FALSE, width=6, height=6)

# plot FDR -------------------------
plot(n.relevant, FDR.combined, ylim=c(0,0.25), type="b", lty=2,
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="Group SLOPE with EigenPrism")

# FDR nominal level
abline(fdr, 0)

legend(9, 0.24, c("gFDR, q=0.1", "Theoretical upper bound"), lty=c(2,1), pch=c(1,NA))

# FDR error bars
FDR.se.combined <- FDR.sd.combined/sqrt(n.runs*n.iter)
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

postscript(file="Brzyski_figure_4b_with_EigenPrism.eps", horizontal=FALSE, width=6, height=6)

# plot power -------------------------
plot(n.relevant, pow.combined, ylim=c(0,1), type="b", col=1, pch=1,
     xlab="Number of relevant groups", ylab="Estimated power", 
     main="Group SLOPE with Eigen Prism")

# Power error bars
pow.se.combined <- pow.sd.combined/sqrt(n.iter)
segments(n.relevant, pow.combined-2*pow.se.combined, n.relevant, 
         pow.combined+2*pow.se.combined, col="blue")
segments(n.relevant-1, pow.combined-2*pow.se.combined, n.relevant+1, 
         pow.combined-2*pow.se.combined, col="blue")
segments(n.relevant-1, pow.combined+2*pow.se.combined, n.relevant+1, 
         pow.combined+2*pow.se.combined, col="blue")

dev.off()
