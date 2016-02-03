###############################################################################
###############################################################################
###############################################################################
# This code aims to reproduce the data generation, analysis and plots from
# D. Brzyski, W. Su, M. Bogdan (2015), "Group SLOPE - adaptive selection of 
# groups of predictors" (http://arxiv.org/abs/1511.09078)
###############################################################################
###############################################################################
###############################################################################

library(grpSLOPE)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=3)

set.seed(1)

############
# Figure 4
############

fdr <- 0.1
n.iter <- 200

n.group <- 1000
group.length <- rbinom(n.group, 1000, 0.008)
group <- c()
for (i in 1:n.group) {
  group <- c(group, rep(i, group.length[i]))
}
group.id <- getGroupID(group)

n <- 5000
p <- length(group)
X <- matrix(rnorm(n*p, sd=1/sqrt(n)), n, p)
X <- scale(X, center=TRUE, scale=FALSE)
X <- apply(X, 2, function(x) x/sqrt(sum(x^2)) )

Bfun <- function(l) {
  B <- 4*log(n.group) / (1 - n.group^(-2/l)) - l
  return(sqrt(B))
}
signal.strength <- sum(Bfun(group.length)) / n.group

n.relevant <- floor(seq(1, 60, length=7))

FDR    <- rep(NA, length(n.relevant))
FDR.sd <- rep(NA, length(n.relevant))
pow    <- rep(NA, length(n.relevant))
pow.sd <- rep(NA, length(n.relevant))

one.iteration <- function(n.signif){
  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) {
    X1 <- apply(as.matrix(X[ , group.id[[j]]]), 1, sum)
    b[group.id[[j]]] <- (signal.strength / norm(as.matrix(X1), "f")) * rep(1, group.length[j])
  }

  # generate the response vector
  y <- X %*% b + rnorm(n, sd=1)

  # get Group SLOPE solution
  b.grpSLOPE <- grpSLOPE(X=X, y=y, group=group, fdr=fdr,
                         lambda="chiMean", verbose=FALSE,
                         orthogonalize=FALSE, normalize=FALSE)

  # FDR and power
  n.selected <- length(b.grpSLOPE$selected)
  true.relevant <- names(group.id)[ind.relevant]
  truepos <- intersect(b.grpSLOPE$selected, true.relevant)
  n.truepos <- length(truepos)
  n.falsepos <- n.selected - n.truepos
  FDR <- n.falsepos / max(1, n.selected)
  pow <- n.truepos / length(true.relevant)

  return(list(FDR=FDR, pow=pow))
}

for (k in 1:length(n.relevant)) {
  parallel.results <- foreach(i=1:n.iter) %dopar% {
    one.iteration(n.relevant[k])
  }

  FDR.vec  <- rep(NA, n.iter)
  pow.vec  <- rep(NA, n.iter)

  for (j in 1:n.iter) {
    FDR.vec[j]  <- parallel.results[[j]]$FDR
    pow.vec[j]  <- parallel.results[[j]]$pow
  }

  FDR[k] <- mean(FDR.vec)
  FDR.sd[k] <- sd(FDR.vec)
  pow[k] <- mean(pow.vec)
  pow.sd[k] <- sd(pow.vec)
  
  print(paste(k, "sparsity levels completed"))
}

#####################################################
# Figure 4 (a) - q=0.1, Brzyski et. al. (2015)
#####################################################

# plot FDR -------------------------
plot(n.relevant, FDR, ylim=c(0,0.25), type="b", lty=2,
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="corrected lambdas + sigma estimation")

# FDR nominal level
abline(fdr, 0)

legend(9, 0.24, c("gFDR, q=0.1", "Theoretical upper bound"), lty=c(2,1), pch=c(1,NA))

# FDR error bars
FDR.se <- FDR.sd/sqrt(n.iter)
segments(n.relevant, FDR-2*FDR.se, n.relevant, FDR+2*FDR.se, col="blue")
segments(n.relevant-1, FDR-2*FDR.se, n.relevant+1, FDR-2*FDR.se, col="blue")
segments(n.relevant-1, FDR+2*FDR.se, n.relevant+1, FDR+2*FDR.se, col="blue")

####################################################################
# Figure 4 (b) - q=0.1, Brzyski et. al. (2015)
####################################################################

# plot power -------------------------
plot(n.relevant, pow, ylim=c(0,1), type="b", col=1, pch=1,
     xlab="Number of relevant groups", ylab="Estimated power", 
     main="Power")

# Power error bars
pow.se <- pow.sd/sqrt(n.iter)
segments(n.relevant, pow-2*pow.se, n.relevant, pow+2*pow.se, col="blue")
segments(n.relevant-1, pow-2*pow.se, n.relevant+1, pow-2*pow.se, col="blue")
segments(n.relevant-1, pow+2*pow.se, n.relevant+1, pow+2*pow.se, col="blue")

####################################################################
# Figure 4 (c) - q=0.1, Brzyski et. al. (2015)
####################################################################

hist(group.length, xlab = "Group size",
     main = "Histogram of group sizes")
