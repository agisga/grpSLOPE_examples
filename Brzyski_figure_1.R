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
registerDoParallel(cores=6)

############
# Figure 1
############

fdr <- 0.1
n.iter <- 300

p <- 5000
X <- diag(rep(1,p))
n.group <- 1000
group <- c(rep(1:200, each=3),
           rep(201:400, each=4),
           rep(401:600, each=5),
           rep(601:800, each=6),
           rep(801:1000, each=7))
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)

Bfun <- function(l) {
  sqrt(4*log(n.group) * (1 - n.group^(-2/l)) - l)
}
a <- sum(Bfun(group.length)) / sum(sqrt(group.length))

n.relevant <- floor(seq(1, 250, length=11))

FDR.max    <- rep(NA, length(n.relevant))
FDR.max.sd <- rep(NA, length(n.relevant))
pow.max    <- rep(NA, length(n.relevant))
pow.max.sd <- rep(NA, length(n.relevant))

FDR.mean    <- rep(NA, length(n.relevant))
FDR.mean.sd <- rep(NA, length(n.relevant))
pow.mean    <- rep(NA, length(n.relevant))
pow.mean.sd <- rep(NA, length(n.relevant))

one.iteration <- function(n.signif){
  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) { b[group.id[[j]]] <-a }

  # generate the response vector
  y <- X %*% b + rnorm(p, sd=1)

  # get Group SLOPE solution
  b.grpSLOPE.max <- grpSLOPE(X=X, y=y, group=group, fdr=fdr,
                             lambda="chiOrthoMax", sigma=1, verbose=FALSE,
                             orthogonalize=FALSE, normalize=FALSE)
  b.grpSLOPE.mean <- grpSLOPE(X=X, y=y, group=group, fdr=fdr,
                              lambda="chiOrthoMean", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)

  # FDR and power
  true.relevant <- names(group.id)[ind.relevant]
  truepos <- length(intersect(b.grpSLOPE.max$selected, true.relevant))
  falsepos <- length(b.grpSLOPE.max$selected) - truepos
  FDR.max <- falsepos / max(1, length(b.grpSLOPE.max$selected))
  pow.max <- truepos / length(true.relevant)

  true.relevant <- names(group.id)[ind.relevant]
  truepos <- length(intersect(b.grpSLOPE.mean$selected, true.relevant))
  falsepos <- length(b.grpSLOPE.mean$selected) - truepos
  FDR.mean <- falsepos / max(1, length(b.grpSLOPE.mean$selected))
  pow.mean <- truepos / length(true.relevant)

  return(list(FDR.max=FDR.max, pow.max=pow.max,
              FDR.mean=FDR.mean, pow.mean=pow.mean))
}

for (k in 1:length(n.relevant)) {
  parallel.results <- foreach(i=1:n.iter) %dopar% {
    one.iteration(n.relevant[k])
  }

  FDR.max.vec  <- rep(NA, n.iter)
  pow.max.vec  <- rep(NA, n.iter)
  FDR.mean.vec <- rep(NA, n.iter)
  pow.mean.vec <- rep(NA, n.iter)

  for (j in 1:n.iter) {
    FDR.max.vec[j]  <- parallel.results[[j]]$FDR.max
    pow.max.vec[j]  <- parallel.results[[j]]$pow.max
    FDR.mean.vec[j] <- parallel.results[[j]]$FDR.mean
    pow.mean.vec[j] <- parallel.results[[j]]$pow.mean
  }

  FDR.max[k] <- mean(FDR.max.vec)
  FDR.max.sd[k] <- sd(FDR.max.vec)
  pow.max[k] <- mean(pow.max.vec)
  pow.max.sd[k] <- sd(pow.max.vec)

  FDR.mean[k] <- mean(FDR.mean.vec)
  FDR.mean.sd[k] <- sd(FDR.mean.vec)
  pow.mean[k] <- mean(pow.mean.vec)
  pow.mean.sd[k] <- sd(pow.mean.vec)
  
  print(paste(k, "sparsity levels completed"))
}

#####################################################
# Figure 1 (a) - q=0.1, Brzyski et. al. (2015)
#####################################################

# plot FDR -------------------------
plot(n.relevant, FDR.max, ylim=c(0,0.15), type="b", lty=2,
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="gFDR for lambda_max")

# FDR nominal level
lines(n.relevant, fdr*(n.group-n.relevant)/n.group)

legend(90, 0.14, c("gFDR, q=0.1", "Theoretical upper bound"), lty=c(2,1), pch=c(1,NA))

# FDR error bars
FDR.max.se <- FDR.max.sd/sqrt(n.iter)
segments(n.relevant, FDR.max-2*FDR.max.se, n.relevant, FDR.max+2*FDR.max.se, col="blue")
segments(n.relevant-1, FDR.max-2*FDR.max.se, n.relevant+1, FDR.max-2*FDR.max.se, col="blue")
segments(n.relevant-1, FDR.max+2*FDR.max.se, n.relevant+1, FDR.max+2*FDR.max.se, col="blue")

#####################################################
# Figure 1 (b) - q=0.1, Brzyski et. al. (2015)
#####################################################

# plot FDR -------------------------
plot(n.relevant, FDR.mean, ylim=c(0,0.15), type="b", lty=2,
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="gFDR for lambda_mean")

# FDR nominal level
lines(n.relevant, fdr*(n.group-n.relevant)/n.group)

legend(90, 0.14, c("gFDR, q=0.1", "Theoretical upper bound"), lty=c(2,1), pch=c(1,NA))

# FDR error bars
FDR.mean.se <- FDR.mean.sd/sqrt(n.iter)
segments(n.relevant, FDR.mean-2*FDR.mean.se, n.relevant, FDR.mean+2*FDR.mean.se, col="blue")
segments(n.relevant-1, FDR.mean-2*FDR.mean.se, n.relevant+1, FDR.mean-2*FDR.mean.se, col="blue")
segments(n.relevant-1, FDR.mean+2*FDR.mean.se, n.relevant+1, FDR.mean+2*FDR.mean.se, col="blue")

####################################################################
# Figure 1 (c) - q=0.1, basic/relaxed lambda, Brzyski et. al. (2015)
####################################################################

# plot power -------------------------
plot(n.relevant, pow.max, ylim=c(0,1), type="b", col=1, pch=1,
     xlab="Number of relevant groups", ylab="Estimated power", 
     main="Power")
lines(n.relevant, pow.mean, type="b", col=2, pch=2)
legend(90, 0.4, c("Basic lambda, q=0.1", "Relaxed lambda, q=0.1"),
       lty=c(1,1), pch=c(1,2), col=c(1,2))
