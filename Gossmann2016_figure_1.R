library(grpSLOPE)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))


one.iteration <- function(n.signif, target.fdr) {
  n.group <- 100

  # group lengths are random
  group.length <- rbinom(n.group, 1000, 0.008) + 1
  group <- c()
  for (i in 1:n.group) {
    group <- c(group, rep(i, group.length[i]))
  }
  group.id <- getGroupID(group)

  n <- sum(group.length) 
  p <- length(group)

  # signal strength
  Bfun <- function(l) {
    B <- 4*log(n.group) / (1 - n.group^(-2/l)) - l
    return(sqrt(B))
  }
  signal.strength <- sum(Bfun(group.length)) / n.group

  # generate an orthogonal model matrix
  tmpmat <- matrix(rnorm(n*p), n, p)
  tmpmat.qr <- qr(tmpmat)
  X <- qr.Q(tmpmat.qr)

  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) {
    signals <- runif(group.length[j])
    X1 <- as.matrix(X[ , group.id[[j]]]) %*% signals
    b[group.id[[j]]] <- (signal.strength / norm(as.matrix(X1), "f")) * signals
  }

  # generate the response vector
  y <- X %*% b + rnorm(n, sd=1)

  # get Group SLOPE solution
  b.grpSLOPE <- grpSLOPE(X=X, y=y, group=group, fdr=target.fdr, sigma=1,
                         lambda="chiOrthoMean", verbose=FALSE,
                         orthogonalize=FALSE, normalize=FALSE)

  # FDR and power
  n.selected <- length(b.grpSLOPE$selected)
  true.relevant <- names(group.id)[ind.relevant]
  truepos <- intersect(b.grpSLOPE$selected, true.relevant)
  n.truepos <- length(truepos)
  n.falsepos <- n.selected - n.truepos
  FDR <- n.falsepos / max(1, n.selected)
  if (length(true.relevant) > 0) {
    pow <- n.truepos / length(true.relevant)
  } else {
    pow <- 1.0
  }

  return(list(FDR=FDR, pow=pow))
}

############
# Figure 1
############

set.seed(1)

n.relevant <- c(0, floor(seq(1, 60, length=7)))
fdr <- c(0.1, 0.05)
n.iter <- 500

FDR    <- list("0.1"=rep(NA, length(n.relevant)), "0.05"=rep(NA, length(n.relevant)))
FDR.sd <- list("0.1"=rep(NA, length(n.relevant)), "0.05"=rep(NA, length(n.relevant)))
pow    <- list("0.1"=rep(NA, length(n.relevant)), "0.05"=rep(NA, length(n.relevant)))
pow.sd <- list("0.1"=rep(NA, length(n.relevant)), "0.05"=rep(NA, length(n.relevant)))

for (target in 1:2) {
  for (k in 1:length(n.relevant)) {
    parallel.results <- foreach(i=1:n.iter) %dopar% {
      one.iteration(n.relevant[k], fdr[target])
    }

    FDR.vec  <- rep(NA, n.iter)
    pow.vec  <- rep(NA, n.iter)

    for (j in 1:n.iter) {
      FDR.vec[j]  <- parallel.results[[j]]$FDR
      pow.vec[j]  <- parallel.results[[j]]$pow
    }

    FDR[[target]][k] <- mean(FDR.vec)
    FDR.sd[[target]][k] <- sd(FDR.vec)
    pow[[target]][k] <- mean(pow.vec)
    pow.sd[[target]][k] <- sd(pow.vec)
    
    print(paste(k, "sparsity levels completed"))
  }
}

#####################################################
# Save
#####################################################

save(list=ls(all.names=TRUE), file="Gossmann2016_figure_1.RData")

#####################################################
# Figure 1 (a) - q=0.1
#####################################################

postscript(file="Gossmann2016_figure_1a.eps", horizontal=FALSE, width=400, height=400)

# plot FDR -------------------------
plot(n.relevant, FDR[[1]], ylim=c(0,0.25), type="b", lty=2,
     xlab="Number of relevant groups", ylab="Estimated gFDR")
lines(n.relevant, FDR[[2]], type="b", lty=3, pch=2)
legend(9, 0.24, c("gFDR, q=0.1", "gFDR, q=0.05"), lty=c(2,3), pch=c(1,2))

# FDR nominal level
abline(fdr[1], 0)
abline(fdr[2], 0)

# FDR error bars
for (target in 1:2) {
  FDR.se <- FDR.sd[[target]] / sqrt(n.iter)
  segments(n.relevant, FDR[[target]]-2*FDR.se, n.relevant, FDR[[target]]+2*FDR.se, col="blue")
  segments(n.relevant-1, FDR[[target]]-2*FDR.se, n.relevant+1, FDR[[target]]-2*FDR.se, col="blue")
  segments(n.relevant-1, FDR[[target]]+2*FDR.se, n.relevant+1, FDR[[target]]+2*FDR.se, col="blue")
}

dev.off()

####################################################################
# Figure 1 (b) - q=0.1
####################################################################

postscript(file="Gossmann_figure_1b.eps", horizontal=FALSE, width=400, height=400)

# plot power -------------------------
plot(n.relevant, pow[[1]], ylim=c(0,1), type="b", lty=1, pch=1,
     xlab="Number of relevant groups", ylab="Estimated power")
lines(n.relevant, pow[[2]], type="b", pch=2, lty=2)
legend(30, 0.4, c("Power, q=0.1", "Power, q=0.05"), lty=c(1,2), pch=c(1,2))

# Power error bars
for (target in 1:2) {
  pow.se <- pow.sd[[target]] / sqrt(n.iter)
  segments(n.relevant, pow[[target]]-2*pow.se, n.relevant, pow[[target]]+2*pow.se, col="blue")
  segments(n.relevant-1, pow[[target]]-2*pow.se, n.relevant+1, pow[[target]]-2*pow.se, col="blue")
  segments(n.relevant-1, pow[[target]]+2*pow.se, n.relevant+1, pow[[target]]+2*pow.se, col="blue")
}

dev.off()
