
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
registerDoParallel(cores=10)

####################################
# Figure 1
####################################

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
wt <- rep(NA, p)
for (j in 1:n.group) {
  wt[group.id[[j]]] <- sqrt(group.length[j])
}
#wt <- sqrt(group.length)

Bfun <- function(l) {
  sqrt(4*log(n.group) * (1 - n.group^(-2/l)) - l)
}
a <- sum(Bfun(group.length)) / sum(sqrt(group.length))

n.relevant <- floor(seq(1, 250, length=11))

# generate lambda
lambda.max <- lambdaGroupSLOPE(fdr=fdr, group=group,
                               wt=sqrt(group.length),
                               method="chiOrthoMax")

FDR    <- rep(NA, length(n.relevant))
FDR.sd <- rep(NA, length(n.relevant))
pow    <- rep(NA, length(n.relevant))
pow.sd <- rep(NA, length(n.relevant))

one.iteration <- function(n.signif){
  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) { b[group.id[[j]]] <-a }

  # generate the response vector
  y <- X %*% b + rnorm(p, sd=1)

  # get Group SLOPE solution
  b.grpSLOPE <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group,
                                                 wt=wt, lambda=lambda.max,
                                                 verbose=FALSE)

  # FDR and power
  nonzero <- rep(NA, n.group)
  for (j in 1:n.group) { nonzero[j] <- (sum(b.grpSLOPE$x[group.id[[j]]]^2) > 0) }
  truepos <- sum(nonzero[ind.relevant])
  falsepos <- sum(nonzero) - truepos
  FDR <- falsepos / max(1, sum(nonzero))
  pow <- truepos / length(ind.relevant)

  return(list(FDR=FDR, pow=pow))
}

for (k in 1:length(n.relevant)) {
  parallel.results <- foreach(i=1:n.iter) %dopar% {
    one.iteration(n.relevant[k])
  }

  FDR.vec <- rep(NA, n.iter)
  pow.vec <- rep(NA, n.iter)

  for (j in 1:n.iter) {
    FDR.vec[j] <- parallel.results[[j]]$FDR
    pow.vec[j] <- parallel.results[[j]]$pow
  }

  FDR[k] <- mean(FDR.vec)
  FDR.sd[k] <- sd(FDR.vec)
  pow[k] <- mean(pow.vec)
  pow.sd[k] <- sd(pow.vec)
  
  print(paste(k, "sparsity levels completed"))
}

#####################################################
# Figure 1 (a) - q=0.1, Brzyski et. al. (2015)
#####################################################

# plot FDR -------------------------
plot(n.relevant, FDR, ylim=c(0,0.15), type="b",
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="gFDR for lambda_max")

# FDR nominal level
lines(n.relevant, fdr*(n.group-n.relevant)/n.group)

# FDR error bars
segments(n.relevant, FDR-FDR.sd, n.relevant, FDR+FDR.sd, col="grey")
segments(n.relevant-1, FDR-FDR.sd, n.relevant+1, FDR-FDR.sd, col="grey")
segments(n.relevant-1, FDR+FDR.sd, n.relevant+1, FDR+FDR.sd, col="grey")

#############################################################
# Figure 1 (c) - q=0.1, basic lambda, Brzyski et. al. (2015)
#############################################################

# plot power -------------------------
plot(n.relevant, pow, ylim=c(0,1), type="b",
     xlab="Number of relevant groups", ylab="Estimated power", 
     main="power for lambda_max")

# pow error bars
segments(n.relevant, pow-pow.sd, n.relevant, pow+pow.sd, col="grey")
segments(n.relevant-1, pow-pow.sd, n.relevant+1, pow-pow.sd, col="grey")
segments(n.relevant-1, pow+pow.sd, n.relevant+1, pow+pow.sd, col="grey")
