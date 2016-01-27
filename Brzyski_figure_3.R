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
# Figure 3
####################################

fdr <- 0.1
n.iter <- 200

p <- 5000
n.group <- 1000
group <- rep(1:1000, each=5)
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
wt <- rep(sqrt(5), p)
B <- sqrt(4*log(n.group) * (1 - n.group^(-2/5)) - 5)
n.relevant <- floor(seq(1, 100, length=10))

# generate lambda
lambda.max <- lambdaGroupSLOPE(fdr=fdr, group=group,
                               wt=sqrt(group.length),
                               method="chiOrthoMax")
lambda.chi <- lambdaGroupSLOPE(fdr=fdr, n.obs=p, group=group, 
                               wt=sqrt(group.length),
                               method="chiEqual")

FDR.chi    <- rep(NA, length(n.relevant))
FDR.chi.sd <- rep(NA, length(n.relevant))
FDR.max    <- rep(NA, length(n.relevant))
FDR.max.sd <- rep(NA, length(n.relevant))

one.iteration <- function(n.signif){
  # model matrix
  X <- matrix(rnorm(p^2, sd=sqrt(1/p)), p, p)

  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) {
    X1 <- apply(X[ , group.id[[j]]], 1, sum)
    b[group.id[[j]]] <- (B / norm(as.matrix(X1), "f")) * rep(1, group.length[j])
  }

  # generate the response vector
  y <- X %*% b + rnorm(p, sd=1)

  # get Group SLOPE solution
  b.grpSLOPE.chi <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group,
                                                     wt=wt, lambda=lambda.chi,
                                                     verbose=FALSE)
  b.grpSLOPE.max <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group,
                                                     wt=wt, lambda=lambda.max,
                                                     verbose=FALSE)

  # FDR and power
  nonzero <- rep(NA, n.group)
  for (j in 1:n.group) { nonzero[j] <- (sum(b.grpSLOPE.chi$x[group.id[[j]]]^2) > 0) }
  truepos <- sum(nonzero[ind.relevant])
  falsepos <- sum(nonzero) - truepos
  FDR.chi <- falsepos / max(1, sum(nonzero))

  nonzero <- rep(NA, n.group)
  for (j in 1:n.group) { nonzero[j] <- (sum(b.grpSLOPE.max$x[group.id[[j]]]^2) > 0) }
  truepos <- sum(nonzero[ind.relevant])
  falsepos <- sum(nonzero) - truepos
  FDR.max <- falsepos / max(1, sum(nonzero))

  return(list(FDR.chi=FDR.chi, FDR.max=FDR.max))
}

for (k in 1:length(n.relevant)) {
  parallel.results <- foreach(i=1:n.iter) %dopar% {
    one.iteration(n.relevant[k])
  }

  FDR.chi.vec <- rep(NA, n.iter)
  FDR.max.vec <- rep(NA, n.iter)

  for (j in 1:n.iter) {
    FDR.chi.vec[j] <- parallel.results[[j]]$FDR.chi
    FDR.max.vec[j] <- parallel.results[[j]]$FDR.max
  }

  FDR.chi[k]    <- mean(FDR.chi.vec)
  FDR.chi.sd[k] <- sd(FDR.chi.vec)
  FDR.max[k]    <- mean(FDR.max.vec)
  FDR.max.sd[k] <- sd(FDR.max.vec)
  
  print(paste(k, "sparsity levels completed"))
}

#####################################################
# Figure 3 (a) - q=0.1, Brzyski et. al. (2015)
#####################################################

# plot FDR -------------------------
plot(n.relevant, FDR.max, ylim=c(0,0.3), type="b",
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="gFDR for lambda_max, q = 0.1")

# FDR nominal level
abline(fdr, 0)

# FDR error bars
FDR.max.se <- FDR.max.sd / sqrt(n.iter)
segments(n.relevant, FDR.max-FDR.max.se, n.relevant,
         FDR.max+FDR.max.se, col="blue")
segments(n.relevant-1, FDR.max-FDR.max.se, n.relevant+1,
         FDR.max-FDR.max.se, col="blue")
segments(n.relevant-1, FDR.max+FDR.max.se, n.relevant+1,
         FDR.max+FDR.max.se, col="blue")

#############################################################
# Figure 3 (b) - q=0.1, Brzyski et. al. (2015)
#############################################################

# plot FDR -------------------------
plot(n.relevant, FDR.chi, ylim=c(0,0.3), type="b",
     xlab="Number of relevant groups", ylab="Estimated gFDR", 
     main="gFDR for corrected lambda, q = 0.1")

# FDR nominal level
abline(fdr, 0)

# FDR error bars
FDR.chi.se <- FDR.chi.sd / sqrt(n.iter)
segments(n.relevant, FDR.chi-FDR.chi.se, n.relevant,
         FDR.chi+FDR.chi.se, col="blue")
segments(n.relevant-1, FDR.chi-FDR.chi.se, n.relevant+1,
         FDR.chi-FDR.chi.se, col="blue")
segments(n.relevant-1, FDR.chi+FDR.chi.se, n.relevant+1,
         FDR.chi+FDR.chi.se, col="blue")

#############################################################
# Figure 3 (c) - q=0.1, Brzyski et. al. (2015)
#############################################################

plot(lambda.max[1:100], type="l", lty=2, ylim=c(1.6,2.4),
     xlab="Index", ylab="Value of coefficient",
     main="Basic and corrected lambdas, q = 0.1")
lines(lambda.chi[1:100], lty=1)
legend(50, 2.3, c("Basic lambda", "Corrected lambda"), lty=c(2,1))
