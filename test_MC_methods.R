library(grpSLOPE)

fdr <- 0.1
p <- 500
n.group <- 100
group <- rep(1:100, each=5)
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
B <- sqrt(4*log(n.group) / (1 - n.group^(-2/5)) - 5)
n.signif <- 20

#######################################################
# Gaussian design
#######################################################

lambda.MC <- lambda.mean <- lambda.max <- list()

for (i in 1:4) {
  n <- i * p

  # model matrix
  X <- matrix(rnorm(n*p, sd=sqrt(1/n)), n, p)

  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) {
    X1 <- apply(X[ , group.id[[j]]], 1, sum)
    b[group.id[[j]]] <- (B / norm(as.matrix(X1), "f")) * rep(1, group.length[j])
  }

  # generate the response vector
  y <- X %*% b + rnorm(n, sd=1)

  # compute the lambda sequences
  lambda.MC[[i]] <- lambdaGroupSLOPE(fdr=fdr, group=group, A=X, y=y, wt=sqrt(group.length),
                                     method="chiMC", n.MC=20, MC.reps=10000)
  lambda.mean[[i]] <- lambdaGroupSLOPE(fdr=fdr, group=group, A=X, y=y, wt=sqrt(group.length),
                                       method="chiMean")
  lambda.max[[i]] <- lambdaGroupSLOPE(fdr=fdr, group=group, A=X, y=y, wt=sqrt(group.length),
                                      method="chiOrthoMax")

  print(paste("iteration", i, "done"))
}

png("./img/test_MC_methods_gaussian_design.png", height=800, width=800)
par(mfrow=c(2,2))

for (i in 1:4) {
  plot(lambda.max[[i]][1:20], type="l", lty=1, ylab="Lambda", main=paste("Gaussian design", i*p, "by", p))
  lines(lambda.mean[[i]][1:20], col=2, lty=1)
  lines(lambda.MC[[i]][1:20], col=3, lty=1)
}

legend(10, 2.0, c("chiOrthoMax", "chiMean", "chiMC"), lty=c(1,1,1), col=1:3)
dev.off()

#######################################################
# Orthogonal design
#######################################################

lambda.MC <- lambda.mean <- lambda.max <- list()

for (i in 1:4) {
  n <- i * p

  # model matrix
  tmpmat <- matrix(rnorm(n*p, sd=sqrt(1/n)), n, p)
  X <- qr.Q(qr(tmpmat))

  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) {
    signals <- runif(group.length[j])
    b[group.id[[j]]] <- (B / norm(as.matrix(signals), "f")) * signals
  }

  # generate the response vector
  y <- X %*% b + rnorm(n, sd=1)

  # compute the lambda sequences
  lambda.MC[[i]] <- lambdaGroupSLOPE(fdr=fdr, group=group, A=X, y=y, wt=sqrt(group.length),
                                     method="chiMC", n.MC=20, MC.reps=10000)
  lambda.mean[[i]] <- lambdaGroupSLOPE(fdr=fdr, group=group, A=X, y=y, wt=sqrt(group.length),
                                       method="chiMean")
  lambda.max[[i]] <- lambdaGroupSLOPE(fdr=fdr, group=group, A=X, y=y, wt=sqrt(group.length),
                                      method="chiOrthoMax")

  print(paste("iteration", i, "done"))
}

png("./img/test_MC_methods_orthogonal_design.png", height=800, width=800)
par(mfrow=c(2,2))

for (i in 1:4) {
  plot(lambda.max[[i]][1:20], type="l", lty=1, ylab="Lambda", main=paste("Orthogonal design", i*p, "by", p))
  lines(lambda.mean[[i]][1:20], col=2, lty=1)
  lines(lambda.MC[[i]][1:20], col=3, lty=1)
}

legend(10, 2.0, c("chiOrthoMax", "chiMean", "chiMC"), lty=c(1,1,1), col=1:3)
dev.off()
