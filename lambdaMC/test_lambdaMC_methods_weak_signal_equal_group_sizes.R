library(grpSLOPE)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=10)

# Orthogonalize each group of columns of a matrix A.
# For i = 1, ..., m let A_i = A[ , group_i] and compute
# A_i[ , P] = Q %*% R, where P is a permutation vector.
#
orthogonalizeGroups <- function(X, group.id) {
  n.group <- length(group.id)

  getGroupQR <- function(ids) {
    submat <- X[ , ids]

    if (length(ids) == 1) {
      return(list(Q=as.matrix(submat), R=1, P=1))
    } else {
      submat.qr <- qr(submat, LAPACK=TRUE)
      return(list(Q=qr.Q(submat.qr),
                  R=qr.R(submat.qr),
                  P=submat.qr$pivot))
    }
  }

  return(lapply(group.id, getGroupQR))
}

TestGroupSLOPE <- function(n.significant, n.subjects, n.predictors, fdr, verbose=FALSE){

  # *All groups are assumed to have size m.*
  m <- 5

  n.group <- n.predictors / m
  # Generate the n.subjects x n.predictors design matrix A
  # from the multivariate normal distribution 
  # with mean 0 and covariance matrix Sigma.
  Sigma <- matrix(0.3, n.group, n.group)
  diag(Sigma) <- rep(0.7, n.group)
  Sigma <- Sigma %x% matrix(1, m, m) 
  diag(Sigma) <- 1

  Sigma.chol <- as.matrix(chol(Sigma))
  A <- matrix(rnorm(n.predictors * n.subjects), n.subjects, n.predictors) %*% Sigma.chol

  # Normalize A
  A <- scale(A, center=TRUE, scale=FALSE)
  A <- apply(A, 2, function(x) x/sqrt(sum(x^2)) )

  # Generate a group structure with groups of length m
  group <- rep(1:n.group, each=m)
  group.id <- getGroupID(group)
  group.length <- sapply(group.id, FUN=length)

  # Create the coefficient vector b
  b <- rep(0, n.predictors)
  block.significant <- sort(sample(n.group, n.significant))

  # set signal strength according to Brzyski et. al. (2015)
  signal.strength <- sqrt( 4*log(n.group) / (1 - n.group^(-2/m)) - m )

  for(i in block.significant){
    b[group.id[[i]]] <- signal.strength
  }

  # Create the response vector y

  ortho <- orthogonalizeGroups(A, group.id)
  A.ortho <- matrix(NA, n.subjects, n.predictors)
  for (i in 1:n.group) {
    A.ortho[ , group.id[[i]]] <- ortho[[i]]$Q
  }

  errorvector <- rnorm(n.subjects, 0, 1)
  y <- A.ortho %*% b + errorvector

  # Compute the solution

  # gaussianMC Group SLOPE
  grpslope.g <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="gaussianMC",
                         sigma=1, n.MC=10, MC.reps=5000,
                         orthogonalize=FALSE, normalize=FALSE)
  # chiMC Group SLOPE
  grpslope.c <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="chiMC",
                         sigma=1, n.MC=10, MC.reps=5000,
                         normalize=FALSE)
  # chiMean Group SLOPE
  grpslope.m <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="chiMean",
                         sigma=1, normalize=FALSE)

  # Compute the number of true and the false discoveries
  # gaussianMC Group SLOPE
  total.discoveries.g <- length(grpslope.g$selected) 
  true.discoveries.g <- length(intersect(grpslope.g$selected, block.significant))
  false.discoveries.g <- total.discoveries.g - true.discoveries.g 

  # chiMC Group SLOPE
  total.discoveries.c <- length(grpslope.c$selected) 
  true.discoveries.c <- length(intersect(grpslope.c$selected, block.significant))
  false.discoveries.c <- total.discoveries.c - true.discoveries.c 

  # chiMean Group SLOPE
  total.discoveries.m <- length(grpslope.m$selected) 
  true.discoveries.m <- length(intersect(grpslope.m$selected, block.significant))
  false.discoveries.m <- total.discoveries.m - true.discoveries.m 

  # Print output to screen
  if(verbose){
    print(paste("Number significant blocks: ", length(block.significant)))
    print(paste("Convergence status optimal: ", grpslope.g$optimal, "(gaussianMC),",
                grpslope.c$optimal, "(chiMC),", grpslope.m$optimal, "(chiMean)"))
    print(paste("Number discoveries: ", total.discoveries.g, "(gaussianMC),", 
                total.discoveries.c, "(chiMC),", total.discoveries.m, "(chiMean)"))
    print(paste("True discoveries: ", true.discoveries.g, "(gaussianMC),",
                true.discoveries.c, "(chiMC),", true.discoveries.m, "(chiMean)"))
    print(paste("False discoveries: ", false.discoveries.g, "(gaussianMC),",
                false.discoveries.c, "(chiMC),", false.discoveries.m, "(chiMean)"))
  }

  return(list(total.discoveries.g=total.discoveries.g, 
              true.discoveries.g=true.discoveries.g, 
              false.discoveries.g=false.discoveries.g, 
              total.discoveries.c=total.discoveries.c, 
              true.discoveries.c=true.discoveries.c, 
              false.discoveries.c=false.discoveries.c, 
              total.discoveries.m=total.discoveries.m, 
              true.discoveries.m=true.discoveries.m, 
              false.discoveries.m=false.discoveries.m, 
              n.significant=n.significant, 
              lambda.g=grpslope.g$lambda,
              lambda.c=grpslope.c$lambda,
              lambda.m=grpslope.m$lambda))
}


# CAREFUL: m is also hardcoded into TestGroupSLOPE()
m <- 5

fdr <- 0.1
n.predictors <- 2000
n.group <- n.predictors / m
n.replications <- 50 
n.subjects <- 2000
sparsity.vec <- seq(0,1,l=12)
total.discoveries.g <- total.discoveries.c <- total.discoveries.m <- vector(mode="list")
true.discoveries.g <- true.discoveries.c <- true.discoveries.m <- vector(mode="list")
false.discoveries.g <- false.discoveries.c <- false.discoveries.m <- vector(mode="list")
lambda.g <- lambda.c <- lambda.m <- vector(mode="list")

for(k in 1:length(sparsity.vec)){
  total.discoveries.g[[k]] <- rep(NA, n.replications)
  true.discoveries.g[[k]] <- rep(NA, n.replications)
  false.discoveries.g[[k]] <- rep(NA, n.replications)
  total.discoveries.c[[k]] <- rep(NA, n.replications)
  true.discoveries.c[[k]] <- rep(NA, n.replications)
  false.discoveries.c[[k]] <- rep(NA, n.replications)
  total.discoveries.m[[k]] <- rep(NA, n.replications)
  true.discoveries.m[[k]] <- rep(NA, n.replications)
  false.discoveries.m[[k]] <- rep(NA, n.replications)
  lambda.g[[k]] <- lambda.c[[k]] <- lambda.m[[k]] <- list()
}

parallel.results <- vector(mode="list")
for(k in 1:length(sparsity.vec)){
  parallel.results[[k]] <- foreach(i=1:n.replications) %dopar% {
    TestGroupSLOPE(n.significant=as.integer(sparsity.vec[k]*n.group), 
                   n.subjects=n.subjects, n.predictors=n.predictors,
                   fdr=fdr, verbose=FALSE)
  }
  print(paste("sparsity level", k, "done"))
}

for(k in 1:length(sparsity.vec)){
  for(i in 1:n.replications){
    tmp <- parallel.results[[k]][[i]]
    total.discoveries.g[[k]][i] <- tmp$total.discoveries.g
    true.discoveries.g[[k]][i] <- tmp$true.discoveries.g
    false.discoveries.g[[k]][i] <- tmp$false.discoveries.g
    total.discoveries.c[[k]][i] <- tmp$total.discoveries.c
    true.discoveries.c[[k]][i] <- tmp$true.discoveries.c
    false.discoveries.c[[k]][i] <- tmp$false.discoveries.c
    total.discoveries.m[[k]][i] <- tmp$total.discoveries.m
    true.discoveries.m[[k]][i] <- tmp$true.discoveries.m
    false.discoveries.m[[k]][i] <- tmp$false.discoveries.m
    lambda.g[[k]][[i]] <- tmp$lambda.g
    lambda.c[[k]][[i]] <- tmp$lambda.c
    lambda.m[[k]][[i]] <- tmp$lambda.m
  }
}

# *** Compute Power and FDR ***

fdrfun <- function(x) { 
  if(x[2] > 0){
    proportion = x[1] / x[2]
  } else{
    proportion = 0.0
  }
  return(proportion)
}

g.fdr <- c.fdr <- m.fdr <- rep(NA, length(sparsity.vec))
g.fdr.sd <- c.fdr.sd <- m.fdr.sd <- rep(NA, length(sparsity.vec))
g.power <- c.power <- m.power <- rep(NA, length(sparsity.vec))
g.power.sd <- c.power.sd <- m.power.sd <- rep(NA, length(sparsity.vec))

for(j in 1:length(sparsity.vec)){
  fdrvec <- apply(cbind(false.discoveries.g[[j]], total.discoveries.g[[j]]), 1, fdrfun)
  g.fdr[j] <- mean(fdrvec)
  g.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.g[[j]]/as.integer(sparsity.vec[j]*n.group)
  g.power[j] <- mean(powervec)
  g.power.sd[j] <- sd(powervec)

  fdrvec <- apply(cbind(false.discoveries.c[[j]], total.discoveries.c[[j]]), 1, fdrfun)
  c.fdr[j] <- mean(fdrvec)
  c.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.c[[j]]/as.integer(sparsity.vec[j]*n.group)
  c.power[j] <- mean(powervec)
  c.power.sd[j] <- sd(powervec)

  fdrvec <- apply(cbind(false.discoveries.m[[j]], total.discoveries.m[[j]]), 1, fdrfun)
  m.fdr[j] <- mean(fdrvec)
  m.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.m[[j]]/as.integer(sparsity.vec[j]*n.group)
  m.power[j] <- mean(powervec)
  m.power.sd[j] <- sd(powervec)
}

#*** Save everything ***

save(list = ls(all.names = TRUE), file = "test_lambdaMC_methods_weak_signal_equal_group_sizes_30x70.RData")

#*** Do Plots ***

num.signif.blocks <- as.integer(sparsity.vec*n.group)

postscript(file="weak_signal_equal_group_sizes_30x70.eps", width=1000, height=400)
par(mfrow=c(1,3))

# chiMean Group SLOPE

plot(num.signif.blocks, m.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="lambda = chiMean", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, m.fdr, type="b", pch=1, lty=1, col="red")
legend(20, 0.4, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,1500),c(0.1,0.1),lty=3, col="black")

# plot fdr error bars
m.fdr.se <- m.fdr.sd/sqrt(n.replications)
segments(num.signif.blocks, m.fdr-2*m.fdr.se, num.signif.blocks, m.fdr+2*m.fdr.se, col="red")
segments(num.signif.blocks-1, m.fdr-2*m.fdr.se, num.signif.blocks+1, m.fdr-2*m.fdr.se, col="red")
segments(num.signif.blocks-1, m.fdr+2*m.fdr.se, num.signif.blocks+1, m.fdr+2*m.fdr.se, col="red")

# plot power error bars
m.power.se <- m.power.sd/sqrt(n.replications)
segments(num.signif.blocks, m.power-2*m.power.se, num.signif.blocks, m.power+2*m.power.se, col="blue")
segments(num.signif.blocks-1, m.power-2*m.power.se, num.signif.blocks+1, m.power-2*m.power.se, col="blue")
segments(num.signif.blocks-1, m.power+2*m.power.se, num.signif.blocks+1, m.power+2*m.power.se, col="blue")

# chiMC Group SLOPE

plot(num.signif.blocks, c.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="lambda = chiMC", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, c.fdr, type="b", pch=1, lty=1, col="red")
legend(20, 0.4, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,1500),c(0.1,0.1),lty=3, col="black")

# plot fdr error bars
c.fdr.se <- c.fdr.sd/sqrt(n.replications)
segments(num.signif.blocks, c.fdr-2*c.fdr.se, num.signif.blocks, c.fdr+2*c.fdr.se, col="red")
segments(num.signif.blocks-1, c.fdr-2*c.fdr.se, num.signif.blocks+1, c.fdr-2*c.fdr.se, col="red")
segments(num.signif.blocks-1, c.fdr+2*c.fdr.se, num.signif.blocks+1, c.fdr+2*c.fdr.se, col="red")

# plot power error bars
c.power.se <- c.power.sd/sqrt(n.replications)
segments(num.signif.blocks, c.power-2*c.power.se, num.signif.blocks, c.power+2*c.power.se, col="blue")
segments(num.signif.blocks-1, c.power-2*c.power.se, num.signif.blocks+1, c.power-2*c.power.se, col="blue")
segments(num.signif.blocks-1, c.power+2*c.power.se, num.signif.blocks+1, c.power+2*c.power.se, col="blue")

# gaussianMC Group SLOPE

plot(num.signif.blocks, g.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="lambda = gaussianMC", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, g.fdr, type="b", pch=1, lty=1, col="red")
legend(20, 0.9, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,1500),c(0.1,0.1),lty=3, col="black")

# plot fdr error bars
g.fdr.se <- g.fdr.sd/sqrt(n.replications)
segments(num.signif.blocks, g.fdr-2*g.fdr.se, num.signif.blocks, g.fdr+2*g.fdr.se, col="red")
segments(num.signif.blocks-1, g.fdr-2*g.fdr.se, num.signif.blocks+1, g.fdr-2*g.fdr.se, col="red")
segments(num.signif.blocks-1, g.fdr+2*g.fdr.se, num.signif.blocks+1, g.fdr+2*g.fdr.se, col="red")

# plot power error bars
g.power.se <- g.power.sd/sqrt(n.replications)
segments(num.signif.blocks, g.power-2*g.power.se, num.signif.blocks, g.power+2*g.power.se, col="blue")
segments(num.signif.blocks-1, g.power-2*g.power.se, num.signif.blocks+1, g.power-2*g.power.se, col="blue")
segments(num.signif.blocks-1, g.power+2*g.power.se, num.signif.blocks+1, g.power+2*g.power.se, col="blue")

####################
dev.off()

# Selected Groups

postscript(file="selected_by_lambdaMC_weak_signal_equal_group_sizes_30x70.eps", horizontal=FALSE, width=400, height=400)

avg.total.discoveries.g <- apply(as.matrix(as.data.frame(total.discoveries.g)), 2, mean)
avg.total.discoveries.c <- apply(as.matrix(as.data.frame(total.discoveries.c)), 2, mean)
avg.total.discoveries.m <- apply(as.matrix(as.data.frame(total.discoveries.m)), 2, mean)

plot(as.integer(sparsity.vec*n.group), avg.total.discoveries.g, type="b",
     xlab="Number of significant groups", ylab="Groups selected",
     main="Number of selected groups", pch = 1, lty = 1, col="red",
     ylim=c(0, max(c(avg.total.discoveries.g, avg.total.discoveries.m, avg.total.discoveries.c))))
lines(as.integer(sparsity.vec*n.group), avg.total.discoveries.c, 
      type="b", pch = 2, lty=2, col="blue")
lines(as.integer(sparsity.vec*n.group), avg.total.discoveries.m, 
      type="b", pch = 3, lty=3, col="green")
legend(30, 20, c("gaussianMC", "chiMC", "chiMean"), pch=1:3, lty=1:3, col=c("red", "blue", "green"))
abline(0, 1, col="grey")

dev.off()
