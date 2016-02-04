library(grpSLOPE)
library(grpreg)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=10)

TestGroupSLOPE <- function(n.significant.blocks, n.subjects, signal, fdr, verbose=FALSE){

  # (1) Generate the nsubjects x 1050 design matrix X from the multivariate normal distribution 
  # with mean 0 and 1050x1050 covariance matrix Sigma,
  # where Sigma is block diagonal with blocks Gamma of sizes 
  # 5x5, 10x10, and 20x20 (30 blocks of each size),
  # such that Gamma[i,i] = 1 if i=j and =0.7 otherwise.
  nsubjects <- n.subjects
  blockcovmat5 <- diag(1,5,5)
  blockcovmat5[upper.tri(blockcovmat5,diag=F)] <- blockcovmat5[lower.tri(blockcovmat5, diag=F)] <- 0.7
  blockcovmat10 <- diag(1,10,10)
  blockcovmat10[upper.tri(blockcovmat10,diag=F)] <- blockcovmat10[lower.tri(blockcovmat10, diag=F)] <- 0.7
  blockcovmat20 <- diag(1,20,20)
  blockcovmat20[upper.tri(blockcovmat20,diag=F)] <- blockcovmat20[lower.tri(blockcovmat20, diag=F)] <- 0.7
  Sigma <- bdiag(rep(list(blockcovmat5, blockcovmat10, blockcovmat20), 30))
  withinblock.ind <- as.matrix(Sigma) # for (1.2) below
  withinblock.ind[upper.tri(withinblock.ind,diag=T)] <- 0 
  #increase between group cor
  Sigma <- as.matrix(Sigma)
  Sigma[which(Sigma==0)] <- 0.1
  #---
  Sigma.chol <- as.matrix(chol(Sigma))
  A <- matrix(rnorm(1050*nsubjects),nsubjects, 1050) %*% Sigma.chol

  # Normalize A
  A <- scale(A, center=TRUE, scale=FALSE)
  A <- apply(A, 2, function(x) x/sqrt(sum(x^2)) )

  # (1.2) Compute the average within and between blocks correlations
  corA <- abs(cor(A))
  withinblock.cor <- list("mean"=mean(c(corA[which(withinblock.ind!=0)])),
                          "var"=var(c(corA[which(withinblock.ind!=0)])))
  betweenblock.cor <- list("mean"=mean(c(corA[which(withinblock.ind==0)])),
                           "var"=var(c(corA[which(withinblock.ind==0)])))

  #(2) Create the vector b consisting of 90 blocks of sizes 5, 10, and 20 
  b <- rep(0, 1050)
  # indices of the first entries of each block
  block.start <- rep(NA, 90)
  block.count <- 1
  for (i in 1:30){
    block.start[3*i-2] <- block.count
    block.start[3*i-1] <- block.count + 5
    block.start[3*i] <- block.count + 15 
    block.count <- block.count + 35
  }
  # choose significant blocks
  nsignificant <- n.significant.blocks #number of significant blocks
  block.significant <- sort(sample(90,nsignificant))
  for(i in block.significant){
    if (i == 90){
      b[block.start[90]:1050] <- signal
    }
    else{
      b[block.start[i]:(block.start[i+1]-1)] <- (-1)^i * signal 
    }
  }

  #(3) Create the response vector y
  errorvector <- rnorm(nsubjects, 0, 1)
  y <- A %*% b + errorvector

  # (4) Compute the solution
  group <- vector()
  for (i in 1:30) {
    tmp <- rep((i-1)*3+c(1,2,3), c(5,10,20))
    group <- c(group, tmp)
  }
  m <- length(getGroupID(group))

  # (4.1) gaussianMC Group SLOPE
  grpslope.g <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="gaussianMC",
                         sigma=1, n.MC=20, MC.reps=5000, verbose=verbose,
                         orthogonalize=FALSE, normalize=FALSE)
  # (4.2) chiMC Group SLOPE
  grpslope.c <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="chiMC",
                         sigma=1, n.MC=20, MC.reps=5000, verbose=verbose,
                         normalize=FALSE)
  # (4.3) chiMean Group SLOPE
  grpslope.m <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="chiMean",
                         sigma=1, verbose=verbose, normalize=FALSE)

  #(5) Compute the number of true and the false discoveries
  #(5.1) gaussianMC Group SLOPE
  total.discoveries.g <- length(grpslope.g$selected) 
  true.discoveries.g <- length(intersect(grpslope.g$selected, block.significant))
  false.discoveries.g <- total.discoveries.g - true.discoveries.g 

  #(5.2) chiMC Group SLOPE
  total.discoveries.c <- length(grpslope.c$selected) 
  true.discoveries.c <- length(intersect(grpslope.c$selected, block.significant))
  false.discoveries.c <- total.discoveries.c - true.discoveries.c 

  #(5.3) chiMean Group SLOPE
  total.discoveries.m <- length(grpslope.m$selected) 
  true.discoveries.m <- length(intersect(grpslope.m$selected, block.significant))
  false.discoveries.m <- total.discoveries.m - true.discoveries.m 

  # (6) Print the output to screen
  if(verbose){
    print(paste("Number significant blocks: ", length(block.significant)))
    print(paste("Number discoveries: ", total.discoveries.g, "(gaussianMC)", 
                total.discoveries.c, "(chiMC)", total.discoveries.m, "(chiMean)"))
    print(paste("True discoveries: ", true.discoveries.g, "(gaussianMC)",
                true.discoveries.c, "(chiMC)", true.discoveries.m, "(chiMean)"))
    print(paste("False discoveries: ", false.discoveries.g, "(gaussianMC)",
                false.discoveries.c, "(chiMC)", false.discoveries.m, "(chiMean)"))
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
              nsignificant=nsignificant, 
              withinblock.cor=withinblock.cor,
              betweenblock.cor=betweenblock.cor,
              lambda.g=grpslope.g$lambda,
              lambda.c=grpslope.c$lambda,
              lambda.m=grpslope.m$lambda))
}


fdr <- 0.1
numblocks <- 90
n.replications <- 100
n.subjects <- 1000
sparsity.vec <- seq(0.05,0.8,l=10)
total.discoveries.g <- total.discoveries.c <- total.discoveries.m <- vector(mode="list")
true.discoveries.g <- true.discoveries.c <- true.discoveries.m <- vector(mode="list")
false.discoveries.g <- false.discoveries.c <- false.discoveries.m <- vector(mode="list")
lambda.g <- lambda.c <- lambda.m <- vector(mode="list")
withinblock.cor <- vector(mode="list")
betweenblock.cor <- vector(mode="list")

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
  withinblock.cor[[k]] <- vector(mode="list")
  betweenblock.cor[[k]] <- vector(mode="list")
}

parallel.results <- vector(mode="list")
for(k in 1:length(sparsity.vec)){
  parallel.results[[k]] <- foreach(i=1:n.replications) %dopar% {
    TestGroupSLOPE(n.significant.blocks=as.integer(sparsity.vec[k]*numblocks), 
                   n.subjects=n.subjects, signal=1, fdr=fdr, verbose=TRUE)
  }
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
    withinblock.cor[[k]][[i]] <- tmp$withinblock.cor
    betweenblock.cor[[k]][[i]] <- tmp$betweenblock.cor
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
  powervec <- true.discoveries.g[[j]]/as.integer(sparsity.vec[j]*numblocks)
  g.power[j] <- mean(powervec)
  g.power.sd[j] <- sd(powervec)

  fdrvec <- apply(cbind(false.discoveries.c[[j]], total.discoveries.c[[j]]), 1, fdrfun)
  c.fdr[j] <- mean(fdrvec)
  c.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.c[[j]]/as.integer(sparsity.vec[j]*numblocks)
  c.power[j] <- mean(powervec)
  c.power.sd[j] <- sd(powervec)

  fdrvec <- apply(cbind(false.discoveries.m[[j]], total.discoveries.m[[j]]), 1, fdrfun)
  m.fdr[j] <- mean(fdrvec)
  m.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.m[[j]]/as.integer(sparsity.vec[j]*numblocks)
  m.power[j] <- mean(powervec)
  m.power.sd[j] <- sd(powervec)
}

#*** Save everything ***

save(list = ls(all.names = TRUE), file = "test_lambdaMC_methods.RData")

#*** Do Plots ***

num.signif.blocks <- as.integer(sparsity.vec*numblocks)

# gaussianMC Group SLOPE

plot(num.signif.blocks, g.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="Group SLOPE", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, g.fdr, type="b", pch=1, lty=1, col="red")
legend(10, 0.4, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,100),c(0.1,0.1),lty=3, col="grey")

# plot fdr error bars
segments(num.signif.blocks, g.fdr-g.fdr.sd, num.signif.blocks, g.fdr+g.fdr.sd, col="red")
segments(num.signif.blocks-1, g.fdr-g.fdr.sd, num.signif.blocks+1, g.fdr-g.fdr.sd, col="red")
segments(num.signif.blocks-1, g.fdr+g.fdr.sd, num.signif.blocks+1, g.fdr+g.fdr.sd, col="red")

# plot power error bars
segments(num.signif.blocks, g.power-g.power.sd, num.signif.blocks, g.power+g.power.sd, col="blue")
segments(num.signif.blocks-1, g.power-g.power.sd, num.signif.blocks+1, g.power-g.power.sd, col="blue")
segments(num.signif.blocks-1, g.power+g.power.sd, num.signif.blocks+1, g.power+g.power.sd, col="blue")

# chiMC Group SLOPE

plot(num.signif.blocks, c.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="Group SLOPE", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, c.fdr, type="b", pch=1, lty=1, col="red")
legend(10, 0.4, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,100),c(0.1,0.1),lty=3, col="grey")

# plot fdr error bars
segments(num.signif.blocks, c.fdr-c.fdr.sd, num.signif.blocks, c.fdr+c.fdr.sd, col="red")
segments(num.signif.blocks-1, c.fdr-c.fdr.sd, num.signif.blocks+1, c.fdr-c.fdr.sd, col="red")
segments(num.signif.blocks-1, c.fdr+c.fdr.sd, num.signif.blocks+1, c.fdr+c.fdr.sd, col="red")

# plot power error bars
segments(num.signif.blocks, c.power-c.power.sd, num.signif.blocks, c.power+c.power.sd, col="blue")
segments(num.signif.blocks-1, c.power-c.power.sd, num.signif.blocks+1, c.power-c.power.sd, col="blue")
segments(num.signif.blocks-1, c.power+c.power.sd, num.signif.blocks+1, c.power+c.power.sd, col="blue")

# chiMean Group SLOPE

plot(num.signif.blocks, m.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="Group SLOPE", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, m.fdr, type="b", pch=1, lty=1, col="red")
legend(10, 0.4, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,100),c(0.1,0.1),lty=3, col="grey")

# plot fdr error bars
segments(num.signif.blocks, m.fdr-m.fdr.sd, num.signif.blocks, m.fdr+m.fdr.sd, col="red")
segments(num.signif.blocks-1, m.fdr-m.fdr.sd, num.signif.blocks+1, m.fdr-m.fdr.sd, col="red")
segments(num.signif.blocks-1, m.fdr+m.fdr.sd, num.signif.blocks+1, m.fdr+m.fdr.sd, col="red")

# plot power error bars
segments(num.signif.blocks, m.power-m.power.sd, num.signif.blocks, m.power+m.power.sd, col="blue")
segments(num.signif.blocks-1, m.power-m.power.sd, num.signif.blocks+1, m.power-m.power.sd, col="blue")
segments(num.signif.blocks-1, m.power+m.power.sd, num.signif.blocks+1, m.power+m.power.sd, col="blue")

# Selected Groups

avg.total.discoveries.g <- apply(as.matrix(as.data.frame(total.discoveries.g)), 2, mean)
avg.total.discoveries.c <- apply(as.matrix(as.data.frame(total.discoveries.c)), 2, mean)
avg.total.discoveries.m <- apply(as.matrix(as.data.frame(total.discoveries.m)), 2, mean)

plot(as.integer(sparsity.vec*numblocks), avg.total.discoveries.g, type="b",
     xlab="Number of significant groups", ylab="Groups selected",
     main="Number of selected groups", pch = 1, lty = 1, col="red")
lines(as.integer(sparsity.vec*numblocks), avg.total.discoveries.c, 
      type="b", pch = 2, lty=2, col="blue")
lines(as.integer(sparsity.vec*numblocks), avg.total.discoveries.m, 
      type="b", pch = 3, lty=3, col="green")
legend(30, 20, c("gaussianMC", "chiMC", "chiMean"), pch=1:3, lty=1:3, col=c("red", "blue", "green"))
