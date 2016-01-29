###############################################################################
###############################################################################
###############################################################################
# Data generation, analysis and plots for figures 4 (a), (b) and 6 (a) in 
# A. Gossmann, S. Cao, Y.-P. Wang (2015), "Identification of Significant
# Genetic Variants via SLOPE, and Its Extension to Group SLOPE" 
# (http://dx.doi.org/10.1145/2808719.2808743).
#
# Data generation, analysis and plots for figures 5 (a), (b) and 6 (b) are
# obtained by the same code with a change of parameters
###############################################################################
###############################################################################
###############################################################################

library(grpSLOPE)
library(grpreg)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=10)

TestGroupSLOPE <- function(n.significant.blocks, n.subjects, signal, verbose=FALSE){

  # (1) Generate the nsubjects x 1050 design matrix X from the multivariate normal distribution 
  # with mean 0 and 1050x1050 covariance matrix Sigma,
  # where Sigma is block diagonal with blocks Gamma of sizes 
  # 5x5, 10x10, and 20x20 (30 blocks of each size),
  # such that Gamma[i,i] = 1 if i=j and =0.99 otherwise.
  nsubjects <- n.subjects
  blockcovmat5 <- diag(1,5,5)
  blockcovmat5[upper.tri(blockcovmat5,diag=F)] <- blockcovmat5[lower.tri(blockcovmat5, diag=F)] <- 0.99
  blockcovmat10 <- diag(1,10,10)
  blockcovmat10[upper.tri(blockcovmat10,diag=F)] <- blockcovmat10[lower.tri(blockcovmat10, diag=F)] <- 0.99
  blockcovmat20 <- diag(1,20,20)
  blockcovmat20[upper.tri(blockcovmat20,diag=F)] <- blockcovmat20[lower.tri(blockcovmat20, diag=F)] <- 0.99
  Sigma <- bdiag(rep(list(blockcovmat5, blockcovmat10, blockcovmat20), 30))
  withinblock.ind <- as.matrix(Sigma) # for (1.2) below
  withinblock.ind[upper.tri(withinblock.ind,diag=T)] <- 0 
  #increase between group cor
  #Sigma <- as.matrix(Sigma)
  #Sigma[which(Sigma==0)] <- 0.3
  #---
  Sigma.chol <- as.matrix(chol(Sigma))
  A <- matrix(rnorm(mean=0, sd=1/1050, 1050*nsubjects),nsubjects, 1050) %*% Sigma.chol
  # A needs colnorms 1
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
  errorvector <- rnorm(nsubjects,0,1)
  y <- A %*% b + errorvector

  # (4) Compute the solution
  # (4.1) Group SLOPE
  fdr <- 0.10
  group <- vector()
  for (i in 1:30) {
    tmp <- rep((i-1)*3+c(1,2,3), c(5,10,20))
    group <- c(group, tmp)
  }
  m <- length(getGroupID(group))

  grpslope <- grpSLOPE(X=A, y=y, group=group, fdr=fdr, lambda="gaussianMC",
                       sigma=1, n.MC=20, MC.reps=5000, verbose=verbose)
  # (4.2) CV Group LASSO
  cvgrplasso <- cv.grpreg(A, y, group, penalty="grLasso", family="gaussian")

  #(5) Compute the true and the false discovery rates
  #(5.1) Group SLOPE
  total.discoveries.slope <- 0
  true.discoveries.slope <- 0
  false.discoveries.slope <- 0
  for(i in 1:(m-1)){
    if(sqrt(sum(grpslope$beta[block.start[i]:(block.start[i+1]-1)]^2)) > 1e-2){
      total.discoveries.slope <- total.discoveries.slope + 1
      if(i %in% block.significant) true.discoveries.slope <- true.discoveries.slope + 1
      else false.discoveries.slope <- false.discoveries.slope + 1
    }
  }
  if(sqrt(sum(grpslope$beta[block.start[m]:1050]^2)) > 1e-2){
    total.discoveries.slope <- total.discoveries.slope + 1
    if(m %in% block.significant) true.discoveries.slope <- true.discoveries.slope + 1
    else false.discoveries.slope <- false.discoveries.slope + 1
  }
  #(5.2) CV Group LASSO
  cvlasso.coef <- coef(cvgrplasso)[-1]
  total.discoveries.cvlasso <- 0
  true.discoveries.cvlasso <- 0
  false.discoveries.cvlasso <- 0
  for(i in 1:(m-1)){
    if(sqrt(sum(cvlasso.coef[block.start[i]:(block.start[i+1]-1)]^2)) > 1e-2){
      total.discoveries.cvlasso <- total.discoveries.cvlasso + 1
      if(i %in% block.significant) true.discoveries.cvlasso <- true.discoveries.cvlasso + 1
      else false.discoveries.cvlasso <- false.discoveries.cvlasso + 1
    }
  }
  if(sqrt(sum(cvlasso.coef[block.start[m]:1050]^2)) > 1e-2){
    total.discoveries.cvlasso <- total.discoveries.cvlasso + 1
    if(m %in% block.significant) true.discoveries.cvlasso <- true.discoveries.cvlasso + 1
    else false.discoveries.cvlasso <- false.discoveries.cvlasso + 1
  }


  # (6) Print the output to screen
  if(verbose){
    print(paste("Number significant blocks: ", length(block.significant)))
    print(paste("Number discoveries: ", total.discoveries.slope, "(GrpSLOPE)", 
                total.discoveries.cvlasso, "(CV GrpLASSO)"))
    print(paste("True discoveries: ", true.discoveries.slope, "(GrpSLOPE)",
                true.discoveries.cvlasso, "(CV GrpLASSO)"))
    print(paste("False discoveries: ", false.discoveries.slope, "(GrpSLOPE)",
                false.discoveries.cvlasso, "(CV GrpLASSO)"))
    print(paste("Number of iterations and runtime: ", grpslope$iter, grpslope$runtime)) 
  }

  return(list(SLOPEsol=grpslope$x, 
              total.discoveries.slope=total.discoveries.slope, 
              true.discoveries.slope=true.discoveries.slope, 
              false.discoveries.slope=false.discoveries.slope, 
              total.discoveries.cvlasso=total.discoveries.cvlasso, 
              true.discoveries.cvlasso=true.discoveries.cvlasso, 
              false.discoveries.cvlasso=false.discoveries.cvlasso, 
              nsignificant=nsignificant, 
              iter=grpslope$iter,
              withinblock.cor=withinblock.cor, betweenblock.cor=betweenblock.cor))
}


numblocks <- 90
n.replications <- 1000
n.subjects <- 200
sparsity.vec <- seq(0.05,0.8,l=10)
total.discoveries.slope <- vector(mode="list")
true.discoveries.slope <- vector(mode="list")
false.discoveries.slope <- vector(mode="list")
total.discoveries.cvlasso <- vector(mode="list")
true.discoveries.cvlasso <- vector(mode="list")
false.discoveries.cvlasso <- vector(mode="list")
iter <- vector(mode="list")
withinblock.cor <- vector(mode="list")
betweenblock.cor <- vector(mode="list")

for(k in 1:length(sparsity.vec)){
  total.discoveries.slope[[k]] <- rep(NA, n.replications)
  true.discoveries.slope[[k]] <- rep(NA, n.replications)
  false.discoveries.slope[[k]] <- rep(NA, n.replications)
  total.discoveries.cvlasso[[k]] <- rep(NA, n.replications)
  true.discoveries.cvlasso[[k]] <- rep(NA, n.replications)
  false.discoveries.cvlasso[[k]] <- rep(NA, n.replications)
  iter[[k]] <- rep(NA, n.replications)
  withinblock.cor[[k]] <- vector(mode="list")
  betweenblock.cor[[k]] <- vector(mode="list")
}

parallel.results <- vector(mode="list")
for(k in 1:length(sparsity.vec)){
  parallel.results[[k]] <- foreach(i=1:n.replications) %dopar% {
    TestGroupSLOPE(n.significant.blocks=as.integer(sparsity.vec[k]*numblocks), 
                   n.subjects=n.subjects, signal=1, verbose=TRUE)
  }
}

for(k in 1:length(sparsity.vec)){
  for(i in 1:n.replications){
    tmp <- parallel.results[[k]][[i]]
    total.discoveries.slope[[k]][i] <- tmp$total.discoveries.slope
    true.discoveries.slope[[k]][i] <- tmp$true.discoveries.slope
    false.discoveries.slope[[k]][i] <- tmp$false.discoveries.slope
    total.discoveries.cvlasso[[k]][i] <- tmp$total.discoveries.cvlasso
    true.discoveries.cvlasso[[k]][i] <- tmp$true.discoveries.cvlasso
    false.discoveries.cvlasso[[k]][i] <- tmp$false.discoveries.cvlasso
    iter[[k]][i] <- tmp$iter
    withinblock.cor[[k]][[i]] <- tmp$withinblock.cor
    betweenblock.cor[[k]][[i]] <- tmp$betweenblock.cor
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

slope.fdr <- rep(NA, length(sparsity.vec))
slope.fdr.sd <- rep(NA, length(sparsity.vec))
slope.power <- rep(NA, length(sparsity.vec))
slope.power.sd <- rep(NA, length(sparsity.vec))
cvlasso.fdr <- rep(NA, length(sparsity.vec))
cvlasso.fdr.sd <- rep(NA, length(sparsity.vec))
cvlasso.power <- rep(NA, length(sparsity.vec))
cvlasso.power.sd <- rep(NA, length(sparsity.vec))

for(j in 1:length(sparsity.vec)){
  fdrvec <- apply(cbind(false.discoveries.slope[[j]], total.discoveries.slope[[j]]), 1, fdrfun)
  slope.fdr[j] <- mean(fdrvec)
  slope.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.slope[[j]]/as.integer(sparsity.vec[j]*numblocks)
  slope.power[j] <- mean(powervec)
  slope.power.sd[j] <- sd(powervec)
  
  fdrvec <- apply(cbind(false.discoveries.cvlasso[[j]], total.discoveries.cvlasso[[j]]), 1, fdrfun)
  cvlasso.fdr[j] <- mean(fdrvec)
  cvlasso.fdr.sd[j] <- sd(fdrvec)
  powervec <- true.discoveries.cvlasso[[j]]/as.integer(sparsity.vec[j]*numblocks)
  cvlasso.power[j] <- mean(powervec)
  cvlasso.power.sd[j] <- sd(powervec)
}

#*** Do Plots ***

num.signif.blocks <- as.integer(sparsity.vec*numblocks)

# SLOPE

plot(num.signif.blocks, slope.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="Group SLOPE", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, slope.fdr, type="b", pch=1, lty=1, col="red")
legend(10, 0.4, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))
lines(c(0,100),c(0.1,0.1),lty=3, col="grey")

# plot fdr error bars
segments(num.signif.blocks, slope.fdr-slope.fdr.sd, num.signif.blocks, slope.fdr+slope.fdr.sd, col="red")
segments(num.signif.blocks-1, slope.fdr-slope.fdr.sd, num.signif.blocks+1, slope.fdr-slope.fdr.sd, col="red")
segments(num.signif.blocks-1, slope.fdr+slope.fdr.sd, num.signif.blocks+1, slope.fdr+slope.fdr.sd, col="red")

# plot power error bars
segments(num.signif.blocks, slope.power-slope.power.sd, num.signif.blocks, slope.power+slope.power.sd, col="blue")
segments(num.signif.blocks-1, slope.power-slope.power.sd, num.signif.blocks+1, slope.power-slope.power.sd, col="blue")
segments(num.signif.blocks-1, slope.power+slope.power.sd, num.signif.blocks+1, slope.power+slope.power.sd, col="blue")

# CV LASSO

plot(num.signif.blocks, cvlasso.power, type="b", xlab="Number of significant groups", 
     ylab="FDR and Power", main="Group LASSO with cross-validation", ylim=c(0,1), pch=2, lty=2, col="blue")
lines(num.signif.blocks, cvlasso.fdr, type="b", pch=1, lty=1, col="red")
legend(10, 0.25, c("Power", "False discovery rate"), pch=2:1, lty=2:1, col=c("blue", "red"))

# plot fdr error bars
segments(num.signif.blocks, cvlasso.fdr-cvlasso.fdr.sd, num.signif.blocks, cvlasso.fdr+cvlasso.fdr.sd, col="red")
segments(num.signif.blocks-1, cvlasso.fdr-cvlasso.fdr.sd, num.signif.blocks+1, cvlasso.fdr-cvlasso.fdr.sd, col="red")
segments(num.signif.blocks-1, cvlasso.fdr+cvlasso.fdr.sd, num.signif.blocks+1, cvlasso.fdr+cvlasso.fdr.sd, col="red")

# plot power error bars
segments(num.signif.blocks, cvlasso.power-cvlasso.power.sd, num.signif.blocks, cvlasso.power+cvlasso.power.sd, col="blue")
segments(num.signif.blocks-1, cvlasso.power-cvlasso.power.sd, num.signif.blocks+1, cvlasso.power-cvlasso.power.sd, col="blue")
segments(num.signif.blocks-1, cvlasso.power+cvlasso.power.sd, num.signif.blocks+1, cvlasso.power+cvlasso.power.sd, col="blue")

# Selected Groups

avg.total.discoveries.slope <- apply(as.matrix(as.data.frame(total.discoveries.slope)), 2, mean)
avg.total.discoveries.cvlasso <- apply(as.matrix(as.data.frame(total.discoveries.cvlasso)), 2, mean)

plot(as.integer(sparsity.vec*numblocks), avg.total.discoveries.slope, type="b",
     xlab="Number of significant groups", ylab="Groups selected",
     main="Number of selected groups", pch = 1, lty = 1, col="red")
lines(as.integer(sparsity.vec*numblocks), avg.total.discoveries.cvlasso, 
      type="b", pch = 2, lty=2, col="blue")
legend(30, 20, c("Group SLOPE", "CV Group LASSO"), pch=1:2, lty=1:2, col=c("red", "blue"))
