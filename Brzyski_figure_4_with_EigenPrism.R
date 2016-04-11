###############################################################################
###############################################################################
###############################################################################
# This code performs a simulation under the same settings as in Figure 4 of 
# D. Brzyski, W. Su, M. Bogdan (2015), "Group SLOPE - adaptive selection of 
# groups of predictors" (http://arxiv.org/abs/1511.09078), but we use the
# EigenPrism procedure in order to get a point estimate for the noise level 
# sigma here. For detail on EigenPrism see L. Janson, R. Foygel Barber, 
# and E. Candès (2015), "EigenPrism: Inference for High-Dimensional Signal-to-Noise 
# Ratios" (http://arxiv.org/abs/1505.02097).
###############################################################################
###############################################################################
###############################################################################

library(grpSLOPE)
library(CVXfromR)

set.seed(1)

fdr <- 0.1
n.iter <- 2#200

n.group <- 1000
group.length <- rbinom(n.group, 1000, 0.008)
group <- c()
for (i in 1:n.group) {
  group <- c(group, rep(i, group.length[i]))
}
group.id <- getGroupID(group)

n <- 5000
p <- length(group)

Bfun <- function(l) {
  B <- 4*log(n.group) / (1 - n.group^(-2/l)) - l
  return(sqrt(B))
}
signal.strength <- sum(Bfun(group.length)) / n.group

n.relevant <- floor(seq(1, 60, length=7))

FDR      <- rep(NA, length(n.relevant))
FDR.sd   <- rep(NA, length(n.relevant))
pow      <- rep(NA, length(n.relevant))
pow.sd   <- rep(NA, length(n.relevant))
sigma    <- rep(NA, length(n.relevant))
sigma.sd <- rep(NA, length(n.relevant))

# this list holds TRUE/FALSE values, which signify when the EigenPrism estimate
# of sigma is positive or negotive
is.sigma.positive <- list()

one.iteration <- function(n.signif){
  # generate and normalize the model matrix
  X <- matrix(rnorm(n*p), n, p)
  X <- scale(X, center=TRUE, scale=FALSE)
  X <- apply(X, 2, function(x) x/sqrt(sum(x^2)) )

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

  # SVD
  X.svd <- svd(X, nv=0)
  EigenPrism.lambda <- X.svd$d^2 / p

  # solve the EigenPrism optimization problem 
  setup.dir <- "/home/agossman/cvx"
  cvx.code <- paste("variable t",
                    "variable w(n)",
                    "minimize t",
                    "subject to",
                    "sum(w) == 1",
                    "sum(w . * lambda) == 0",
                    "norm([w; (t/2-1)/2])  <= (t/2+1)/2",
                    "norm([w . * lambda; (t/2-1)/2]) <= (t/2+1)/2",
                    sep=";")
  EigenPrism <- CallCVX(cvx.code, const.vars=list(n=n, lambda=EigenPrism.lambda),
                        opt.var.names=c("t", "w"), setup.dir=setup.dir,
                        cvx.modifiers="quiet")

  # get EigenPrism point estimate for sigma
  z <- t(X.svd$u) %*% y
  EigenPrism.sigma <- sqrt(t(z^2) %*% EigenPrism$w)

  # check if point estimate is positive;
  # if it is not, then use sample SD of y as estimate for sigma
  positive <- TRUE
  if (EigenPrism.sigma > 0) {
    sigma <- EigenPrism.sigma
  } else {
    sigma <- sd(y)
    positive <- FALSE
  }

  # get Group SLOPE solution
  b.grpSLOPE <- grpSLOPE(X=X, y=y, group=group, fdr=fdr,
                         lambda="chiMean", sigma=sigma, verbose=FALSE,
                         orthogonalize=FALSE, normalize=FALSE)

  # FDR and power
  n.selected <- length(b.grpSLOPE$selected)
  true.relevant <- names(group.id)[ind.relevant]
  truepos <- intersect(b.grpSLOPE$selected, true.relevant)
  n.truepos <- length(truepos)
  n.falsepos <- n.selected - n.truepos
  FDR <- n.falsepos / max(1, n.selected)
  pow <- n.truepos / length(true.relevant)

  return(list(FDR=FDR, pow=pow, sigma=sigma, positive=positive))
}

for (k in 1:length(n.relevant)) {
  parallel.results <- foreach(i=1:n.iter) %dopar% {
    one.iteration(n.relevant[k])
  }

  FDR.vec   <- rep(NA, n.iter)
  pow.vec   <- rep(NA, n.iter)
  sigma.vec <- rep(NA, n.iter)
  is.sigma.positive[[k]] <- rep(NA, n.iter)

  for (j in 1:n.iter) {
    FDR.vec[j]   <- parallel.results[[j]]$FDR
    pow.vec[j]   <- parallel.results[[j]]$pow
    sigma.vec[j] <- parallel.results[[j]]$sigma
    is.sigma.positive[[k]][j] <- parallel.results[[j]]$positive
  }

  FDR[k] <- mean(FDR.vec)
  FDR.sd[k] <- sd(FDR.vec)
  pow[k] <- mean(pow.vec)
  pow.sd[k] <- sd(pow.vec)
  sigma[k] <- mean(sigma.vec)
  sigma.sd[k] <- sd(sigma.vec)
  
  print(paste(k, "sparsity levels completed"))
}

#####################################################
# Save
#####################################################

save(list=ls(all.names=TRUE), file="Brzyski_figure_4_with_EigenPrism.RData")

#####################################################
# Analog to Figure 4 (a) - q=0.1, Brzyski et. al. (2015)
#####################################################

postscript(file="Brzyski_figure_4a_with_EigenPrism.eps", 
           horizontal=FALSE, width=400, height=400)

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

dev.off()

####################################################################
# Analog to Figure 4 (b) - q=0.1, Brzyski et. al. (2015)
####################################################################

postscript(file="Brzyski_figure_4b_with_EigenPrism.eps", 
           horizontal=FALSE, width=400, height=400)

# plot power -------------------------
plot(n.relevant, pow, ylim=c(0,1), type="b", col=1, pch=1,
     xlab="Number of relevant groups", ylab="Estimated power", 
     main="Power")

# Power error bars
pow.se <- pow.sd/sqrt(n.iter)
segments(n.relevant, pow-2*pow.se, n.relevant, pow+2*pow.se, col="blue")
segments(n.relevant-1, pow-2*pow.se, n.relevant+1, pow-2*pow.se, col="blue")
segments(n.relevant-1, pow+2*pow.se, n.relevant+1, pow+2*pow.se, col="blue")

dev.off()

####################################################################
# Analog to Figure 4 (c) - q=0.1, Brzyski et. al. (2015)
####################################################################

postscript(file="Brzyski_figure_4c_with_EigenPrism.eps", 
           horizontal=FALSE, width=400, height=400)

hist(group.length, xlab = "Group size",
     main = "Histogram of group sizes")

dev.off()
