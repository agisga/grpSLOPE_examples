library(microbenchmark)
library(grpreg)
# install the version of grpSLOPE that has an ADMM optimizer
devtools::install_github("agisga/grpSLOPE", ref = "optimization_and_ADMM")
library(grpSLOPE)

#--- function definitions

generate_data <- function(n = NA, p = NA, n_grp = NA) {
  A <- matrix(rnorm(n * p), n, p)
  # group membership per predictor
  grp <- sample(x = 1:n_grp, size = p, replace = TRUE)
  # weights: sqrt of the group size as weight for each group
  wt_per_grp <- sqrt(c(table(grp)))
  wt_per_predictor <- wt_per_grp[as.character(grp)]
  # a sparse vector of coefficients per predictor (i.e., the true solution)
  grp_signif <- 1:10
  ind_signif <- which(grp %in% grp_signif)
  x <- rep(0, p)
  x[ind_signif] <- rnorm(length(ind_signif), mean = 5, sd = 1)
  # response variable
  y <- A %*% x + rnorm(n)
  # normalize
  col_norms <- apply(A, 2, function(a) sqrt(sum((a - mean(a))^2)))
  A <- scale(A, center = TRUE, scale = col_norms)
  y <- scale(y, center = TRUE, scale = FALSE)
  # regularizing parameters
  lambda <- lambdaGroupSLOPE(method = "corrected", fdr = 0.1, n.obs = n,
                             group = grp, wt = wt_per_grp)
  return(list("y" = y, "A" = A, "grp" = grp, "lambda" = lambda,
              "wt_per_predictor" = wt_per_predictor))
}

admm_sim <- function(n = NULL, p = NULL, n_grp = NULL, rho = NULL,
                     init = "zeros", seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dat <- generate_data(n = n, p = p, n_grp = n_grp)
  if (init == "gLASSO") {
    cvgrplasso <- cv.grpreg(dat$A, dat$y, dat$grp,
                            penalty="grLasso", family="gaussian")
    x_init <- coef(cvgrplasso)[-1]
  } else {
    x_init <- NULL
  }
  result <- admmSolverGroupSLOPE(y = dat$y, A = dat$A, group = dat$grp,
                                 wt = dat$wt_per_predictor, lambda = dat$lambda,
                                 absolute.tol = 1e-6, relative.tol=1e-6,
                                 rho = rho, max.iter = 100000,
                                 z.init = x_init, verbose = FALSE)
  return(result$x)
}

fista_sim <- function(n = NULL, p = NULL, n_grp = NULL,
                      init = "zeros", seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dat <- generate_data(n = n, p = p, n_grp = n_grp)
  if (init == "gLASSO") {
    cvgrplasso <- cv.grpreg(dat$A, dat$y, dat$grp,
                            penalty="grLasso", family="gaussian")
    grpreg
    x_init <- coef(cvgrplasso)[-1]
  } else {
    x_init <- NULL
  }
  result <- proximalGradientSolverGroupSLOPE(y = dat$y, A = dat$A, group = dat$grp,
                                             wt = dat$wt_per_predictor,
                                             dual.gap.tol = 1e-6, infeas.tol = 1e-6,
                                             lambda = dat$lambda,
                                             x.init = x_init, verbose = FALSE)
  return(result$x)
}

#--- simulation

n <- 200
p <- 400
n_grp <- 100
set.seed(2018)

check_for_equal_coefs <- function(values, tol = 1e-3) {
  l <- length(values)
  max_error <- 0
  for (i in 1:(l-1)) {
    for (j in (i+1):l) {
      max_error_prev <- max_error
      max_error <- max(c(abs(values[[i]] - values[[j]])))
      max_error <- max(max_error_prev, max_error)
    }
  }
  max_error < tol
}

mbm <- microbenchmark(
  "FISTA: 0s init." = {
    b <- fista_sim(n = n, p = p, n_grp = n_grp)
  },
  "ADMM: rho=1e-2, 0s init." = {
    b <- admm_sim(rho = 1e-2, n = n, p = p, n_grp = n_grp)
  },
  "ADMM: rho=1e-1, 0s init." = {
    b <- admm_sim(rho = 1e-1, n = n, p = p, n_grp = n_grp)
  },
  "ADMM: rho=1, 0s init." = {
    b <- admm_sim(rho = 1, n = n, p = p, n_grp = n_grp)
  },
  "ADMM: rho=1e1, 0s init." = {
    b <- admm_sim(rho = 1e1, n = n, p = p, n_grp = n_grp)
  },
  "ADMM: rho=1e2, 0s init." = {
    b <- admm_sim(rho = 1e2, n = n, p = p, n_grp = n_grp)
  },
  "FISTA: gLASSO init." = {
    b <- fista_sim(n = n, p = p, n_grp = n_grp, init = "gLASSO")
  },
  "ADMM: rho=1e-2, gLASSO init." = {
    b <- admm_sim(rho = 1e-2, n = n, p = p, n_grp = n_grp, init = "gLASSO")
  },
  "ADMM: rho=1e-1, gLASSO init." = {
    b <- admm_sim(rho = 1e-1, n = n, p = p, n_grp = n_grp, init = "gLASSO")
  },
  "ADMM: rho=1, gLASSO init." = {
    b <- admm_sim(rho = 1, n = n, p = p, n_grp = n_grp, init = "gLASSO")
  },
  "ADMM: rho=1e1, gLASSO init." = {
    b <- admm_sim(rho = 1e1, n = n, p = p, n_grp = n_grp, init = "gLASSO")
  },
  "ADMM: rho=1e2, gLASSO init." = {
    b <- admm_sim(rho = 1e2, n = n, p = p, n_grp = n_grp, init = "gLASSO")
  },
  #check = check_for_equal_coefs,
  # (the check passes only if data are generated with same
  # seed within all compared functions)
  times = 100, unit = "s")

save(mbm, file = "mbm.RData")

#print(mbm)
#library(ggplot2)
#autoplot(mbm)
