#--------------------------------------------------------------------------
# This code is based on the simulation study of Figure 4 in D. Brzyski, 
# W. Su, M. Bogdan (2015), "Group SLOPE - adaptive selection of 
# groups of predictors" (http://arxiv.org/abs/1511.09078)
#--------------------------------------------------------------------------

library(grpSLOPE)
library(dplyr)
library(tidyr)
library(ggplot2)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# design matrix dimensions
n <- 1000
p <- 2000

#--- Set up global parameters for the simulation

set.seed(20161006)

# auxilliary function to get (group-wise) FDP and power from one grpSLOPE solution object
get_FDP_and_power <- function(result, true.relevant){
  truepos <- length(intersect(result$selected, true.relevant))
  falsepos <- length(result$selected) - truepos
  FDP <- falsepos / max(1, length(result$selected))
  pow <- truepos / length(true.relevant)
  return(c("FDP" = FDP, "pow" = pow))
}

# set the grouping structure
len <- 200
n.group <- 10 
group <- rep(1:n.group, each = len)
group.length <- rep(len, n.group)
group.id <- getGroupID(group)

# determine signal strength, such as used in Figure 1 in Brzyski et. al. (2015)
signal.strength <- sqrt(4*log(n.group) / (1 - n.group^(-2/len)) - len)

# considered numbers of truly relevant groups
sparsity <- seq(0, 0.9, length = 10)
n.relevant <- unique(floor(sparsity * n.group))

# how many times the simulation is repeated
n.replications <- 80 


#--- Simulation main loop

# run the simulations n.replications times at each level of n.relevant,
# and with each combination of lambda="max", lambda="mean", fdr=0.1, fdr=0.05
parallel.results <- vector(mode="list")
for (k in 1:length(n.relevant)) {
  cat(paste("sparsity level", k, "started"))

  parallel.results[[k]] <- foreach(i=1:n.replications, .combine = rbind) %dopar% {
    cat(".")

    # generate and normalize the model matrix
    X <- matrix(rnorm(n*p), n, p)
    X <- scale(X, center=TRUE, scale=FALSE)
    X <- apply(X, 2, function(x) x/sqrt(sum(x^2)) )

    # generate coeffient vector, pick relevant groups at random
    b <- rep(0, p)
    n.signif <- n.relevant[k]
    ind.relevant <- sample(1:n.group, n.signif)
    for (j in ind.relevant) {
      signals <- runif(group.length[j]) + 0.1
      X1 <- as.matrix(X[ , group.id[[j]]]) %*% signals
      b[group.id[[j]]] <- (signal.strength / norm(as.matrix(X1), "f")) * signals
    }

    # generate the response vector
    y <- X %*% b + rnorm(n)

    # get Group SLOPE solutions with different lambda and fdr values
    # (this has the same gFDR controlling properties with orthogonalize=FALSE too)
    lambda.1 <- grpSLOPE(X = X, y = y, group = group, fdr = 0.1, sigma = 1)
    lambda.05 <- grpSLOPE(X = X, y = y, group = group, fdr = 0.05, sigma = 1)

    # get the FDPs and powers of the grpSLOPE solutions
    true.relevant <- names(group.id)[ind.relevant]
    FDPs.and.powers <- rep(NA, 6)
    names(FDPs.and.powers) <- c("n.relevant", "replication",
                                "lambda.1.FDP", "lambda.1.power",
                                "lambda.05.FDP", "lambda.05.power")
    FDPs.and.powers["n.relevant"] <- length(true.relevant)
    FDPs.and.powers["replication"] <- i
    FDPs.and.powers[c("lambda.1.FDP", "lambda.1.power")] <- get_FDP_and_power(lambda.1, true.relevant)
    FDPs.and.powers[c("lambda.05.FDP", "lambda.05.power")] <- get_FDP_and_power(lambda.05, true.relevant)

    FDPs.and.powers
  }

  cat("done\n")
}


#--- Collect and summarize simulation results

# collect results in a data frame
results <- data.frame() %>% tbl_df
for(k in 1:length(n.relevant)) {
  results <- bind_rows(results, tbl_df(parallel.results[[k]]))
}

# get means and error bars for FDP and power
results.summary <- results %>% select(-replication) %>% group_by(n.relevant) %>% 
  summarize_all(funs(mean, lwr = (mean(.) - 2*sd(.)/sqrt(n.replications)),
                     upr = (mean(.) + 2*sd(.)/sqrt(n.replications))))


#--- Plot a figure analogous to Figure 4a in Brzyski et. al. (2015)

# pre-process for plotting
FDR.results <- results.summary %>% select(n.relevant, ends_with("FDP_mean")) %>% 
  gather(scenario, FDR, -n.relevant) %>%
  mutate(scenario = gsub("\\.FDP_mean", "", scenario))
lwr <- results.summary %>% select(n.relevant, ends_with("FDP_lwr")) %>% 
  gather(scenario, lwr, -n.relevant) %>%
  mutate(scenario = gsub("\\.FDP_lwr", "", scenario))
upr <- results.summary %>% select(n.relevant, ends_with("FDP_upr")) %>% 
  gather(scenario, upr, -n.relevant) %>%
  mutate(scenario = gsub("\\.FDP_upr", "", scenario))
FDR.results <- FDR.results %>% left_join(lwr) %>% left_join(upr)

# plot estimated FDR with error bars 
filename <- "./figures/gFDR_large_groups.pdf"
pdf(file = filename)

xend <- tail(n.relevant, 1)
ggplot(FDR.results) +
  geom_segment(mapping = aes(x = 0, xend = xend, y = 0.1, yend = 0.1), 
               linetype = 2, color = "black") +
  geom_segment(mapping = aes(x = 0, xend = xend, y = 0.05, yend = 0.05), 
               linetype = 2, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = FDR, color = scenario)) + 
  geom_point(mapping = aes(x = n.relevant, y = FDR, color = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, 
                            fill = scenario), alpha=0.25) +
  xlab("Number of significant groups") + ylab("Estimated gFDR +/- 2SE") +
  scale_color_discrete(name = "target gFDR", 
                       breaks = c("lambda.1", "lambda.05"),
                       labels=c("0.1", "0.05")) +
  scale_fill_discrete(name = "target gFDR", 
                      breaks = c("lambda.1", "lambda.05"),
                      labels=c("0.1", "0.05")) +
  coord_cartesian(ylim = c(0, 0.2)) +
  ggtitle(paste("n =", n, ", p=", p, "group size l =", len))

dev.off()

