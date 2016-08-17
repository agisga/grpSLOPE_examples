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


#--- Set up global parameters for the simulation

# auxilliary function to get (group-wise) FDP and power from one grpSLOPE solution object
get_FDP_and_power <- function(result, true.relevant){
  truepos <- length(intersect(result$selected, true.relevant))
  falsepos <- length(result$selected) - truepos
  FDP <- falsepos / max(1, length(result$selected))
  pow <- truepos / length(true.relevant)
  return(c("FDP" = FDP, "pow" = pow))
}

# set the grouping structure
n.group <- 1000
group.length <- rbinom(n.group, 1000, 0.008)
group <- c()
for (i in 1:n.group) {
  group <- c(group, rep(i, group.length[i]))
}
group.id <- getGroupID(group)

# design matrix dimensions
n <- 5000
p <- length(group)

# determine signal strength, such as used in Figure 1 in Brzyski et. al. (2015)
Bfun <- function(l) {
  sqrt(4*log(n.group) / (1 - n.group^(-2/l)) - l)
}
signal.strength <- sum(Bfun(group.length)) / n.group

# considered numbers of truly relevant groups
n.relevant <- floor(seq(1, 60, length=7))

# how many times the simulation is repeated
n.replications <- 200 


#--- Simulation main loop

# run the simulations n.replications times at each level of n.relevant,
# and with each combination of lambda="max", lambda="mean", fdr=0.1, fdr=0.05
set.seed(20160807)
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
      signals <- runif(group.length[j])
      X1 <- as.matrix(X[ , group.id[[j]]]) %*% signals
      b[group.id[[j]]] <- (signal.strength / norm(as.matrix(X1), "f")) * signals
    }

    # generate the response vector
    y <- X %*% b + rnorm(n)

    # get Group SLOPE solutions with different lambda and fdr values
    lambda.1 <- grpSLOPE(X = X, y = y, group = group, fdr = 0.1,
                         lambda = "corrected", verbose = FALSE,
                         orthogonalize = FALSE, normalize = FALSE)
    lambda.05 <- grpSLOPE(X = X, y = y, group = group, fdr = 0.05,
                          lambda = "corrected", verbose = FALSE,
                          orthogonalize = FALSE, normalize = FALSE)

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
save(parallel.results, file = "./RData/Brzyski_4_results.RData") # keep a backup


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
png(file = "./figures/Brzyski_4a.png", width = 600, height = 480)

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
  coord_cartesian(ylim = c(0, 0.2))

dev.off()
 

#--- Plot a figure analogous to Figure 1b in Brzyski et. al. (2015)

# pre-process for plotting
power.results <- results.summary %>% select(n.relevant, ends_with("power_mean")) %>% 
  gather(scenario, power, -n.relevant) %>%
  mutate(scenario = gsub("\\.power_mean", "", scenario))
lwr <- results.summary %>% select(n.relevant, ends_with("power_lwr")) %>% 
  gather(scenario, lwr, -n.relevant) %>%
  mutate(scenario = gsub("\\.power_lwr", "", scenario))
upr <- results.summary %>% select(n.relevant, ends_with("power_upr")) %>% 
  gather(scenario, upr, -n.relevant) %>%
  mutate(scenario = gsub("\\.power_upr", "", scenario))
power.results <- power.results %>% left_join(lwr) %>% left_join(upr)

# plot estimated power with error bars 
png(file = "./figures/Brzyski_4b.png", width = 600, height = 480)

ggplot(power.results) +
  geom_line(mapping = aes(x = n.relevant, y = power, color = scenario)) + 
  geom_point(mapping = aes(x = n.relevant, y = power, color = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, 
                            fill = scenario), alpha=0.25) +
  xlab("Number of significant groups") + ylab("Estimated power +/- 2SE") +
  scale_color_discrete(name = "target gFDR", 
                       breaks = c("lambda.1", "lambda.05"),
                       labels=c("0.1", "0.05")) +
  scale_fill_discrete(name = "target gFDR", 
                      breaks = c("lambda.1", "lambda.05"),
                      labels=c("0.1", "0.05")) +
  coord_cartesian(ylim = c(0, 0.8))

dev.off()

#--- Plot a figure analogous to Figure 1c in Brzyski et. al. (2015)

# plot histogram of group sizes 
png(file = "./figures/Brzyski_4c.png", width = 600, height = 480)

data.frame(len = group.length) %>%
  ggplot(aes(len)) + geom_histogram(binwidth = 1, color = "grey", fill = "purple") +
  xlab("Group size") + ylab("Frequency")

dev.off()
