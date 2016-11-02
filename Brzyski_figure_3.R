#--------------------------------------------------------------------------
# This code is based on the simulation study of Figure 3 in
# D. Brzyski, W. Su, M. Bogdan (2015), "Group SLOPE - adaptive selection 
# of groups of predictors" (http://arxiv.org/abs/1511.09078)
#--------------------------------------------------------------------------

library(grpSLOPE)
library(dplyr)
library(tidyr)
library(ggplot2)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores=4)


#--- Set up global parameters for the simulation

# auxilliary function to get (group-wise) FDP and power from one grpSLOPE solution object
get_FDP_and_power <- function(result, true.relevant){
  truepos <- length(intersect(result$selected, true.relevant))
  falsepos <- length(result$selected) - truepos
  FDP <- falsepos / max(1, length(result$selected))
  pow <- truepos / length(true.relevant)
  return(c("FDP" = FDP, "pow" = pow))
}

# number of predictors
p <- 5000

# set the grouping structure
n.group <- 1000
group <- rep(1:1000, each=5)
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)

# determine signal strength, such as used in Figure 1 in Brzyski et. al. (2015)
signal.strength <- sqrt(4*log(n.group) / (1 - n.group^(-2/5)) - 5)

# considered numbers of truly relevant groups
n.relevant <- floor(seq(1, 100, length=10))

# how many times the simulation is repeated
n.replications <- 200 


#--- Simulation main loop

# run the simulations n.replications times at each level of n.relevant,
# and with each combination of lambda="max", lambda="mean", fdr=0.1, fdr=0.05
set.seed(20160817)
parallel.results <- vector(mode="list")
for (k in 1:length(n.relevant)) {
  cat(paste("sparsity level", k, "started"))

  parallel.results[[k]] <- foreach(i=1:n.replications, .combine = rbind) %dopar% {
    cat(".")

    # model matrix
    X <- matrix(rnorm(p^2, sd=sqrt(1/p)), p, p)

    # generate coeffient vector, pick relevant groups at random
    b <- rep(0, p)
    n.signif <- n.relevant[k]
    ind.relevant <- sample(1:n.group, n.signif)
    for (j in ind.relevant) { b[group.id[[j]]] <- signal.strength }
    for (j in ind.relevant) {
      X1 <- apply(X[ , group.id[[j]]], 1, sum)
      b[group.id[[j]]] <- (signal.strength / norm(as.matrix(X1), "f")) * rep(1, group.length[j])
    }

    # generate the response vector
    y <- X %*% b + rnorm(p, sd=1)

    # get Group SLOPE solutions with different lambda and fdr values
    lambda.max.1 <- grpSLOPE(X=X, y=y, group=group, fdr=0.1,
                             lambda="max", sigma=1, verbose=FALSE,
                             orthogonalize=FALSE, normalize=FALSE)
    lambda.max.05 <- grpSLOPE(X=X, y=y, group=group, fdr=0.05,
                              lambda="max", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)
    lambda.corrected.1 <- grpSLOPE(X=X, y=y, group=group, fdr=0.1,
                                   lambda="corrected", sigma=1, 
                                   verbose=FALSE, orthogonalize=FALSE, 
                                   normalize=FALSE)
    lambda.corrected.05 <- grpSLOPE(X=X, y=y, group=group, fdr=0.05,
                                    lambda="corrected", sigma=1, 
                                    verbose=FALSE, orthogonalize=FALSE, 
                                    normalize=FALSE)

    # get the FDPs and powers of the grpSLOPE solutions
    true.relevant <- names(group.id)[ind.relevant]
    FDPs.and.powers <- rep(NA, 10)
    names(FDPs.and.powers) <- c("n.relevant", "replication",
                                "lambda.max.1.FDP", "lambda.max.1.power",
                                "lambda.max.05.FDP", "lambda.max.05.power",
                                "lambda.corrected.1.FDP", 
                                "lambda.corrected.1.power",
                                "lambda.corrected.05.FDP", 
                                "lambda.corrected.05.power")
    FDPs.and.powers["n.relevant"] <- length(true.relevant)
    FDPs.and.powers["replication"] <- i
    FDPs.and.powers[c("lambda.max.1.FDP", "lambda.max.1.power")] <- get_FDP_and_power(lambda.max.1, true.relevant)
    FDPs.and.powers[c("lambda.max.05.FDP", "lambda.max.05.power")] <- get_FDP_and_power(lambda.max.05, true.relevant)
    FDPs.and.powers[c("lambda.corrected.1.FDP", "lambda.corrected.1.power")] <- get_FDP_and_power(lambda.corrected.1, true.relevant)
    FDPs.and.powers[c("lambda.corrected.05.FDP", "lambda.corrected.05.power")] <- get_FDP_and_power(lambda.corrected.05, true.relevant)

    FDPs.and.powers
  }

  cat("done\n")
}
save(parallel.results, file = "./RData/Brzyski_3_results.RData") # keep a backup


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


#--- Plot a figure analogous to Figure 3a-b in Brzyski et. al. (2015)

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
png(file = "./figures/Brzyski_3ab.png", width = 600, height = 480)

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
  scale_color_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                       breaks = c("lambda.max.05", "lambda.max.1",
                                  "lambda.corrected.05", 
                                  "lambda.corrected.1"),
                       labels=c("'max', 0.05", "'max', 0.1", 
                                "'corrected', 0.05", "'corrected', 0.1")) +
  scale_fill_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                      breaks = c("lambda.max.05", "lambda.max.1",
                                 "lambda.corrected.05", 
                                 "lambda.corrected.1"),
                      labels=c("'max', 0.05", "'max', 0.1", 
                               "'corrected', 0.05", "'corrected', 0.1")) +
  coord_cartesian(ylim = c(0, 0.3))

dev.off()
 
#--- Plot a figure analogous to Figure 3c in Brzyski et. al. (2015)

# plot lambda sequences
png(file = "./figures/Brzyski_3c.png", width = 600, height = 480)

weight <- sqrt(group.length)
lambda.corrected.1 <- lambdaGroupSLOPE(method = "corrected", fdr = 0.1, 
                                       group = group, wt = weight, n.obs = p)
lambda.corrected.05 <- lambdaGroupSLOPE(method = "corrected", fdr = 0.05, 
                                        group = group, wt = weight, n.obs = p)
lambda.max.1 <- lambdaGroupSLOPE(method = "max", fdr = 0.1, group = group, wt = weight)
lambda.max.05 <- lambdaGroupSLOPE(method = "max", fdr = 0.05, group = group, wt = weight)

lambdas.1 <- data.frame("gFDR" = rep(0.1, n.group), "Index" = 1:n.group, 
                        "corrected" = lambda.corrected.1, "max" = lambda.max.1) %>% 
             tbl_df %>% gather(Method, Lambda, -gFDR, -Index)
lambdas.05 <- data.frame("gFDR" = rep(0.05, n.group), "Index" = 1:n.group, 
                         "corrected" = lambda.corrected.05, "max" = lambda.max.05) %>% 
              tbl_df %>% gather(Method, Lambda, -gFDR, -Index)
lambdas <- bind_rows(lambdas.1, lambdas.05) %>% mutate(gFDR = as.factor(gFDR))

lambdas %>% filter(Index <= 100) %>%
  ggplot(aes(x = Index, y = Lambda, color = Method, linetype = gFDR)) + geom_line()

dev.off()


#--- Plot estimated power for a simulation analogous to Figure 3 in Brzyski et. al. (2015)

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
png(file = "./figures/Brzyski_3_power.png", width = 600, height = 480)

ggplot(power.results) +
  geom_line(mapping = aes(x = n.relevant, y = power, color = scenario)) + 
  geom_point(mapping = aes(x = n.relevant, y = power, color = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, 
                            fill = scenario), alpha=0.25) +
  xlab("Number of significant groups") + ylab("Estimated power +/- 2SE") +
  scale_color_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                       breaks = c("lambda.max.05", "lambda.max.1",
                                  "lambda.corrected.05", 
                                  "lambda.corrected.1"),
                       labels=c("'max', 0.05", "'max', 0.1", 
                                "'corrected', 0.05", "'corrected', 0.1")) +
  scale_fill_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                      breaks = c("lambda.max.05", "lambda.max.1",
                                 "lambda.corrected.05", 
                                 "lambda.corrected.1"),
                      labels=c("'max', 0.05", "'max', 0.1", 
                               "'corrected', 0.05", 
                               "'corrected', 0.1")) +
  coord_cartesian(ylim = c(0.5, 1))

dev.off()
