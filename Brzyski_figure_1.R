#---------------------------------------------------------------------------
# This code is based on the simulation study of Figure 1 in
# D. Brzyski, W. Su, M. Bogdan (2015), "Group SLOPE - adaptive selection of 
# groups of predictors" (http://arxiv.org/abs/1511.09078)
#---------------------------------------------------------------------------

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

# set X to identity matrix, b/c gSLOPE with orthogonal design is equivalent to a problem with identity mat.
p <- 5000
X <- diag(rep(1,p))

# set the grouping structure
group <- c(rep(1:200, each=3),
           rep(201:400, each=4),
           rep(401:600, each=5),
           rep(601:800, each=6),
           rep(801:1000, each=7))
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
n.group <- length(group.id)

# determine signal strength, such as used in Figure 1 in Brzyski et. al. (2015)
Bfun <- function(l) {
  sqrt(4*log(n.group) / (1 - n.group^(-2/l)) - l)
}
signal.strength <- sum(Bfun(group.length)) / sum(sqrt(group.length))

# considered numbers of truly relevant groups
n.relevant <- floor(seq(1, 250, length=11))

# how many times the simulation is repeated
n.replications <- 300 


#--- Simulation main loop

# run the simulations n.replications times at each level of n.relevant,
# and with each combination of lambda="max", lambda="mean", fdr=0.1, fdr=0.05
set.seed(20160807)
parallel.results <- vector(mode="list")
for (k in 1:length(n.relevant)) {
  cat(paste("sparsity level", k, "started"))

  parallel.results[[k]] <- foreach(i=1:n.replications, .combine = rbind) %dopar% {
    cat(".")

    # generate coeffient vector, pick relevant groups at random
    b <- rep(0, p)
    n.signif <- n.relevant[k]
    ind.relevant <- sample(1:n.group, n.signif)
    for (j in ind.relevant) { b[group.id[[j]]] <- signal.strength }

    # generate the response vector
    y <- X %*% b + rnorm(p, sd=1)

    # get Group SLOPE solutions with different lambda and fdr values
    lambda.max.1 <- grpSLOPE(X=X, y=y, group=group, fdr=0.1,
                             lambda="max", sigma=1, verbose=FALSE,
                             orthogonalize=FALSE, normalize=FALSE)
    lambda.max.05 <- grpSLOPE(X=X, y=y, group=group, fdr=0.05,
                              lambda="max", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)
    lambda.mean.1 <- grpSLOPE(X=X, y=y, group=group, fdr=0.1,
                              lambda="mean", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)
    lambda.mean.05 <- grpSLOPE(X=X, y=y, group=group, fdr=0.05,
                               lambda="mean", sigma=1, verbose=FALSE,
                               orthogonalize=FALSE, normalize=FALSE)

    # get the FDPs and powers of the grpSLOPE solutions
    true.relevant <- names(group.id)[ind.relevant]
    FDPs.and.powers <- rep(NA, 10)
    names(FDPs.and.powers) <- c("n.relevant", "replication",
                                "lambda.max.1.FDP", "lambda.max.1.power",
                                "lambda.max.05.FDP", "lambda.max.05.power",
                                "lambda.mean.1.FDP", "lambda.mean.1.power",
                                "lambda.mean.05.FDP", "lambda.mean.05.power")
    FDPs.and.powers["n.relevant"] <- length(true.relevant)
    FDPs.and.powers["replication"] <- i
    FDPs.and.powers[c("lambda.max.1.FDP", "lambda.max.1.power")] <- get_FDP_and_power(lambda.max.1, true.relevant)
    FDPs.and.powers[c("lambda.max.05.FDP", "lambda.max.05.power")] <- get_FDP_and_power(lambda.max.05, true.relevant)
    FDPs.and.powers[c("lambda.mean.1.FDP", "lambda.mean.1.power")] <- get_FDP_and_power(lambda.mean.1, true.relevant)
    FDPs.and.powers[c("lambda.mean.05.FDP", "lambda.mean.05.power")] <- get_FDP_and_power(lambda.mean.05, true.relevant)

    FDPs.and.powers
  }

  cat("done\n")
}
save(parallel.results, file = "./RData/Brzyski_1_results.RData") # keep a backup


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


#--- Plot a figure analogous to Figure 1a-b in Brzyski et. al. (2015)

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
pdf(file = "./figures/Brzyski_1ab.pdf")

xend <- tail(n.relevant, 1)
ggplot(FDR.results) +
  geom_segment(mapping = aes(x = 0, xend = xend, y = 0.1, 
                             yend = 0.1*(n.group - xend)/n.group), 
               linetype = 2, color = "black") +
  geom_segment(mapping = aes(x = 0, xend = xend, y = 0.05, 
                             yend = 0.05*(n.group - xend)/n.group), 
               linetype = 2, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = FDR, color = scenario)) + 
  geom_point(mapping = aes(x = n.relevant, y = FDR, color = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, 
                            fill = scenario), alpha=0.25) +
  xlab("Number of significant groups") + ylab("Estimated gFDR +/- 2SE") +
  scale_color_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                       breaks = c("lambda.max.05", "lambda.max.1",
                                  "lambda.mean.05", "lambda.mean.1"),
                       labels=c("'max', 0.05", "'max', 0.1", 
                                "'mean', 0.05", "'mean', 0.1")) +
  scale_fill_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                      breaks = c("lambda.max.05", "lambda.max.1",
                                 "lambda.mean.05", "lambda.mean.1"),
                      labels=c("'max', 0.05", "'max', 0.1", 
                               "'mean', 0.05", "'mean', 0.1")) +
  coord_cartesian(ylim = c(0, 0.15))

dev.off()
 

#--- Plot a figure analogous to Figure 1c in Brzyski et. al. (2015)

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
pdf(file = "./figures/Brzyski_1c.pdf")

xend <- tail(n.relevant, 1)
ggplot(power.results) +
  geom_line(mapping = aes(x = n.relevant, y = power, color = scenario)) + 
  geom_point(mapping = aes(x = n.relevant, y = power, color = scenario)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, 
                            fill = scenario), alpha=0.25) +
  xlab("Number of significant groups") + ylab("Estimated power +/- 2SE") +
  scale_color_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                       breaks = c("lambda.max.05", "lambda.max.1",
                                  "lambda.mean.05", "lambda.mean.1"),
                       labels=c("'max', 0.05", "'max', 0.1", 
                                "'mean', 0.05", "'mean', 0.1")) +
  scale_fill_discrete(name = "Parameters:\nlambda,\ntarget gFDR", 
                      breaks = c("lambda.max.05", "lambda.max.1",
                                 "lambda.mean.05", "lambda.mean.1"),
                      labels=c("'max', 0.05", "'max', 0.1", 
                               "'mean', 0.05", "'mean', 0.1")) +
  coord_cartesian(ylim = c(0.3, 1))

dev.off()
