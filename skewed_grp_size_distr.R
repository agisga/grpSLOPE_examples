# This code aims to answer the following question:
#
# > Have you checked the performance of lambda_mean in case when the numbers of groups of different sizes are substantially different ?
# > Say we have 400 groups of size 3, 200 groups of size 4, 100 groups of size 5,  50 groups of size 6 and 10 groups of size 7.
#

library(grpSLOPE)
library(dplyr)
library(tidyr)
library(ggplot2)

# Adjust the number of cores to the particular system
library(doParallel)
registerDoParallel(cores = as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")))

#--- one simulation (to be repeated many times at various sparsity levels)
one.iteration <- function(n.signif, X, group, fdr, signal.strength){
  # generate coeffient vector, pick relevant groups at random
  b <- rep(0, p)
  ind.relevant <- sample(1:n.group, n.signif)
  for (j in ind.relevant) { b[group.id[[j]]] <- a }

  # generate the response vector
  y <- X %*% b + rnorm(p, sd=1)

  # get Group SLOPE solution
  b.grpSLOPE.mean <- grpSLOPE(X=X, y=y, group=group, fdr=fdr,
                              lambda="chiOrthoMean", sigma=1, verbose=FALSE,
                              orthogonalize=FALSE, normalize=FALSE)

  # FDR and power
  true.relevant <- names(group.id)[ind.relevant]
  truepos <- length(intersect(b.grpSLOPE.mean$selected, true.relevant))
  falsepos <- length(b.grpSLOPE.mean$selected) - truepos
  FDP <- falsepos / max(1, length(b.grpSLOPE.mean$selected))
  pow <- truepos / length(true.relevant)

  return(c("FDP" = FDP, "pow" = pow))
}
#---


# set X to identity matrix, b/c gSLOPE with orthogonal design is equivalent to a problem with identity mat.
p <- 2870 # 400*3 + 200*4 + 100*5 + 50*6 + 10*7
X <- diag(rep(1,p))

# set the grouping structure
group <- c(rep(1:400, each=3),
           rep(401:600, each=4),
           rep(601:700, each=5),
           rep(701:750, each=6),
           rep(751:760, each=7))
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
n.group <- length(group.id)

# determine signal strength, such as used in Figure 1 of Brzyski et. al. (2015)
Bfun <- function(l) {
  sqrt(4*log(n.group) / (1 - n.group^(-2/l)) - l)
}
a <- sum(Bfun(group.length)) / sum(sqrt(group.length))

# target FDR
fdr <- 0.1

# how many times the simulation is repeated
n.replications <- 50 

# considered numbers of truly relevant groups
n.relevant <- floor(seq(1, 250, length=11))

# run the simulations
set.seed(20160804)
parallel.results <- vector(mode="list")
for (k in 1:length(n.relevant)) {
  print(paste("sparsity level", k, "started..."))
  parallel.results[[k]] <- foreach(i=1:n.replications, .combine = rbind) %dopar% {
    one.iteration(n.relevant[k], X, group, fdr, signal.strength=a)
  }
  print(paste("...sparsity level", k, "done"))
}

# collect results in a data frame
results <- data.frame() %>% tbl_df
for(k in 1:length(n.relevant)) {
  results.k <- tbl_df(parallel.results[[k]])
  results.k.sub <- results.k %>% select(FDP, pow) %>%
    mutate(n.relevant = rep(n.relevant[k], n.replications), repl = 1:n.replications)
  results <- bind_rows(results, results.k.sub)
}

# get means and error bounds for FDR and power
results.summary <- results %>% group_by(n.relevant) %>% 
  summarize(FDR = mean(FDP), avg.pow = mean(pow), FDR.se = (sd(FDP)/sqrt(n.replications)),
            pow.se = (sd(pow)/sqrt(n.replications))) %>%
  mutate(FDR.lower = (FDR - 2*FDR.se), FDR.upper = (FDR + 2*FDR.se),
         pow.lower = (avg.pow - 2*pow.se), pow.upper = (avg.pow + 2*pow.se))

# pre-process for plotting
FDR.and.pow <- results.summary %>% select(n.relevant, FDR, avg.pow) %>% 
  rename(Power = avg.pow) %>% gather(key, value, -n.relevant)
lwr <- results.summary %>% select(n.relevant, FDR.lower, pow.lower) %>% 
  rename(FDR = FDR.lower, Power = pow.lower) %>% gather(key, lwr, -n.relevant)
upr <- results.summary %>% select(n.relevant, FDR.upper, pow.upper) %>% 
  rename(FDR = FDR.upper, Power = pow.upper) %>% gather(key, upr, -n.relevant)
FDR.and.pow <- FDR.and.pow %>% left_join(lwr) %>% left_join(upr)

# save everything
save(list=ls(all.names=TRUE), file="skewed_grp_size_distr.RData")

# plot estimated FDR and power
pdf(file="FDR_and_pow.pdf")

ggplot(FDR.and.pow) +
  geom_segment(mapping = aes(x = 0, xend = tail(n.relevant, 1), y = fdr, 
                             yend = fdr*(n.group - tail(n.relevant, 1))/n.group), 
               linetype = 2, color = "black") +
  geom_line(mapping = aes(x = n.relevant, y = value, color = key)) + 
  geom_point(mapping = aes(x = n.relevant, y = value, color = key)) +
  geom_ribbon(mapping = aes(x = n.relevant, ymin = lwr, ymax = upr, fill = key), alpha=0.25) +
  xlab("Num. significant groups") + ylab("Estimated FDR and Power") +
  guides(fill=guide_legend(title="+-2*SE"), color=guide_legend(title="")) +
  coord_cartesian(ylim = c(0, 1))

dev.off()
