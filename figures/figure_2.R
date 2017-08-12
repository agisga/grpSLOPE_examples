#---------------------------------------------------------------------------
# This code reproduces the figures from simulation study of Figure 2 in
# D. Brzyski, A.Gossmann, W. Su, M. Bogdan (2016), "Group SLOPE -
# - adaptive selection of groups of predictors"
# Preprint available at https://arxiv.org/abs/1610.04960
#---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# Load the simulation results for figures 2b-e
load("../RData/Figure_2bcde_3a_results.RData")

# get means and error bars for FDP and power
results.summary <- results %>%
  select(-replication, -starts_with("grp_size_")) %>%
  group_by(n.relevant) %>%
  summarize_all(funs(mean, lwr = (mean(.) - 2*sd(.)/sqrt(n.replications)),
                     upr = (mean(.) + 2*sd(.)/sqrt(n.replications))))


#--- Figure 2b-c

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

ggsave(file = "Figure_2bc.png", width = 12, height = 9.6, units = "cm")

#--- Figure 2d

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

ggsave(file = "Figure_2d.png", width = 12, height = 9.6, units = "cm")
