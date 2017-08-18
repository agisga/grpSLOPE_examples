#---------------------------------------------------------------------------
# This code reproduces the figures from simulation study of Figures 3a-c in
# D. Brzyski, A.Gossmann, W. Su, M. Bogdan (2016), "Group SLOPE -
# - adaptive selection of groups of predictors"
# Preprint available at https://arxiv.org/abs/1610.04960
#---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# palette to be used in all plots
library(RColorBrewer)
pal = brewer.pal(3, "Set1")

# TODO: legend labels to be used in all plots
math.labs <- list(bquote(list(q == 0.05, lambda^.("max"))),
                  bquote(list(q == 0.1, lambda^.("max"))),
                  bquote(list(q == 0.05, lambda^.("mean"))),
                  bquote(list(q == 0.1, lambda^.("mean"))))

#--- Figures 3 a,b,c

# load the simulation results
load("../RData/figure_3_simulation_results.RData")
results <- rename(results, n_relevant = n.relevant)

# prepare data for plotting
fractions_df <- results %>%
  group_by(n_relevant) %>%
  summarize_all(sum) %>% tbl_df() %>%
  gather(grp_size, count, -n_relevant) %>%
  mutate(weight = gsub("^(.+)_wt_\\d$", "\\1", grp_size)) %>%
  mutate(grp_size = gsub("^.+_wt_(\\d)$", "\\1", grp_size)) %>%
  group_by(n_relevant, weight) %>%
  mutate(total_count_by_n_relevant = sum(count)) %>%
  # TODO: figure out what's going on with `len`...
  filter(weight != "len") %>%
  mutate(fraction = count / total_count_by_n_relevant)

ggplot(fractions_df, aes(n_relevant, fraction, color = grp_size)) +
  geom_line() + geom_point() + facet_wrap(~weight) + ylim(0, 0.45)
