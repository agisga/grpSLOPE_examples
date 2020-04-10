library(tidyverse)
library(microbenchmark)
library(xtable)

load("mbm.RData")

ls()
print(mbm, unit = "s")
xtable(print(mbm, unit = "s"))

autoplot(mbm) + theme_minimal()
ggsave("runtimes.pdf")
