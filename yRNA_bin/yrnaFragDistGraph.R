# Load necessary libraries
library(ggplot2)
library(plyr)
library(gridExtra)

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Load yRNA fragment distribution data
df <- read.csv(args[1])

# Normalize so the location with the highest # of counts is 1
yRNAs = unique(df$yRNA)
for (sample in unique(df$sample)) {
  for (yRNA in yRNAs) {
    max_c = max(df[df$yRNA == yRNA & df$sample == sample, 'count'])
    df$count = ifelse(df$yRNA == yRNA & df$sample == sample, df$count / max_c, df$count)
  }
}

# Graph fragment data with ggplot2

height_multi = length(unique(df$sample))

png("yRNA_fragment_distribution.png", width = 10, height = 3 * height_multi, units = "in", res = 200)
gg <- ggplot(df, aes(loc, count)) +
  geom_bar(stat='identity', width = 1) + 
  scale_x_continuous(name = 'Nucleotide position') +
  scale_y_continuous(name = 'Fragment loc',
                     breaks = c(0,.25,.5,.75,1)) +
  facet_wrap(~yRNA, nrow = 1, scales = "free_x") +
  theme_bw()

lp <- dlply(df, "sample", function(d) gg %+% d + ggtitle(unique(d$sample)))
grid.arrange(grobs=lp, ncol=1)
invisible(dev.off())
