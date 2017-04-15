# Load libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)

# Get command line argument
args <- commandArgs(TRUE)

# Set working directory
setwd(dirname(args[1]))

# Load Data
df <- read.csv(args[1], header = T, sep = ",")

# Get number of samples
num_samp = length(unique(df$sample))

### Create yRNA length distribution on sample by sample basis

# Number of columns for facet_wrap
fw_col <- length(unique(df$yRNA))

# Make graph
png("yRNA_length_distribution.png", width = 10, height = 3 * num_samp, units = "in", res = 200)
gg <- ggplot(df, aes(x = length, y = len_per, fill = prevalence)) +
        geom_bar(stat="identity") +
        ylim(0,100) +
        facet_wrap(~yRNA, ncol = fw_col) +
        scale_fill_gradientn(colours = c('#289fd6', '#f4ee42', '#f45041'), limits=c(0, 100)) +
        theme_bw() +
        labs(list(title="Length Distribution", x="Read length", y="Ratio of total reads")) 

lp <- dlply(df, "sample", function(d) gg %+% d + ggtitle(unique(d$sample)))
grid.arrange(grobs=lp, ncol=1)

invisible(dev.off())

### Create yRNA distribution across samples

# Get colors
color_set = c('#FF8360', '#63A375', '#E8E288', '#3CDBD3')

# Make graph
png('yRNA_distribution.png', unit='in', width = num_samp + 2, height = 6, res = 200)
ggplot(df, aes(sample, prevalence, fill = yRNA)) +
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_manual(values=color_set) +
  labs(list(title="yRNA Distribution", x="Sample", y="Ratio of total yRNA reads")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
invisible(dev.off())

