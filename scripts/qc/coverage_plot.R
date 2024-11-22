#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)
log <- file(args[4], open="wt")
sink(log, type="message")

data_dir <- args[1]
plot_dist_single <- args[2]
plot_dist_grouped <- args[3]

suppressPackageStartupMessages(library(tidyverse))

setwd(data_dir)

con <- list.files(pattern = '.mosdepth.global.dist.txt')

df1<-data.frame()
for (i in con) {
tmp<- read_tsv(i,col_names = FALSE,show_col_types = FALSE)
tmp_sum<- tmp %>%
group_by(X1) %>%
summarise(average_coverage = mean(X3))
tmp_sum$sample<-gsub('.mosdepth.global.dist.txt','',i)
df1<-rbind(df1,tmp_sum)
}

df2<-df1[!df1$X1 %in% c('chrX','chrY','chrM','total'),]
df2$X1 <- factor(df2$X1, levels = paste0("chr", 1:22))
df2$type<- gsub('.*-','',df2$sample)

p_grouped <- ggplot(df2, aes(x = X1, y = average_coverage, color = sample, group = sample)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 2) + # Add points to highlight values
  theme_minimal() +      # Minimal theme for a clean look
  labs(
    title = "Coverage per Contig",
    x = "Region",
    y = "Average Coverage",
    color = "Sample"      # Legend title
  ) +
theme(
    strip.text = element_text(size = 12, face = "bold"), # Bold and large facet title
    axis.title = element_text(size = 12, face = "bold") # Bold and large x and y axis labels
  )

p_single<- ggplot(df2, aes(x = X1, y = average_coverage, color = sample, group = sample)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 2) + # Add points to highlight values
  theme_minimal() +      # Minimal theme for a clean look
  labs(
    title = "Coverage per Contig",
    x = "Region",
    y = "Average Coverage",
    color = "Sample"      # Legend title
  ) +
  facet_wrap(~ type, ncol = 1) +
theme(
    strip.text = element_text(size = 12, face = "bold"), # Bold and large facet title
    axis.title = element_text(size = 12, face = "bold") # Bold and large x and y axis labels
  )

tiff(plot_dist_single, width = 12, height = 7, units = 'in', res = 300, compression = 'lzw')
p_single
dev.off()

tiff(plot_dist_grouped, width = 12, height = 7, units = 'in', res = 300, compression = 'lzw')
p_grouped
dev.off()