#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

set.seed(1225)
log <- file(args[4], open="wt")
sink(log, type="message")

data_dir <- args[1]
fragment_size <- args[2]
median_size <- args[3]

suppressPackageStartupMessages(library(cfDNAPro))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

setwd(data_dir)
grp_list <- as.list(dir())
names(grp_list) <- dir()

result<-sapply(grp_list, function(x){
  result <-callSize(path = data_dir) %>%
    dplyr::filter(group==as.character(x)) %>%
    plotSingleGroup()
}, simplify = FALSE)

plasma_l<-length(list.files(pattern = '.txt',path= paste0(data_dir,'/plasma')))
tumor_l<-length(list.files(pattern = '.txt',path= paste0(data_dir,'/tumor')))
control_l<-length(list.files(pattern = '.txt',path= paste0(data_dir,'/control')))

expr <- parse(text = paste0('result$plasma','$prop_plot'))
plasma <- eval(expr)
expr <- parse(text = paste0('result$tumor','$prop_plot'))
tumor <- eval(expr)
expr <- parse(text = paste0('result$control','$prop_plot'))
control <- eval(expr)

expr <- parse(text = paste0('result$plasma','$cdf_plot'))
plasma_cdf <- eval(expr)
expr <- parse(text = paste0('result$tumor','$cdf_plot'))
tumor_cdf <- eval(expr)
expr <- parse(text = paste0('result$control','$cdf_plot'))
control_cdf <- eval(expr)

suppressWarnings(
  multiplex <-
    ggarrange(plasma +
              theme(axis.title.x = element_blank()),
            tumor +
              theme(axis.title = element_blank()),
            control +
              theme(axis.title = element_blank()),
            plasma_cdf,
            tumor_cdf,
            control_cdf  +
              theme(axis.title.y = element_blank()),
            labels = c(paste0("Plasma cohort (n=",plasma_l,")"),
                       paste0("Tumor cohort (n=",tumor_l,")"),
                       paste0("Control cohort (n=",control_l,")")),
            label.x = 0.2,
            ncol = 3,
            nrow = 2))

tiff(fragment_size, width = 12, height = 7, units = 'in', res = 300, compression = 'lzw')
multiplex
dev.off()

# Set an order for those groups (i.e. the levels of factors).
order <- dir()
compare_grps<-callMetrics(data_dir) %>% plotMetrics(order=order)
#> setting default input_type to picard.

# Modify plots.
p1<-compare_grps$median_prop_plot +
  ylim(c(0, 0.028)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank())

p2<-compare_grps$median_cdf_plot +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank())

# Finalize plots.
suppressWarnings(
  median_grps<-ggpubr::ggarrange(p1,
                       p2,
                       label.x = 0.3,
                       ncol = 1,
                       nrow = 2
                       ))

tiff(median_size, width = 12, height = 7, units = 'in', res = 300, compression = 'lzw')
median_grps
dev.off()