#!/usr/bin/env Rscript
## Purpose of script: alpha diversity boxplot
## Date created: 22.07.2020
## Author: luigui gallardo-becerra (bfllg77@gmail.com)
# Package installation
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("ggpubr")
library("dplyr")
library("ggplot2")
library("ggpubr")
# Data input
alpha_data <- read.table(file = "gut_all_metrics.txt", sep = "\t", header = T)
# Define colors
colors <- c("royalblue", "orangered")
p_value <- list(c("Genetics 1", "Genetics 2"))
long_format <- stack(tab)
levels(long_format[,2]) <- c("Genetics 1","Genetics 2")
boxplot <- ggplot(long_format, aes(x = ind, y = values, color = ind)) + geom_boxplot() + theme_test() + theme(legend.position = "none") + labs(x = "", y = "Chao1 Index") + scale_color_manual(values = colors) + theme(text = element_text(size = 12), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14)) + stat_compare_means(comparisons = p_value, method = "wilcox.test")
boxplot
# Creation of pdf
pdf(file = "chao1_gut_genetics.pdf", width = 5, height = 5)
print (boxplot)
dev.off()






# Luigui's plot
### Shannon ###
boxplot <- ggplot(alpha_data, aes(x = group, y = shannon)) +
geom_boxplot(aes(color = group)) +
theme_test() +
theme(legend.position = "none") +
labs(x = "",
    y = "Shannon Index") +
scale_color_manual(values = colors) +
theme(text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14)) +
stat_compare_means(comparisons = p_value, 
#    label = "p.signif",
    method = "wilcox.test")
# Creation of pdf
pdf(file = "shannon_gut_genetics.pdf", width = 5, height = 5)
print (boxplot)
dev.off()

# Some other plot from the internet that may be useful:
set.seed(123)
#test df
library(tidyverse)
library(rstatix)
mydf <- data.frame(ID=paste(sample(LETTERS, 163, replace=TRUE), sample(1:1000, 163, replace=FALSE), sep=''),
                   Group=c(rep('C',10),rep('FH',10),rep('I',19),rep('IF',42),rep('NA',14),rep('NF',42),rep('NI',15),rep('NS',10),rep('PGMC4',1)),
                   Value=c(runif(n=100), runif(63,max= 0.5)))

stat_pvalue <- mydf %>% 
 rstatix::wilcox_test(Value ~ Group) %>%
 filter(p < 0.05) %>% 
 rstatix::add_significance("p") %>% 
 rstatix::add_y_position() %>% 
 mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

ggplot(mydf, aes(x=Group, y=Value)) + geom_boxplot() +
  ggpubr::stat_pvalue_manual(stat_pvalue, label = "p.signif") +
  theme_bw(base_size = 16)
