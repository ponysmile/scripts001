library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
load("Data/mhh_RNASeq_DEG.Rdata")
data=DEG_edgeR 
p3 <- ggplot(data, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  theme_bw() +
  ylab(expression("-log"[10]*"PValue")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p3

a=c("PIR","FAS","DDIT3","CASP8","TNFRSF10B","STAT3")
top_genes <- bind_rows(
  data %>%
    filter(SYMBOL %in% a) %>%
    arrange(PValue, desc(abs(logFC))) %>%
    head(top))
p6 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log(PValue,10), label = SYMBOL),
                   size = 2)

