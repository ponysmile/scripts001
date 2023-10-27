library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
load("Data/mhh_RNASeq_DEG.Rdata")
data=DEG_edgeR #limma DEG for example
# data$Genes=rownames(data)
colnames(data)
# "logFC" "AveExpr" "t" "P.Value" "adj.P.Val" "B" "change" "Genes" 
# colnames(data)[4]="PValue"
# colnames(data)[5]="FDR"
# colnames(data)
p1 <- ggplot(data, aes(logFC, -log(PValue,10))) + # -log10 conversion  
  theme_bw() +
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"PValue"))
p1
data$Expression = case_when(data$logFC >= log(2) & data$PValue <= 0.05 ~ "Up-regulated",
                            data$logFC <= -log(2) & data$PValue <= 0.05 ~ "Down-regulated",
                            TRUE ~ "Unchanged")
p2 <- ggplot(data, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  theme_bw() +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"PValue")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

data$Significance = case_when(
  abs(data$logFC) >= log(2) & data$PValue <= 0.05 & data$PValue > 0.01 ~ "PValue 0.05", 
  abs(data$logFC) >= log(2) & data$PValue <= 0.01 & data$PValue > 0.001 ~ "PValue 0.01",
  abs(data$logFC) >= log(2) & data$PValue <= 0.001 ~ "PValue 0.001", 
  TRUE ~ "Unchanged")

p3 <- ggplot(data, aes(logFC, -log(PValue,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  theme_bw() +
  ylab(expression("-log"[10]*"PValue")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3
options(ggrepel.max.overlaps = Inf)
top <- 10
top_genes <- bind_rows(
  data %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(PValue, desc(abs(logFC))) %>% 
    head(top),
  data %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(PValue, desc(abs(logFC))) %>% 
    head(top)
)
p4 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log(PValue,10), label = SYMBOL),
                   size = 2)
p4
pdf("Results/mhh_volcano.pdf",width = 6,height = 4)
print(p4)
dev.off()

top <- 1
# top_genes <- bind_rows(
#   data %>% 
#     filter(Expression == 'Up-regulated') %>% 
#     arrange(PValue, desc(abs(logFC))) %>% 
#     head(top))

top_genes <- bind_rows(
  data %>%
    filter(SYMBOL == 'PIR') %>%
    arrange(PValue, desc(abs(logFC))) %>%
    head(top))
p5 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log(PValue,10), label = SYMBOL),
                   size = 2)
p5

pdf("Results/mhh_volcano1.pdf",width = 5,height = 4)
print(p5)
dev.off()

top <- 6
# top_genes <- bind_rows(
#   data %>% 
#     filter(Expression == 'Up-regulated') %>% 
#     arrange(PValue, desc(abs(logFC))) %>% 
#     head(top))
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
p6

pdf("Results/mhh_volcano2.pdf",width = 5,height = 4)
print(p6)
dev.off()
