library(tidyverse)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(rtracklayer)
rm(list=ls())
if(F){
### ssGSEA ######
## table S1 - https://doi.org/10.1016/j.immuni.2013.10.003 
## pdf -> table -> read
immunity <- read.csv("immune_mhh/immunity-cell-gene.csv", header = T)
colnames(immunity)
immunity=immunity[,c(1,5,9,13,20)]
immunity$ENTREZ_GENE_ID=stringr::str_split(immunity$EntrezGene..Name," ",simplify = T)[,1]
immunity$Gene.Symbol=immunity$Symbol
head(immunity)
#   CellType AffymetrixID Symbol Gene.Symbol ENTREZ_GENE_ID
# 1      aDC    205569_at  LAMP3       LAMP3          27074
# 2      aDC    207533_at   CCL1        CCL1           6346
# 3      aDC    210029_at   INDO        IDO1           3620
# 4      aDC    218400_at   OAS3        OAS3           4940
# 5      aDC    219424_at   EBI3        EBI3          10148
# 6  B cells    204836_at   GLDC        GLDC           2731

idx <- !immunity$CellType %in% c("Blood vessels", "Normal mucosa", "SW480 cancer cells", "Lymph vessels")
immunity <- immunity[idx,]
immunity <- immunity %>%
  split(., .$CellType) %>%
  lapply(., function(x)(x$ENTREZ_GENE_ID))
immunity <- lapply(immunity, unique)

## Ensembl download
anno <- import('immune_mhh/Homo_sapiens.GRCh38.108.gtf')
anno <- as.data.frame(anno)
anno <- anno[!duplicated(anno$gene_id),]

library("AnnotationDbi")
library("org.Hs.eg.db")
gene_symbol = anno[,c("gene_id","gene_name")]
columns(org.Hs.eg.db)
gene_symbol$ENTREZID <-mapIds(org.Hs.eg.db, keys = gene_symbol$gene_id, keytype="ENSEMBL", column = "ENTREZID")
# gene_symbol$SYMBOL <-mapIds(org.Hs.eg.db, keys = gene_symbol$gene_id, keytype="ENSEMBL", column = "SYMBOL")

anno <- merge(anno, gene_symbol, by = "gene_id")
# anno <- rbind(anno, data.frame(gene_name = c("KIAA1324", "IGHA1"),
#                                gene_id = c("ENSG00000116299", "ENSG00000211895"),
#                                ENTREZID = c("57535", "3492")))
# anno <- anno[!duplicated(anno$gene_id),] ### 37417
anno <- anno[, c("gene_id", "ENTREZID","gene_name.x")]
colnames(anno) <- c("gene_id", "ENTREZID","SYMBOL")
save(anno,gene_symbol,immunity,file = "immune_mhh/immune.Rdata")
}
rm(list=ls())
load("immune_mhh/immune.Rdata")
load("data/explore.Rdata")
# data <- fread("~/tpm.txt") %>% 
#   rename("gene_id" = "V1") %>% 
#   left_join(., anno, by = "gene_id") %>% 
#   filter(!is.na(ENTREZID)) %>% 
#   select(-gene_id) %>% 
#   column_to_rownames("ENTREZID")
data <- exp
data=as.data.frame(data)
data$SYMBOL=rownames(data)
data <-  data %>%
  left_join(., anno, by = "SYMBOL") %>% 
  filter(!is.na(ENTREZID)) %>% 
  filter(!duplicated(ENTREZID)) %>% 
  dplyr::select(-gene_id) %>% 
  dplyr::select(-SYMBOL) %>% 
  column_to_rownames("ENTREZID")


data <- log2(data + 1) 

immu_cell <-  as.data.frame(gsva(as.matrix(data), immunity, method = "ssgsea"))


### group plot 
# data <- read.table("~/file.txt", header = T)
#        group       aDC
# 1 1.13092315 0.4709550
# 3 0.55644003 0.1800251
# 4 0.44696904 0.3350859
# 5 0.05474605 0.1191767
# 7 0.61364297 0.1563856
# 8 0.41079217 0.4588979
# colnames(data) <- c("x", "y")
data <- as.data.frame(exp["PIR",])
colnames(data)="PIR"
data$group <- ifelse(data$PIR >= mean(data$PIR), "High", "Low")
data$group <- factor(data$group, levels = c("Low", "High"))
table(data$group)
immune=t(immu_cell)
identical(rownames(immune),rownames(data))
data=cbind(data,immune)
colnames(data)
ggplot(data, aes(x = group, y = `CD8 T cells`, color = group, fill = group)) +
  geom_boxplot(alpha = 0.2) +
  geom_point(position = position_jitter(0.3)) +
  theme_bw()

colnames(data)
library("ggpubr")
p <- ggboxplot(data, x = "group", y = "Cytotoxic cells",
               color = "group", palette =c("#00AFBB", "#E7B800"), # "#FC4E07"
               add = "jitter") #boxplot #jitter
p
my_comparisons <- list(  c("1", "2") )
p1=p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.3)                   # Add global p-value
p1
colnames(data)
if(F){
  for (i in 3:ncol(data)) {
    a=colnames(data)[i];a
    p <- ggboxplot(data, x = "group", y = a,
                   color = "group", palette =c("#00AFBB", "#E7B800"), # "#FC4E07"
                   add = "jitter") #boxplot #jitter
    p
    my_comparisons <- list(  c("1", "2") )
    p1=p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means()                   # Add global p-value
    p1
    pdf(file = paste0("immune_mhh/Results/PIR/",a,".pdf"),width=4,height=4)
    print(p1)
    dev.off()
  }
  
}

###FAS
dat <- as.data.frame(exp["FAS",])
colnames(dat)="FAS"
dat$group <- ifelse(dat$FAS >= mean(dat$FAS), "High", "Low")
dat$group <- factor(dat$group, levels = c("Low", "High"))
table(dat$group)
immune=t(immu_cell)
identical(rownames(immune),rownames(dat))
dat=cbind(dat,immune)
colnames(dat)
ggplot(dat, aes(x = group, y = `CD8 T cells`, color = group, fill = group)) +
  geom_boxplot(alpha = 0.2) +
  geom_point(position = position_jitter(0.3)) +
  theme_bw()

colnames(data)
library("ggpubr")
p <- ggboxplot(dat, x = "group", y = "Cytotoxic cells",
               color = "group", palette =c("#00AFBB", "#E7B800"), # "#FC4E07"
               add = "jitter") #boxplot #jitter
p
my_comparisons <- list(  c("1", "2") )
p1=p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.3)                   # Add global p-value
p1
colnames(dat)
if(F){
  for (i in 3:ncol(dat)) {
    a=colnames(dat)[i];a
    p <- ggboxplot(dat, x = "group", y = a,
                   color = "group", palette =c("#00AFBB", "#E7B800"), # "#FC4E07"
                   add = "jitter") #boxplot #jitter
    p
    my_comparisons <- list(  c("1", "2") )
    p1=p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means()                   # Add global p-value
    p1
    pdf(file = paste0("immune_mhh/Results/FAS/",a,".pdf"),width=4,height=4)
    print(p1)
    dev.off()
  }
  
}

if(F){
  dat$PIR=exp["PIR",]
  dat$group1 <- ifelse(dat$PIR >= mean(dat$PIR), "High", "Low")
  dat$group1 <- factor(dat$group1, levels = c("Low", "High"))
  table(dat$group1)
  table(dat$group)
  dat$Group="NOT"
  dat$Group <- ifelse(dat$group1=="High"&dat$group=="High","HH",
                      ifelse(dat$group1=="High"&dat$group=="Low", "HL", 
                             ifelse(dat$group1=="Low"&dat$group=="High", "LH", "LL")))
  table(dat$Group)
  library("ggpubr")
  p <- ggboxplot(dat, x = "Group", y = "Cytotoxic cells",
                 color = "Group", palette =c("#00AFBB", "#E7B800","#FC4E07","#07b5fc"), # "#FC4E07"
                 add = "jitter") #boxplot #jitter
  p
  my_comparisons <- list(c("LL","HH"),c("HL","LH") )
  p1=p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 0.3)                   # Add global p-value
  p1
  p1=p + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  p1
  colnames(dat)
  for (i in 3:27) {
    a=colnames(dat)[i];a
    p <- ggboxplot(dat, x = "Group", y = a,
                   color = "Group", palette =c("#00AFBB", "#E7B800","#FC4E07","#07b5fc"), # "#FC4E07"
                   add = "jitter") #boxplot #jitter
    p
    my_comparisons <- list(c("LL","HH"),c("HL","LH"))
    p1=p + stat_compare_means(comparisons = my_comparisons)# Add pairwise comparisons p-value
    pdf(file = paste0("immune_mhh/Results/PIR_FAS/",a,".pdf"),width=4,height=4)
    print(p1)
    dev.off()
  }
}
