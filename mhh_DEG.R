### DEG
library(edgeR)
rm(list=ls())
load("data/01.data_prepare.Rdata")
Group
comp <- unlist(strsplit("shPIR_vs_ctrl",split = "_vs_"))
Group <- factor(Group,levels = c("shPIR","ctrl"))
Group
table(Group)
design <- model.matrix(~0+Group)
rownames(design) <- colnames(exp_filter)
colnames(design) <- levels(factor(comp))
design
DEG <- DGEList(counts=exp_filter, 
               group=Group)
DEG <- calcNormFactors(DEG)

DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

fit <- glmFit(DEG, design)
lrt <- glmLRT(fit, contrast=c(1,-1)) 
DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

fc_cutoff <- 2
pvalue <- 0.05

DEG_edgeR$regulated <- "stable"

loc_up <- intersect(which( DEG_edgeR$logFC > log2(fc_cutoff) ),
                    which( DEG_edgeR$PValue < pvalue) )

loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)
DEG_edgeR$SYMBOL=rownames(DEG_edgeR)
DEG_edgeR <- na.omit(DEG_edgeR)
head(DEG_edgeR)
table(DEG_edgeR$regulated)
save(DEG_edgeR,file = "Data/mhh_RNASeq_DEG.Rdata")
