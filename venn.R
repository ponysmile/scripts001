a1=read.csv("Results/apoptosis_pvalue<0.05.csv")
a2=read.csv("../GSE16798_PIR_U937/Results/apoptosis_pvalue<0.01.csv")
library(VennDiagram)
A=a1$symbol
B=a2$symbol
venn.plot <- venn.diagram(
  list(A=A,B=B),
  filename = "apoptosis.tiff",
  lty = 1,
  lwd = 1,
  col = "black", 
  fill = c("coral", "cornflowerblue"),
  alpha = 0.60,
  cat.col = "black",
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
  )
a=intersect(A,B)
a3=a1[a1$symbol%in%a,]
a4=a2[a2$symbol%in%a,]
