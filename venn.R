a1=read.csv("Results/apoptosis_pvalue<0.05.csv")
a2=read.csv("../GSE16798_PIR_U937/Results/apoptosis_pvalue<0.01.csv")
library(VennDiagram)
A=a1$symbol
B=a2$symbol
venn.plot <- venn.diagram(
  list(A=A,B=B),
  filename = "apoptosis.tiff",##韦恩图的名字
  lty = 1,
  lwd = 1,
  col = "black",  ##圈的颜色
  fill = c("coral", "cornflowerblue"),##对应每个圈的颜色，有几个数据集，就需要有相应数量的颜色
  alpha = 0.60,
  cat.col = "black",##此处设置每个数据集的名称颜色，也可以使用c（）函数输入三种颜色
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
  )
a=intersect(A,B)
a3=a1[a1$symbol%in%a,]
a4=a2[a2$symbol%in%a,]
