#################################################
#  File Name:volcano.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sun 30 May 2021 11:31:54 AM UTC
#################################################

library(ggpubr)
library(ggthemes)
args = commandArgs(T)
deg.data <- read.table(args[1],header=T,sep="\t")
deg.data$logP <- -log10(deg.data$padj)
#ggscatter(deg.data,x="logFC",y="logP")+theme_base()
deg.data$Group = "not-significant"
deg.data$Group[which((deg.data$padj < 0.01) & (deg.data$log2FoldChange > 3))] = paste0(args[2])
deg.data$Group[which((deg.data$padj < 0.01) & (deg.data$log2FoldChange < -3))] = paste0(args[3])
deg.data$Group=factor(deg.data$Group,levels = c(args[2],"not-significant",args[3]))
deg.data2 = deg.data[deg.data$Group==args[2]|deg.data$Group==args[3],]
write.table(deg.data2,file="significant_peak_xls",sep="\t",quote=FALSE)
mytable=as.data.frame(table(deg.data$Group))

p=ggscatter(deg.data, x = "log2FoldChange", y = "logP", color = "Group", 
	palette = c("#2f5688", "#BBBBBB", "#CC0000"), size = 1,
	#label = deg.data$Label, font.label = 8, repel = T, 
	xlab = "log2FoldChange", ylab = "-log10(Adjust P-value)",)+
theme_base()+geom_hline(yintercept = 2, linetype="dashed")+
geom_vline(xintercept = c(-3,-3), linetype="dashed")+annotate("text", x = 5, y = 10, label = mytable[1,2])+ 
annotate("text", x = 0, y = 10, label = mytable[2,2])+annotate("text", x = -5, y = 10, label = mytable[3,2])
ggsave("all_peak_diff2.pdf",p,width = 6,height = 4)
save.image('all_peak_diff2.Rdata')
