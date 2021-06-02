#################################################
#  File Name:venn.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Wed 24 Mar 2021 02:02:20 PM UTC
#################################################

library(VennDiagram)
args=commandArgs(T)
name=paste0(args[1],".pdf")
pdf(name,height=4,width=4)
left=as.numeric(args[2])
right=as.numeric(args[3])
ovlp=as.numeric(args[4])
left_name=args[5]
right_name=args[6]
draw.pairwise.venn(left,right,ovlp,c(left_name,right_name),col=c("black","blue"))
dev.off()
