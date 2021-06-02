#
###############################################
#  File Name:corr.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 14 Nov 2020 11:00:07 AM UTC
#################################################
library(edgeR)
args=commandArgs(T)
data = read.table(args[1])
new = data[,4:5]
#filter=new[which(new$encode>20),]
new_norm_temp=cpm(new)
new_norm=log(new_norm_temp+1,base=2)
corr <- paste("Pearson Cor. = ", round(cor(new_norm[,2], new_norm[,1]), 2), sep="")
pdf(paste0(args[1],".pdf"),height=6,width=6)
#plot(new_norm, col = densCols(new_norm), pch = 20,xlim=c(0,8),ylim=c(0,8))
plot(new_norm[,1],new_norm[,2],pch=16,cex=0.7,xlim=c(0,10),ylim=c(0,10),xlab="Rep_1",ylab="Rep_2")
text(x=8,y=3,labels=corr)
title(args[2])
dev.off()
save.image(paste0(args[1],'corr.Rdata'))
