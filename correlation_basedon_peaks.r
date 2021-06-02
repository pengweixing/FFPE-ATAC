#################################################
#  File Name:corr.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 14 Nov 2020 11:00:07 AM UTC
#################################################
library(edgeR)
args=commandArgs(T)
data = read.table(args[1])
mymean1 = apply(data[,4:5],1,mean)
#mymean1 = data[,4:5]
mymean2 = apply(data[,6:7],1,mean)
#mymean2 = data[,6]
new = data.frame(cuttag=mymean1,encode=mymean2)
new_cpm_temp = cpm(new)
new_cpm = log(new_cpm_temp+1,base=2)
corr <- paste("Pearson Cor. = ", round(cor(new_cpm[,2], new_cpm[,1]), 2), sep="")
pdf(paste0(args[1],"density.pdf"),height=6,width=6)
plot(new_cpm, col = densCols(new_cpm), pch = 20,xlim=c(0,10),ylim=c(0,10),xlab=args[2],ylab=args[3])
text(x=7,y=2,labels=corr)
dev.off()

pdf(paste0(args[1],".pdf"),height=6,width=6)
plot(new_cpm[,1],new_cpm[,2],pch=16,cex=0.7,xlim=c(0,10),ylim=c(0,10),xlab=args[2],ylab=args[3])
text(x=7,y=2,labels=corr)
title(args[2])
dev.off()
save.image(paste0(args[1],'corr.Rdata'))
