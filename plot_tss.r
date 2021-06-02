#################################################
#  File Name:plot_tss.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 01 May 2021 10:04:34 PM UTC
#################################################

library(ggplot2)
args = commandArgs(T)
data = read.table(args[1])
a = apply(data,2,sum)
tss = data.frame(score=a)
tss2 = tss/mean(tss[1:3,])
rownames(tss2)=seq(-2000,1980,by=20)
p=ggplot(data=tss2)+geom_line(aes(x=seq(1,200),y=score))+theme_bw()+xlab("Distance from TSS (bp)")+scale_x_continuous(breaks=c(0,50,100,150,200),labels=c(-2000,-1000,0,1000,2000))
ggsave(paste0(args[1],".pdf"),plot=p,width=5,height=4)
