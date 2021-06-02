#################################################
#  File Name:FRiP.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 02 Apr 2021 09:47:02 AM UTC
#################################################
library(ggplot2)
data = read.table("FRiP.stat")
p=ggplot(data=data,aes(x=V1,y=V2))+geom_bar(aes(x=V1,y=V2),stat="identity",fill="#FF8100")+theme_bw()+theme(panel.grid = element_blank())+scale_y_continuous(expand = c(0,0))+xlab("")+ylab("FRiP")+theme(axis.text = element_text(color="black"))+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,color="black"))+geom_text(aes(label=round(V2,3)),position=position_stack(vjust=0.5))
ggsave("FRiP.pdf",plot=p,width=5,height=4)
save.image('my.Rdata')
