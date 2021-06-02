#################################################
#  File Name:relativetotss_barplot_percent2.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Mon 31 May 2021 01:02:36 PM UTC
#################################################

library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
args = commandArgs(T)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
files=list(Overlap=args[1])
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
Differential_dist = peakAnnoList[[1]]@anno@elementMetadata@listData[["distanceToTSS"]]
limit <- c(-100000,-10000,-5000,0, 5000, 10000, 100000)
lbs <- c("<=-100kb","(-100kb,-10kb]","(-10kb,5kb]","(-5kb,0]", "[0,5kb)", "[5kb,10kb)", "(10kb,100kb)", ">=100kb")
Differential_stat = data.frame(range=lbs,value=NA)
Differential_stat$range=factor(Differential_stat$range,levels =c("<=-100kb","(-100kb,-10kb]","(-10kb,5kb]","(-5kb,0]", "[0,5kb)", "[5kb,10kb)", "(10kb,100kb)", ">=100kb"))
diff_overlap = length(which(Differential_dist==0))/2
for (i in 1:8){
    if(i == 1){
        Differential_stat[i,2]=length(which(Differential_dist<=limit[1]))
    }
    else if (i >1 & i<8){
    Differential_stat[i,2]=length(which(Differential_dist>limit[i-1] & Differential_dist< limit[i]))
    }
    else{
    Differential_stat[8,2]=length(which(Differential_dist>=limit[7]))
    }
}
Differential_stat[4,2]=Differential_stat[4,2]+diff_overlap
Differential_stat[5,2]=Differential_stat[5,2]+diff_overlap
Differential_stat[2,2]=Differential_stat[2,2]+length(which(Differential_dist==-10000))
Differential_stat[3,2]=Differential_stat[3,2]+length(which(Differential_dist==-5000))
p_Differential = ggplot(data=Differential_stat,aes(x=range,y=value))+geom_bar(aes(x=range,y=value),stat="identity",fill="#015666")+theme_bw()+theme(axis.text.x = element_text(angle=45,vjust=0.6,color = "black"))+theme(panel.grid.major = element_blank())+ylab("The percentage of peaks (%)")+xlab("")+geom_text(aes(label=value,vjust=-0.3))
ggsave(args[2],plot=p_Differential,height = 4,width = 6)
save.image('my.Rdata')

