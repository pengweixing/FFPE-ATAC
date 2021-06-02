#################################################
#  File Name:relativetotss.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 23 Jan 2021 02:58:06 PM UTC
#################################################
library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
args = commandArgs(T)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#files=list(T7="KM.T7.merge.bed")
files=list(Overlap=args[1],FFPE_exclusive=args[2],THS_exclusive=args[3])
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
overlap_dist = peakAnnoList[[1]]@anno@elementMetadata@listData[["distanceToTSS"]]
FFPE_exclusive_dist = peakAnnoList[[2]]@anno@elementMetadata@listData[["distanceToTSS"]]
THS_exclusive_dist = peakAnnoList[[3]]@anno@elementMetadata@listData[["distanceToTSS"]]
limit <- c(-100000,-10000,-5000,0, 5000, 10000, 100000)
lbs <- c("<=-100kb","(-100kb,-10kb]","(-10kb,5kb]","(-5kb,0]", "[0,5kb)", "[5kb,10kb)", "(10kb,100kb)", ">=100kb")
overlap_stat = data.frame(range=lbs,value=NA)
FFPE_exclusive_stat = data.frame(range=lbs,value=NA)
THS_exclusive_stat = data.frame(range=lbs,value=NA)
overlap_stat$range=factor(overlap_stat$range,levels =c("<=-100kb","(-100kb,-10kb]","(-10kb,5kb]","(-5kb,0]", "[0,5kb)", "[5kb,10kb)", "(10kb,100kb)", ">=100kb"))
FFPE_exclusive_stat$range=factor(FFPE_exclusive_stat$range,levels =c("<=-100kb","(-100kb,-10kb]","(-10kb,5kb]","(-5kb,0]", "[0,5kb)", "[5kb,10kb)", "(10kb,100kb)",">=100kb"))
THS_exclusive_stat$range=factor(THS_exclusive_stat$range,levels =c("<=-100kb","(-100kb,-10kb]","(-10kb,5kb]","(-5kb,0]", "[0,5kb)", "[5kb,10kb)", "(10kb,100kb)",">=100kb"))

zero_overlap = length(which(overlap_dist==0))/2
zero_FFPE = length(which(FFPE_exclusive_dist==0))/2
zero_THS = length(which(FFPE_exclusive_dist==0))/2
for (i in 1:8){
    if(i == 1){
        overlap_stat[i,2]=length(which(overlap_dist<=limit[1]))
    }
    else if (i >1 & i<8){
    overlap_stat[i,2]=length(which(overlap_dist>limit[i-1] & overlap_dist< limit[i]))
    }
    else{
    overlap_stat[8,2]=length(which(overlap_dist>=limit[7]))
    }
}
overlap_stat[4,2]=overlap_stat[4,2]+zero_overlap
overlap_stat[5,2]=overlap_stat[5,2]+zero_overlap
overlap_stat[2,2]=overlap_stat[2,2]+length(which(overlap_dist==-10000))
overlap_stat[3,2]=overlap_stat[3,2]+length(which(overlap_dist==-5000))
overlap_stat[5,2]=overlap_stat[5,2]+length(which(overlap_dist==-5000))
overlap_stat[6,2]=overlap_stat[6,2]+length(which(overlap_dist==-10000))

for (i in 1:8){
    if(i == 1){
        FFPE_exclusive_stat[i,2]=length(which(FFPE_exclusive_dist<limit[1]))
    }
    else if (i >1 & i<8){
    FFPE_exclusive_stat[i,2]=length(which(FFPE_exclusive_dist>limit[i-1] & FFPE_exclusive_dist< limit[i]))
    }
    else{
    FFPE_exclusive_stat[8,2]=length(which(FFPE_exclusive_dist>limit[7]))
    }
}
FFPE_exclusive_stat[4,2] = FFPE_exclusive_stat[4,2]+zero_FFPE
FFPE_exclusive_stat[5,2] = FFPE_exclusive_stat[5,2]+zero_FFPE
FFPE_exclusive_stat[2,2]=FFPE_exclusive_stat[2,2]+length(which(FFPE_exclusive_dist==-10000))
FFPE_exclusive_stat[3,2]=FFPE_exclusive_stat[3,2]+length(which(FFPE_exclusive_dist==-5000))
FFPE_exclusive_stat[5,2]=FFPE_exclusive_stat[5,2]+length(which(FFPE_exclusive_dist==-5000))
FFPE_exclusive_stat[6,2]=FFPE_exclusive_stat[6,2]+length(which(FFPE_exclusive_dist==-10000))
for (i in 1:8){
    if(i == 1){
        THS_exclusive_stat[i,2]=length(which(THS_exclusive_dist<limit[1]))
    }
    else if (i >1 & i<8){
    THS_exclusive_stat[i,2]=length(which(THS_exclusive_dist>limit[i-1] & THS_exclusive_dist< limit[i]))
    }
    else{
    THS_exclusive_stat[8,2]=length(which(THS_exclusive_dist>limit[7]))
    }
}
THS_exclusive_stat[4,2] = THS_exclusive_stat[4,2]+zero_THS
THS_exclusive_stat[5,2] = THS_exclusive_stat[5,2]+zero_THS
THS_exclusive_stat[2,2]=THS_exclusive_stat[2,2]+length(which(THS_exclusive_dist==-10000))
THS_exclusive_stat[3,2]=THS_exclusive_stat[3,2]+length(which(THS_exclusive_dist==-5000))
THS_exclusive_stat[5,2]=THS_exclusive_stat[5,2]+length(which(THS_exclusive_dist==-5000))
THS_exclusive_stat[6,2]=THS_exclusive_stat[6,2]+length(which(THS_exclusive_dist==-10000))

overlap_stat$value=overlap_stat$value/sum(overlap_stat$value)*100
FFPE_exclusive_stat$value=FFPE_exclusive_stat$value/sum(FFPE_exclusive_stat$value)*100
THS_exclusive_stat$value=THS_exclusive_stat$value/sum(THS_exclusive_stat$value)*100

p_overlap = ggplot(data=overlap_stat,aes(x=range,y=value))+geom_bar(aes(x=range,y=value),stat="identity",fill="#015666")+theme_bw()+theme(axis.text.x = element_text(angle=45,vjust=0.6,color = "black"))+theme(panel.grid.major = element_blank())+ylab("The percentage of peaks (%)")+xlab("")+geom_text(aes(label=paste0(round(value,2),"%"),y = value + 2))
p_ffpe = ggplot(data=FFPE_exclusive_stat,aes(x=range,y=value))+geom_bar(aes(x=range,y=value),stat="identity",fill="#015666")+theme_bw()+theme(axis.text.x = element_text(angle=45,vjust=0.6,color = "black"))+theme(panel.grid.major = element_blank())+ylab("The percentage of peaks (%)")+xlab("")+geom_text(aes(label=paste0(round(value,2),"%"),y = value + 2))
p_ths = ggplot(data=THS_exclusive_stat,aes(x=range,y=value))+geom_bar(aes(x=range,y=value),stat="identity",fill="#015666")+theme_bw()+theme(axis.text.x = element_text(angle=45,vjust=0.6,color = "black"))+theme(panel.grid.major = element_blank())+ylab("The percentage of peaks (%)")+xlab("")+geom_text(aes(label=paste0(round(value,2),"%"),y = value + 2))
p_all=ggarrange(p_overlap,p_ffpe,p_ths,ncol=1,nrow=3,labels=c(args[1],args[2],args[3]))
ggsave(paste0(args[4],".pdf"),plot=p_all,height = 7,width = 6)

#pdf(file=paste0("relative_to_tss.pdf"),height= 3, width = 6)
#plotDistToTSS(peakAnnoList)
#dev.off()

save.image("my.Rdata")
