#################################################
#  File Name:SDS_anno2.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Mon 07 Dec 2020 06:03:34 PM UTC
#################################################
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)

FFPE1 = read.table("FFPE_ATAC_liver_1.bed")
colnames(FFPE1) = c("chr","start","end")
FFPE2 = read.table("FFPE_ATAC_liver_2.bed")
colnames(FFPE2) = c("chr","start","end")

THS1 = read.table("THS_Liv_50nM_1.bed")
colnames(THS1) = c("chr","start","end")
THS2 = read.table("THS_Liv_50nM_2.bed")
colnames(THS2) = c("chr","start","end")

SD1 = read.table("SD_ATAC_Liver_1.bed")
colnames(SD1) = c("chr","start","end")
SD2 = read.table("SD_ATAC_Liver_2.bed")
colnames(SD2) = c("chr","start","end")

Normal1 = read.table("NML_1.bed")
colnames(Normal1) = c("chr","start","end")
Normal2 = read.table("NML_2.bed")
colnames(Normal2) = c("chr","start","end")

bbb <- function(data)
{
SDS_data = data
SDS_data=makeGRangesFromDataFrame(SDS_data)
SDS_data_anno = annotatePeak(SDS_data,tssRegion = c(-3000,3000),TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene)
SDS_anno = SDS_data_anno@annoStat
SDS_anno2 = data.frame(Feature=c("Promoter","5'UTR","3'UTR","Exon","Intron","Dwonstream","Distal Intergenic"),Frequency=c(sum(SDS_anno[1:3,2]),SDS_anno[4,2],SDS_anno[5,2],sum(SDS_anno[6:7,2]),sum(SDS_anno[8:9,2]),SDS_anno[10,2],SDS_anno[11,2]))
SDS_anno2$Feature=factor(SDS_anno2$Feature,levels = c("Promoter","Intron","Distal Intergenic","Dwonstream", "Exon","3'UTR" ,"5'UTR"))
SDS_data_anno@annoStat=SDS_anno2
return(SDS_data_anno)
}

FFPE1_anno = bbb(FFPE1)
FFPE2_anno = bbb(FFPE2)

THS1_anno = bbb(THS1)
THS2_anno = bbb(THS2)

SD1_anno = bbb(SD1)
SD2_anno = bbb(SD2)

Normal1_anno = bbb(Normal1)
Normal2_anno = bbb(Normal2)

all = list(FFPE1 =FFPE1_anno,FFPE2 =FFPE2_anno, THS1 = THS1_anno, THS2 = THS2_anno, SD1 = SD1_anno, SD2 = SD2_anno,Normal1 = Normal1_anno, Normal2 = Normal2_anno)

pdf('genomic_Peakanno_all.pdf',width=8,height=6)
plotAnnoBar(all)
dev.off()

