#################################################
#  File Name:relativetotss.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 23 Jan 2021 02:58:06 PM UTC
#################################################

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
args = commandArgs(T)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#files=list(T7="KM.T7.merge.bed")
files=list(FFPE1='FFPE_ATAC_liver_1.bed',FFPE2='FFPE_ATAC_liver_2.bed',THS1='THS_Liv_50nM_1.bed',THS2='THS_Liv_50nM_2.bed',SD1='SD_ATAC_Liver_1.bed',SD2='SD_ATAC_Liver_2.bed',Nomal1='NML_1.bed',Nomal2='NML_2.bed')
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
pdf(file=paste0("relative_to_tss.pdf"),height= 3, width = 6)
plotDistToTSS(peakAnnoList)
dev.off()

