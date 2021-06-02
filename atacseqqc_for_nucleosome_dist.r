#################################################
#  File Name:atacseqqc.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 17 Apr 2021 11:15:21 AM UTC
#################################################

library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(BSgenome.Mmusculus.UCSC.mm9)
library(Rsamtools)
library(ChIPpeakAnno)

args = commandArgs(T)##args[1]:bam
txs <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene)
which <- as(seqinfo(Mmusculus), "GRanges")
gal <- readBamFile(args[1], which=which, asMates=TRUE, bigFile=TRUE)
genome <- Mmusculus
librarySize <- estLibSize(args[1])
cat(librarySize)
NTILE <- 101
dws <- ups <- 1010
gal1 = data_shift_mm <- shiftGAlignmentsList(gal)
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome)
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", "mononucleosome","dinucleosome","trinucleosome")], TSS=TSS, librarySize=librarySize, TSS.filter=0.5, n.tile = NTILE,upstream = ups, downstream = dws)
out <- featureAlignedDistribution(sigs,reCenterPeaks(TSS, width=ups+dws),zeroAt=.5, n.tile=NTILE, type="l",ylab="Averaged coverage")
write.table(out,file=paste0(args[1],'.txt'),sep="\t",quote=F)
pdf(args[2],height=6,width=8)
matplot(out, type="l", xaxt="n", xlab="Position (bp)", ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
writeListOfGAlignments(objs,outPath=args[3])
