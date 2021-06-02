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
library(MotifDb)
library(GenomicAlignments)

args = commandArgs(T)##args[1]:bam
txs <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene)
which <- as(seqinfo(Mmusculus), "GRanges")
genome <- Mmusculus
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2",
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"),
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
bamTop100 <- scanBam(BamFile(args[1], yieldSize = 100),param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
gal <- readBamFile(args[1], which=which, asMates=TRUE, bigFile=TRUE,tag=tags)
gal1 <- shiftGAlignmentsList(gal)
galist=GAlignmentsList(shifted=gal1)
writeListOfGAlignments(galist,outPath=args[3])
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
pdf(paste0(args[3],"/",args[2],"footprint.pdf"),width=6,height=6)
sigs <- factorFootprints(paste0(args[2],"/shifted.bam"), pfm=CTCF[[1]], 
                         genome=genome, ## Don't have a genome? ask ?factorFootprints for help
                         min.score="90%", seqlev = paste0("chr", c(1:19, "X", "Y")),
                         upstream=100, downstream=100)
dev.off()
seqlev= paste0("chr", c(1:19, "X", "Y"))
pdf(paste0(args[3],"/",args[2],"point.pdf"),width=6,height=6)
vp <- vPlot(shiftedBamfile, pfm=CTCF[[1]],
            genome=genome, min.score="90%", seqlev=seqlev,
            upstream=200, downstream=200,
            ylim=c(30, 250), bandwidth=c(2, 1))
dev.off()
