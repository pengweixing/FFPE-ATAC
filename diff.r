#################################################
#  File Name:diff.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sun 10 Jan 2021 02:30:45 PM UTC
#################################################

library(DESeq2)
library(edgeR)
args=commandArgs(T)
data = read.table(args[1],header = F)
name = read.table(args[2],header = F)
name2=c("chr","start","end",as.character(name$V1))
colnames(data)=name2

gene_length=data[,3]-data[,2] 
data2 = data[,4:ncol(data)]
data3=apply(data2,2,as.numeric)
data_norm=rpkm(data3,gene.length = gene_length)
#data_norm=cpm(data3)
#data_norm_1=log(data_norm+1,base=2)
write.table(data_norm,"data_all.rpkm.txt",sep="\t",quote=FALSE,row.names=FALSE)

data_norm2=as.data.frame(apply(data_norm,2,function(x) round(x)))

mynames=paste(rep("peak_",length(data2[,1])),row.names(data2),sep="")
row.names(data_norm2)=mynames
coldata = data.frame(row.names = colnames(data2),factor(as.character(name$V2)))
colnames(coldata)=c("type")
dds <- DESeqDataSetFromMatrix(data_norm2, coldata, design= ~ type)
dds <- DESeq(dds)
res= results(dds)
res = res[order(res$pvalue),]
diff_res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff = as.data.frame(diff_res)
row.names(data)=mynames
peak=data[,1:3]
all=rownames(diff)
all_peak=peak[match(all,rownames(data)),]
diff_all_peak=cbind(diff,all_peak)
write.table(diff_all_peak,"diff_peak_0.05.csv",sep="\t",quote=FALSE)

all_res = as.data.frame(res)
all_names = rownames(all_res)
all_peak=peak[match(all_names,rownames(data)),]
diff_all_peak2 = cbind(all_res,all_peak)
write.table(diff_all_peak2,"all_peak_diff.csv",sep="\t",quote=FALSE)
save.image('my.Rdata')
