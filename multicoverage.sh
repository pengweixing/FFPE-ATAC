#################################################
#  File Name:cov.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 30 Apr 2021 11:37:33 AM UTC
#################################################

#bedtools multicov -bams  500c_rep1.pe.sort.rmdup.bam 500c_rep2.pe.sort.rmdup.bam -bed all.bed > all.bed.cov 
cat *narrowPeak |sort -k1,1 -k2,2n |bedtools merge -i stdin > all.bed
bedtools multicov -bams ../trim_kidney_1.bam.sort.rmdup.bam ../trim_kidney_2.bam.sort.rmdup.bam  500c_rep1.pe.sort.rmdup.0.8down.bam 500c_rep2.pe.sort.rmdup.0.8down.bam -bed all.bed > all.bed.cov
Rscript corr.r all.bed.cov FFPE 500c
