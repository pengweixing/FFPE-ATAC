#################################################
#  File Name:venn.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 21 May 2021 12:51:54 PM UTC
#################################################

FFPE=`cat FPE_ATAC_liver.merge.bed |wc -l`
THS_50nm=`cat THS_Liv_50nM.merge.bed |wc -l`
Omni=`cat Omni_ATAC_Liver.merge.bed|wc -l`
FFPE_vs_THS=`bedtools intersect -a FPE_ATAC_liver.merge.bed -b THS_Liv_50nM.merge.bed -wa |sort -u |wc -l`
FFPE_vs_Omni=`bedtools intersect -b FPE_ATAC_liver.merge.bed -a Omni_ATAC_Liver.merge.bed -wa |sort -u|wc -l`
THS_vs_Omni=`bedtools intersect -b THS_Liv_50nM.merge.bed -a Omni_ATAC_Liver.merge.bed -wa |sort -u |wc -l`
three=`bedtools intersect -a FPE_ATAC_liver.merge.bed -b THS_Liv_50nM.merge.bed -wa |sort -u |bedtools intersect -a stdin -b Omni_ATAC_Liver.merge.bed -wa |sort -u|wc -l`
echo -e "a1\ta2\ta3\ta1_a2\ta2_a3\ta1_a3\ta1_a2_a3"
echo -e "$FFPE\t$THS_50nm\t$Omni\t$FFPE_vs_THS\t$THS_vs_Omni\t$FFPE_vs_Omni\t$three"
Rscript venn.r Kidney_FFPE_THS50nM_Omni $FFPE $THS_50nm $Omni $FFPE_vs_THS $THS_vs_Omni $FFPE_vs_Omni $three FFPE THS_50nm Omni
