#################################################
#  File Name:preseq.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Wed 07 Apr 2021 11:10:46 AM UTC
#################################################
### sh preseq.sh rep1.bam rep2.bam 

out=THS_Liv_50nm
/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq lc_extrap  -o $out.rep1.expect.txt -B $1 
/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq c_curve  -o $out.rep1.curve.txt -B $1 
/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq lc_extrap  -o $out.rep2.expect.txt -B $2 
/disk1/pengweixing/software/preseq/preseq-3.1.2/build/bin/preseq c_curve  -o $out.rep2.curve.txt -B $2 
#python /disk1/pengweixing/pipeline/library_complexity/preseq_plot.py rep1.curve.txt rep1.expect.txt rep2.curve.txt rep2.expect.txt NMK
