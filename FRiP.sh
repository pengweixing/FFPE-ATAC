#################################################
#  File Name:FRiP2.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 22 Jan 2021 06:37:04 PM UTC
#################################################

cat FRiP.list |while read line
do
bam=`echo $line|awk '{print $1}'`
bed=`echo $line|awk '{print $2}'`
out=`basename $bam|sed s/.sort.bam//g`
if [ ! -f $out.withinpeaks.bam ];then
samtools view -@ 10 -L $bed $bam -b -o $out.withinpeaks.bam
fi
if [ ! -f $bam.stat ];then
samtools flagstat $bam > $bam.stat
fi
if [ ! -f $out.withinpeaks.bam.stat ];then
samtools flagstat $out.withinpeaks.bam > $out.withinpeaks.bam.stat
fi
number1=`cat $bam.stat |sed -n '1p' |awk '{print $1}'`
number2=`cat $out.withinpeaks.bam.stat |sed -n '1p' |awk '{print $1}'`
out2=`echo "scale=4;$number2/$number1" |bc`
echo -e "$out\t$out2"
done
