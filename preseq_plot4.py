#################################################
#  File Name:preseq_plot.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Wed 07 Apr 2021 04:47:18 PM UTC
#################################################

import numpy as np
from matplotlib import pyplot as plt
import sys
from matplotlib.lines import Line2D
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42

Fobserved1 = np.loadtxt('trim_liver_1.sort.rep1.curve.txt', skiprows = 1 )
Fpredicted1 = np.loadtxt('trim_liver_1.sort.rep1.expect.txt',skiprows = 1)
Fobserved2 = np.loadtxt('trim_liver_1.sort.rep2.curve.txt',skiprows = 1)
Fpredicted2 = np.loadtxt('trim_liver_1.sort.rep2.expect.txt',skiprows = 1)

Oobserved1 = np.loadtxt('Omni_ATAC_Liv_rep1.curve.txt',skiprows = 1)
Opredicted1 = np.loadtxt('Omni_ATAC_Liv_rep1.expect.txt',skiprows = 1)
Oobserved2 = np.loadtxt('Omni_ATAC_Liv_rep2.curve.txt',skiprows = 1)
Opredicted2 = np.loadtxt('Omni_ATAC_Liv_rep2.expect.txt',skiprows = 1)

THS_observed1 = np.loadtxt('THS_Liv_50nm.rep1.curve.txt', skiprows = 1)
THS_predicted1 = np.loadtxt('THS_Liv_50nm.rep1.expect.txt', skiprows = 1)
THS_observed2 = np.loadtxt('THS_Liv_50nm.rep2.curve.txt', skiprows = 1)
THS_predicted2 = np.loadtxt('THS_Liv_50nm.rep2.expect.txt', skiprows = 1)

Nobserved1 = np.loadtxt('NML.rep1.curve.txt', skiprows = 1)
Npredicted1 = np.loadtxt('NML.rep1.expect.txt', skiprows = 1)
Nobserved2 = np.loadtxt('NML.rep2.curve.txt', skiprows = 1)
Npredicted2 = np.loadtxt('NML.rep2.expect.txt', skiprows = 1)

Sobserved1 = np.loadtxt('SD_ATAC_rep1.curve.txt', skiprows = 1)
Spredicted1 = np.loadtxt('SD_ATAC_rep1.expect.txt', skiprows = 1)
Sobserved2 = np.loadtxt('SD_ATAC_rep2.curve.txt', skiprows = 1)
Spredicted2 = np.loadtxt('SD_ATAC_rep2.expect.txt', skiprows = 1)

fig = plt.figure()
out=open(sys.argv[1]+".txt",'w+')
myname=["FFPE_rep1_x","FFPE_rep1_y","FFPE_rep2_x","FFPE_rep2_y","Omni_rep1_x","Omni_rep1_y","Omni_rep2_x","Omni_rep2_y","THS_rep1_x","THS_rep1_y","THS_rep2_x","THS_rep2_y","Normal_rep1_x","Normal_rep1_y","Normal_rep2_x","Normal_rep2_y","SD_ATAC_1_x","SD_ATAC_1_y","SD_ATAC_2_x","SD_ATAC_2_y"]
print(file=out,*myname,sep="\t")
fp1 = np.round(np.log10(Fpredicted1[1:,0:2]),2)
fp2 = np.round(np.log10(Fpredicted2[1:,0:2]),2)
op1 = np.round(np.log10(Opredicted1[1:,0:2]),2)
op2 = np.round(np.log10(Opredicted2[1:,0:2]),2)
tp1 = np.round(np.log10(THS_predicted1[1:,0:2]),2)
tp2 = np.round(np.log10(THS_predicted2[1:,0:2]),2)
np1 = np.round(np.log10(Npredicted1[1:,0:2]),2)
np2 = np.round(np.log10(Npredicted2[1:,0:2]),2)
sp1 = np.round(np.log10(Spredicted1[1:,0:2]),2)
sp2 = np.round(np.log10(Spredicted2[1:,0:2]),2)
all_list = np.hstack((fp1,fp2,op1,op2,tp1,tp2,np1,np2,sp1,sp2))
np.savetxt(out,all_list,delimiter="\t",fmt="%.2f")

### FFPE###
#rep1
F1, = plt.plot(np.log10(Fobserved1[:,0]),np.log10(Fobserved1[:,1]),'r-')
plt.plot(np.log10(Fpredicted1[:,0]),np.log10(Fpredicted1[:,1]),'r--')

#rep2
F2, = plt.plot(np.log10(Fobserved2[:,0]),np.log10(Fobserved2[:,1]),ls='-',color='orangered')
plt.plot(np.log10(Fpredicted2[:,0]),np.log10(Fpredicted2[:,1]),ls='--',color='orangered')

### Omni###
#rep1
O1, = plt.plot(np.log10(Oobserved1[:,0]),np.log10(Oobserved1[:,1]),'b-')
plt.plot(np.log10(Opredicted1[:,0]),np.log10(Opredicted1[:,1]),'b--')


#rep2
O2, = plt.plot(np.log10(Oobserved2[:,0]),np.log10(Oobserved2[:,1]),ls='-',color='mediumblue')
plt.plot(np.log10(Opredicted2[:,0]),np.log10(Opredicted2[:,1]),ls='--',color='mediumblue')

### Omni###
#rep1
O1_f, = plt.plot(np.log10(THS_observed1[:,0]),np.log10(THS_observed1[:,1]),ls='-',color='violet')
plt.plot(np.log10(THS_predicted1[:,0]),np.log10(THS_predicted1[:,1]),ls='--',color='violet')
#rep2
O2_f, = plt.plot(np.log10(THS_observed2[:,0]),np.log10(THS_observed2[:,1]),ls='-',color='blueviolet')
plt.plot(np.log10(THS_predicted2[:,0]),np.log10(THS_predicted2[:,1]),ls='--',color='blueviolet')

### Normal###
#rep1
N1, = plt.plot(np.log10(Nobserved1[:,0]),np.log10(Nobserved1[:,1]),'g-')
plt.plot(np.log10(Npredicted1[:,0]),np.log10(Npredicted1[:,1]),'g--')
#rep2
N2, = plt.plot(np.log10(Nobserved2[:,0]),np.log10(Nobserved2[:,1]),ls='-',color='lime')
plt.plot(np.log10(Npredicted2[:,0]),np.log10(Npredicted2[:,1]),ls='--',color='lime')

### Standard atac###
#rep1
SN1, = plt.plot(np.log10(Sobserved1[:,0]),np.log10(Sobserved1[:,1]),'-',color='olive')
plt.plot(np.log10(Spredicted1[:,0]),np.log10(Spredicted1[:,1]),ls='--',color='olive')
#rep2
SN2, = plt.plot(np.log10(Sobserved2[:,0]),np.log10(Sobserved2[:,1]),ls='-',color='olivedrab')
plt.plot(np.log10(Spredicted2[:,0]),np.log10(Spredicted2[:,1]),ls='--',color='olivedrab')

#####
plt.title('Preseq estimated yield')
plt.xlabel('Total mapped fragments (log10)')
plt.ylabel('Distinct fragments (log10)')

myline = [F1,F2,O1,O2,O1_f,O2_f,N1,N2,SN1,SN2,Line2D([0],[0],color="black",lw=2,ls="-"),Line2D([0],[10],color="black",lw=2,ls="--")]

plt.legend(myline,["FFPE_rep1","FFPE_rep2","Omni_rep1","Omni_rep2","THS_rep1","THS_rep2","Normal_rep1","Normal_rep2","SD_ATAC1","SD_ATAC2","Observed","Predicted"],handlelength=3,bbox_to_anchor=(1.45, 1),loc='upper right')

plot_img = sys.argv[1]+'.pdf'
fig.savefig(plot_img,format='pdf',bbox_inches='tight')
