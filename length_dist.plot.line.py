#################################################
#  File Name:length_dist.plot.single.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Thu 06 May 2021 10:55:56 AM UTC
#################################################
import numpy as np
from matplotlib import pyplot as plt

from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

NMK1 = np.loadtxt("NMK_1_L001.frag.sort.txt")
NMK2 = np.loadtxt("NMK_2_L001.frag.sort.txt")
Omni1 = np.loadtxt("Kidney_omni_1_total.frag.sort.txt")
Omni2 = np.loadtxt("Kidney_omni_2.frag.sort.txt")

fig = plt.figure(figsize=(12,12))
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)
N1, = ax1.plot(NMK1[:, 1], NMK1[:, 0], color='r')
N2, = ax2.plot(NMK2[:, 1], NMK2[:, 0], color='darkorange')
O1, = ax3.plot(Omni1[:, 1], Omni1[:, 0], color='b')
O2, = ax4.plot(Omni2[:, 1], Omni2[:, 0], color='navy')
ax1.set_xlabel("NMK_1")
ax2.set_xlabel("NMK_2")
ax3.set_xlabel("Omni_1")
ax4.set_xlabel("Omni_1")
ax1.set_ylabel("Read Counts")
ax2.set_ylabel("Read Counts")
ax3.set_ylabel("Read Counts")
ax4.set_ylabel("Read Counts")
fig.savefig("length_distribution.independent.pdf", format='pdf', bbox_inches='tight')
