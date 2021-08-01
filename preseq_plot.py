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

Fobserved1 = np.loadtxt('GM-T7-ac-1_R1.curve.txt', skiprows = 1 )
Fpredicted1 = np.loadtxt('GM-T7-ac-1_R1.expect.txt',skiprows = 1)
Fobserved2 = np.loadtxt('GM-T7-ac-2_R1.curve.txt',skiprows = 1)
Fpredicted2 = np.loadtxt('GM-T7-ac-2_R1.expect.txt',skiprows = 1)

Oobserved1 = np.loadtxt('GM-normal-ac-1_R1.curve.txt',skiprows = 1)
Opredicted1 = np.loadtxt('GM-normal-ac-1_R1.expect.txt',skiprows = 1)
Oobserved2 = np.loadtxt('GM-normal-ac-2_R1.curve.txt',skiprows = 1)
Opredicted2 = np.loadtxt('GM-normal-ac-2_R1.expect.txt',skiprows = 1)

THS_observed1 = np.loadtxt('GM_H3K27ac_ENCODE_rep1.curve.txt', skiprows = 1)
THS_predicted1 = np.loadtxt('GM_H3K27ac_ENCODE_rep1.expect.txt', skiprows = 1)
THS_observed2 = np.loadtxt('GM_H3K27ac_ENCODE_rep2.curve.txt', skiprows = 1)
THS_predicted2 = np.loadtxt('GM_H3K27ac_ENCODE_rep2.expect.txt', skiprows = 1)

fig = plt.figure()
out=open(sys.argv[1]+".txt",'w+')
myname=["FFPE_rep1_x","FFPE_rep1_y","FFPE_rep2_x","FFPE_rep2_y","Omni_rep1_x","Omni_rep1_y","Omni_rep2_x","Omni_rep2_y","THS_rep1_x","THS_rep1_y","THS_rep2_x","THS_rep2_y"]
print(file=out,*myname,sep="\t")
fp1 = np.round(np.log10(Fpredicted1[1:,0:2]),2)
fp2 = np.round(np.log10(Fpredicted2[1:,0:2]),2)
op1 = np.round(np.log10(Opredicted1[1:,0:2]),2)
op2 = np.round(np.log10(Opredicted2[1:,0:2]),2)
tp1 = np.round(np.log10(THS_predicted1[1:,0:2]),2)
tp2 = np.round(np.log10(THS_predicted2[1:,0:2]),2)
all_list = np.hstack((fp1,fp2,op1,op2,tp1,tp2))
np.savetxt(out,all_list,delimiter="\t",fmt="%.2f")

### FFPE###
#rep1
F1, = plt.plot(np.log10(Fobserved1[:,0]),np.log10(Fobserved1[:,1]),'r-')
plt.plot(np.log10(Fpredicted1[:,0]),np.log10(Fpredicted1[:,1]),'r--')
#rep2
F2, = plt.plot(np.log10(Fobserved2[:,0]),np.log10(Fobserved2[:,1]),ls='-',color='darkorange')
plt.plot(np.log10(Fpredicted2[:,0]),np.log10(Fpredicted2[:,1]),ls='--',color='darkorange')

### Omni###
#rep1
O1, = plt.plot(np.log10(Oobserved1[:,0]),np.log10(Oobserved1[:,1]),'b-')
plt.plot(np.log10(Opredicted1[:,0]),np.log10(Opredicted1[:,1]),'b--')
#rep2
O2, = plt.plot(np.log10(Oobserved2[:,0]),np.log10(Oobserved2[:,1]),ls='-',color='mediumblue')
plt.plot(np.log10(Opredicted2[:,0]),np.log10(Opredicted2[:,1]),ls='--',color='mediumblue')

### THS###
#rep1
O1_f, = plt.plot(np.log10(THS_observed1[:,0]),np.log10(THS_observed1[:,1]),ls='-',color='g')
plt.plot(np.log10(THS_predicted1[:,0]),np.log10(THS_predicted1[:,1]),ls='--',color='g')
#rep2
O2_f, = plt.plot(np.log10(THS_observed2[:,0]),np.log10(THS_observed2[:,1]),ls='-',color='lime')
plt.plot(np.log10(THS_predicted2[:,0]),np.log10(THS_predicted2[:,1]),ls='--',color='lime')

#####
plt.title('Preseq estimated yield')
plt.xlabel('Total mapped fragments (log10)')
plt.ylabel('Distinct fragments (log10)')

myline = [F1,F2,O1,O2,O1_f,O2_f,Line2D([0],[0],color="black",lw=2,ls="-"),Line2D([0],[10],color="black",lw=2,ls="--")]

plt.legend(myline,["FACT_rep1","FACT_rep2","Normal_rep1","Normal_rep2","ENCODE_rep1","ENCODE_rep2","Observed","Predicted"],handlelength=3,bbox_to_anchor=(1.45, 1),loc='upper right')

plot_img = sys.argv[1]+'.pdf'
fig.savefig(plot_img,format='pdf',bbox_inches='tight')
