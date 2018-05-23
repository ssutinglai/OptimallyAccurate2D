#!/usr/bin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot;
from pylab import genfromtxt; 
from scipy import integrate;
from scipy.integrate import odeint;

fig,ax = plt.subplots()

# Load the data
m1 = genfromtxt("../../../Katayama/Trace1.txt");
a1 = m1[:,1];
t1 = m1[:,0];
m2 = genfromtxt("../../../Katayama/Trace2.txt");
a2 = m2[:,1];
t2 = m2[:,0];
m3 = genfromtxt("../../../Katayama/Trace3.txt");
a3 = m3[:,1];
t3 = m3[:,0];
m4 = genfromtxt("../../../Katayama/Trace4.txt");
a4 = m4[:,1];
t4 = m4[:,0];
m5 = genfromtxt("../../../Katayama/Trace5.txt");
a5 = m5[:,1];
t5 = m5[:,0];
m6 = genfromtxt("../../../Katayama/Trace6.txt");
a6 = m6[:,1];
t6 = m6[:,0];
m7 = genfromtxt("../../../Katayama/Trace7.txt");
a7 = m7[:,1];
t7 = m7[:,0];
m8 = genfromtxt("../../../Katayama/Trace8.txt");
a8 = m8[:,1];
t8 = m8[:,0];
m9 = genfromtxt("../../../Katayama/Trace9.txt");
a9 = m9[:,1];
t9 = m9[:,0];
m10 = genfromtxt("../../../Katayama/Trace10.txt");
a10 = m10[:,1];
t10 = m10[:,0];
m11 = genfromtxt("../../../Katayama/Trace11.txt");
a11 = m11[:,1];
t11 = m11[:,0];
m12 = genfromtxt("../../../Katayama/Trace12.txt");
a12 = m12[:,1];
t12 = m12[:,0];
m13 = genfromtxt("../../../Katayama/Trace13.txt");
a13 = m13[:,1];
t13 = m13[:,0];
m14 = genfromtxt("../../../Katayama/Trace14.txt");
a14 = m14[:,1];
t14 = m14[:,0];
m15 = genfromtxt("../../../Katayama/Trace15.txt");
a15 = m15[:,1];
t15 = m15[:,0];
m16 = genfromtxt("../../../Katayama/Trace16.txt");
a16 = m16[:,1];
t16 = m16[:,0];
m17 = genfromtxt("../../../Katayama/Trace17.txt");
a17 = m17[:,1];
t17 = m17[:,0];
m18 = genfromtxt("../../../Katayama/Trace18.txt");
a18 = m18[:,1];
t18 = m18[:,0];
m19 = genfromtxt("../../../Katayama/Trace19.txt");
a19 = m19[:,1];
t19 = m19[:,0];
# m20 = genfromtxt("../../../Katayama/Trace1.txt");
# a20 = m20[:,1];
# t20 = m20[:,0];

#Load original data
mm1 = genfromtxt("../../../Katayama/iva1296/rtek00277.txt");
aa1 = mm1[:,1];
tt1 = mm1[:,0];
mm2 = genfromtxt("../../../Katayama/iva1296/rtek00279.txt");
aa2 = mm2[:,1];
tt2 = mm2[:,0];
mm3 = genfromtxt("../../../Katayama/iva1296/rtek00281.txt");
aa3 = mm3[:,1];
tt3 = mm3[:,0];
mm4 = genfromtxt("../../../Katayama/iva1296/rtek00283.txt");
aa4 = mm4[:,1];
tt4 = mm4[:,0];
mm5 = genfromtxt("../../../Katayama/iva1296/rtek00285.txt");
aa5 = mm5[:,1];
tt5 = mm5[:,0];
mm6 = genfromtxt("../../../Katayama/iva1296/rtek00287.txt");
aa6 = mm6[:,1];
tt6 = mm6[:,0];
mm7 = genfromtxt("../../../Katayama/iva1296/rtek00289.txt");
aa7 = mm7[:,1];
tt7 = mm7[:,0];
mm8 = genfromtxt("../../../Katayama/iva1296/rtek00291.txt");
aa8 = mm8[:,1];
tt8 = mm8[:,0];
mm9 = genfromtxt("../../../Katayama/iva1296/rtek00293.txt");
aa9 = mm9[:,1];
tt9 = mm9[:,0];
mm10 = genfromtxt("../../../Katayama/iva1296/rtek00295.txt");
aa10 = mm10[:,1];
tt10 = mm10[:,0];
mm11 = genfromtxt("../../../Katayama/iva1296/rtek00297.txt");
aa11 = mm11[:,1]
tt11 = mm11[:,0];
mm12 = genfromtxt("../../../Katayama/iva1296/rtek00299.txt");
aa12 = mm12[:,1];
tt12 = mm12[:,0];
mm13 = genfromtxt("../../../Katayama/iva1296/rtek00301.txt");
aa13 = mm13[:,1];
tt13 = mm13[:,0];
mm14 = genfromtxt("../../../Katayama/iva1296/rtek00303.txt");
aa14 = mm14[:,1];
tt14 = mm14[:,0];
mm15 = genfromtxt("../../../Katayama/iva1296/rtek00305.txt");
aa15 = mm15[:,1];
tt15 = mm15[:,0];
mm16 = genfromtxt("../../../Katayama/iva1296/rtek00307.txt");
aa16 = mm16[:,1];
tt16 = mm16[:,0];
mm17 = genfromtxt("../../../Katayama/iva1296/rtek00309.txt");
aa17 = mm17[:,1];
tt17 = mm17[:,0];
mm18 = genfromtxt("../../../Katayama/iva1296/rtek00311.txt");
aa18 = mm18[:,1];
tt18 = mm18[:,0];
mm19 = genfromtxt("../../../Katayama/iva1296/rtek00313.txt");
aa19 = mm19[:,1];
tt19 = mm19[:,0];
mm20 = genfromtxt("../../../Katayama/iva1296/rtek00315.txt");
aa20 = mm20[:,1];
tt20 = mm20[:,0];
mm21 = genfromtxt("../../../Katayama/iva1296/rtek00317.txt");
aa21 = mm21[:,1];
tt21 = mm21[:,0];

#Set offset and load parameters
div=50;

offset1=0;
x1=t1;
y1=offset1-a1/div;
offset2=0.04;
x2=t2;
y2=offset2-a2/div;
offset3=0.08;
x3=t3;
y3=offset3-a3/div;
offset4=0.12;
x4=t4;
y4=offset4-a4/div;
offset5=0.16;
x5=t5;
y5=offset5-a5/div;
offset6=0.20;
x6=t6;
y6=offset6-a6/div;
offset7=0.24;
x7=t7;
y7=offset7-a7/div;
offset8=0.28;
x8=t8;
y8=offset8-a8/div;
offset9=0.32;
x9=t9;
y9=offset9-a9/div;
offset10=0.36;
x10=t10;
y10=offset10-a10/div;
offset11=0.40;
x11=t11;
y11=offset11-a11/div;
offset12=0.44;
x12=t12;
y12=offset12-a12/div;
offset13=0.48;
x13=t13;
y13=offset13-a13/div;
offset14=0.52;
x14=t14;
y14=offset14-a14/div;
offset15=0.56;
x15=t15;
y15=offset15-a15/div;
offset16=0.60;
x16=t16;
y16=offset16-a16/div;
offset17=0.64;
x17=t17;
y17=offset17-a17/div;
offset18=0.68;
x18=t18;
y18=offset18-a18/div;
offset19=0.72;
x19=t19;
y19=offset19-a19/div;
# offset20=0.76;
# x20=t20;
# y20=offset20-a20/div;

#Set the parameters for the real observed data

ddiv=1;
offset1=0;
xx1=tt1;
yy1=offset1-aa1/ddiv;
offset2=0.04;
xx2=tt2;
yy2=offset2-aa2/ddiv;
offset3=0.08;
xx3=tt3;
yy3=offset3-aa3/ddiv;
offset4=0.12;
xx4=tt4;
yy4=offset4-aa4/ddiv;
offset5=0.16;
xx5=tt5;
yy5=offset5-aa5/ddiv;
offset6=0.20;
xx6=tt6;
yy6=offset6-aa6/ddiv;
offset7=0.24;
xx7=tt7;
yy7=offset7-aa7/ddiv;
offset8=0.28;
xx8=tt8;
yy8=offset8-aa8/ddiv;
offset9=0.32;
xx9=tt9;
yy9=offset9-aa9/ddiv;
offset10=0.36;
xx10=tt10;
yy10=offset10-aa10/ddiv;
offset11=0.40;
xx11=tt11;
yy11=offset11-aa11/ddiv;
offset12=0.44;
xx12=tt12;
yy12=offset12-aa12/ddiv;
offset13=0.48;
xx13=tt13;
yy13=offset13-aa13/ddiv;
offset14=0.52;
xx14=tt14;
yy14=offset14-aa14/ddiv;
offset15=0.56;
xx15=tt15;
yy15=offset15-aa15/ddiv;
offset16=0.60;
xx16=tt16;
yy16=offset16-aa16/ddiv;
offset17=0.64;
xx17=tt17;
yy17=offset17-aa17/ddiv;
offset18=0.68;
xx18=tt18;
yy18=offset18-aa18/ddiv;
offset19=0.72;
xx19=tt19;
yy19=offset19-aa19/ddiv;
offset20=0.76;
xx20=tt20;
yy20=offset20-aa20/ddiv;
offset21=0.80;
xx21=tt21;
yy21=offset21-aa21/ddiv;



#Plot
ax.plot(x1,y1,'k-',x2,y2,'k-',x3,y3,'k-',x4,y4,'k-',x5,y5,'k-',
	    x6,y6,'k-',x7,y7,'k-',x8,y8,'k-',x9,y9,'k-',x10,y10,'k-',
	    x11,y11,'k-',x12,y12,'k-',x13,y13,'k-',x14,y14,'k-',
	    x15,y15,'k-',x16,y16,'k-',x17,y17,'k-',x18,y18,'k-',
	    x19,y19,'k-',
	    # x20,y20,'k-',
	    xx1,yy1,'r:',xx2,yy2,'r:',xx3,yy3,'r:',xx4,yy4,'r:',
	    xx5,yy5,'r:',xx6,yy6,'r:',xx7,yy7,'r:',xx8,yy8,'r:',
	    xx9,yy9,'r:',xx10,yy10,'r:',xx11,yy11,'r:',xx12,yy12,'r:',
	    xx13,yy13,'r:',xx14,yy14,'r:',xx15,yy15,'r:',xx16,yy16,'r:',
	    xx17,yy17,'r:',xx18,yy18,'r:',xx19,yy19,'r:',xx20,yy20,'r:',
	    xx21,yy21,'r:')
# ax.grid(color='klinestyle=':', linewidth=1)
plt.axis([2*10**(-6), 12*10**(-6), -0.06, 0.86])
plt.xticks(np.arange(2e-6, 12e-6, 1e-6))
plt.yticks(np.arange(0, 1, 0.2))

ax = plt.gca() # grab the current axis for left side
ax.set_xticks([2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,9e-6,10e-6,11e-6,12e-6]) # choose which x locations to have ticks
ax.set_xticklabels([2,3,4,5,6,7,8,9,10,11,12]) 



plt.xlabel('Time(microsecond)')
plt.ylabel('Amplitude(V)')
plt.title('P-wave', x=0.94,y=0.99)
# plt.grid()
plt.show()