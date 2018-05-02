#!/usr/bin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot;
from pylab import genfromtxt; 
from scipy import integrate;
from scipy.integrate import odeint;

fig,ax = plt.subplots()

# Load the data
m1 = genfromtxt("100418_1_Con/Receiver01P.txt");
a1 = m1[:,1];
t1 = m1[:,0];
m2 = genfromtxt("110418_3_Con(2)/Receiver01PP.txt");
a2 = m2[:,1];
t2 = m2[:,0];
m3 = genfromtxt("Con(2.5)/Receiver01PPP.txt");
a3 = m3[:,1];
t3 = m3[:,0];
m4 = genfromtxt("110418_4_Con(3)/Receiver01PPP.txt");
a4 = m4[:,1];
t4 = m4[:,0];
m5 = genfromtxt("Con(3.5.2)/Receiver01PPP.txt");
a5 = m5[:,1];
t5 = m5[:,0];
m6 = genfromtxt("Con(3.5.1)/Receiver01PPP.txt");
a6 = m6[:,1];
t6 = m6[:,0];
m7 = genfromtxt("Con(3.5.3)/Receiver01PPP.txt");
a7 = m7[:,1];
t7 = m7[:,0];
m8 = genfromtxt("Con(3.5.3)/Receiver01PPP.txt");
a8 = m8[:,1];
t8 = m8[:,0];
m9 = genfromtxt("Con(3.5.1)/Receiver01PPP.txt");
a9 = m9[:,1];
t9 = m9[:,0];
m10 = genfromtxt("110418_4_Con(3)/Receiver01PPP.txt");
a10 = m10[:,1];
t10 = m10[:,0];
m11 = genfromtxt("Con(1.5)/Receiver01PPP.txt");
a11 = m11[:,1];
t11 = m11[:,0];
m12 = genfromtxt("Con(23)/Receiver01PPP.txt");
a12 = m12[:,1];
t12 = m12[:,0];
m13 = genfromtxt("Con(28)/Receiver01PPP.txt");
a13 = m13[:,1];
t13 = m13[:,0];
m14 = genfromtxt("Con(32)/Receiver01PPP.txt");
a14 = m14[:,1];
t14 = m14[:,0];
#m15 = genfromtxt("Con(15)/Receiver01PPP.txt");
#a15 = m15[:,1];
#t15 = m15[:,0];
#m16 = genfromtxt("Con(16)/Receiver01PPP.txt");
#a16 = m16[:,1];
#t16 = m16[:,0];
#m17 = genfromtxt("Con(17)/Receiver01PPP.txt");
#a17 = m17[:,1];
#t17 = m17[:,0];
#m18 = genfromtxt("Con(18)/Receiver01PPP.txt");
#a18 = m18[:,1];
#t18 = m18[:,0];
#m19 = genfromtxt("Con(19)/Receiver01PPP.txt");
#a19 = m19[:,1];
#t19 = m19[:,0];
#m20 = genfromtxt("Con(20)/Receiver01PPP.txt");
#a20 = m20[:,1];
#t20 = m20[:,0];
#m21 = genfromtxt("Con(21)/Receiver01PPP.txt");
#a21 = m21[:,1];
#t21 = m21[:,0];


#Set offset and load parameters
div=65;

offset1=0.11;
x1=t1;
y1=offset1-a1/div;
offset2=0.13;
x2=t2;
y2=offset2-a2/div;
offset3=0.15;
x3=t3;
y3=offset3-a3/div;
offset4=0.17;
x4=t4;
y4=offset4-a4/div;
offset5=0.19;
x5=t5;
y5=offset5-a5/div;
offset6=0.21;
x6=t6;
y6=offset6-a6/div;
offset7=0.23;
x7=t7;
y7=offset7-a7/div;
offset8=0.25;
x8=t8;
y8=offset8-a8/div;
offset9=0.27;
x9=t9;
y9=offset9-a9/div;
offset10=0.29;
x10=t10;
y10=offset10-a10/div;
offset11=0.31;
x11=t11;
y11=offset11-a11/div;
offset12=0.33;
x12=t12;
y12=offset12-a12/div;
offset13=0.35;
x13=t13;
y13=offset13-a13/div;
offset14=0.37;
x14=t14;
y14=offset14-a14/div;
#offset15=0.39;
#x15=t15;
#y15=offset15-a15/div;
#offset16=0.41;
#x16=t16;
#y16=offset16-a16/div;
#offset17=0.43;
#x17=t17;
#y17=offset17-a17/div;
#offset18=0.45;
#x18=t18;
#y18=offset18-a18/div;
#offset19=0.47;
#x19=t19;
#y19=offset19-a19/div;
#offset20=0.49;
#x20=t20;
#y20=offset20-a20/div;
#offset21=0.51;
#x21=t21;
#y21=offset21-a21/div;


#Plot
ax.plot(x1,y1,'k-',x2,y2,'k-',x3,y3,'k-',x4,y4,'k-',x5,y5,'k-',
	    x6,y6,'k-',x7,y7,'k-',x8,y8,'k-',x9,y9,'k-',x10,y10,'k-',
	    x11,y11,'k-',x12,y12,'k-',x13,y13,'k-',x14,y14,'k-')
#	    x15,y15,'k-',x16,y16,'k-',x17,y17,'k-',x18,y18,'k-',
#	    x19,y19,'k-',x20,y20,'k-',x21,y21,'k-')
# ax.grid(color='k', linestyle=':', linewidth=1)
plt.axis([2*10**(-6), 12*10**(-6), 0.06, 0.54])
plt.xticks(np.arange(2e-6, 12e-6, 1e-6))
plt.yticks(np.arange(0.1, 0.4, 0.1))

ax = plt.gca() # grab the current axis for left side
ax.set_xticks([2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,9e-6,10e-6,11e-6,12e-6]) # choose which x locations to have ticks
ax.set_xticklabels([2,3,4,5,6,7,8,9,10,11,12]) 



plt.xlabel('Time(microsecond)')
plt.ylabel('Amplitude(V)')
plt.title('P-wave', x=0.94,y=0.99)
# plt.grid()
plt.show()