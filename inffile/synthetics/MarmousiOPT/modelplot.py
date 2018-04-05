#!/usr/bin
import numpy as np
import matplotlib.pyplot as plt
import pylab
import obspy
from obspy.core import read;
from matplotlib import pyplot;
from pylab import genfromtxt; 
from scipy import integrate;
from scipy.integrate import odeint;


# m01 = genfromtxt("Receiver01.txt");
# m01z = genfromtxt("Receiver01Z.txt");

# m02 = genfromtxt("Receiver02.txt");
# m02z = genfromtxt("Receiver02Z.txt");

# m03 = genfromtxt("Receiver03.txt");
# m03z = genfromtxt("Receiver03Z.txt");

# m39 = genfromtxt("Receiver39.txt");
# m39z = genfromtxt("Receiver39Z.txt");

mp = genfromtxt("Receiver01P.txt");
ms = genfromtxt("Receiver01S.txt");

# a1 = m01[:,1];
# t1 = m01[:,0];
# a1z = m01z[:,1];
# t1z = m01z[:,0]; 


# a2 = m02[:,1];
# t2 = m02[:,0];
# a2z = m02z[:,1];
# t2z = m02z[:,0];

# a3 = m03[:,1];
# t3 = m03[:,0];
# a3z = m03z[:,1];
# t3z = m03z[:,0];

# a39 = m39[:,1];
# t39 = m39[:,0];
# a39z = m39z[:,1];
# t39z = m39z[:,0];

mp = mp[:,1];
tp = mp[:,0];
ms = ms[:,1];
ts = ms[:,0]; 

fig = plt.figure(1)

# la1 = np.amax(a1)
# lla1 = np.amin(a1)
# la1z = np.amax(a1z)
# lla1z = np.amin(a1z)
# la39 = np.amax(a39)
# lla39 = np.amin(a39)
# la39z = np.amax(a39z)
# lla39z = np.amin(a39z)

la1 = np.amax(mp)
lla1 = np.amin(mp)
la1z = np.amax(ms)
lla1z = np.amin(ms)

lt1 = 3500*1*(10**(-8))

fig.suptitle('Receiver in X and Z component', fontsize=12)
# plt.subplot(221)
# pyplot.plot(t1,a1, label = "Receiver01_X");
# pyplot.legend();
# plt.axis([0, 5*10**(-5) , lla1, la1])
# plt.subplot(222)
# pyplot.plot(t1z,a1z, label = "Receiver01_Z");
# plt.axis([0, 5*10**(-5), lla1z, la1z])
# pyplot.legend();
# plt.subplot(223)
# pyplot.plot(t39,a39, label = "Receiver39_X");
# plt.axis([0, 5*10**(-5), la39, lla39])
# pyplot.legend();
# plt.subplot(224)
# pyplot.plot(t39z,a39z, label = "Receiver39_Z");
# plt.axis([0, 5*10**(-5), lla39z, la39z])
# pyplot.legend();

plt.subplot(221)
pyplot.plot(tp,mp, label = "UP");
pyplot.legend();
plt.axis([0, 5*10**(-5) , lla1, la1])
plt.subplot(222)
pyplot.plot(ts,ms, label = "US");
plt.axis([0, 5*10**(-5), lla1z, la1z])
pyplot.legend();

pyplot.show();