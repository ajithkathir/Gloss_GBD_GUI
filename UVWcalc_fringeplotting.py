#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 18:01:12 2019

@author: ajith
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from itertools import permutations
from scipy import signal
import csv

zero=np.zeros(32)

gauss=signal.gaussian(3600,std=600)

omega=0.000073   #frequncy of earth rotation

f=80*(10**6)
lat=np.radians(14)

array=np.loadtxt(open("arraycoordGRAPH.csv", "rb"), delimiter=",", skiprows=1)
#print(np.shape(array))


a1=(int(input('Enter ant1: '))-1)
a2=(int(input('Enter ant2: '))-1)


antcoord1=np.array(array[a1])
antcoord2=np.array(array[a2])


E1=antcoord1[0]
E2=antcoord2[0]
N1=antcoord1[1]
N2=antcoord2[1]
U1=antcoord1[2]
U2=antcoord2[2]


c=300000000    
lamda=c/f  #put lamda=1 for plotting in UVdist in units of metre
d=np.radians(0)
hd=np.linspace(-7.5,7.5,3600)
h=np.radians(hd)
    

#tau_E=[]
#tau_N=[]

#for i in range(len(perm_N)):

#basl=list(perm_N[i])

#x1=(np.negative(basl[0]))*(np.sin(lat))
#x2=(np.negative(basl[1]))*(np.sin(lat))
#y1=0
#y2=0

#z1=(np.array(basl[0]))*(np.cos(lat))
#z2=(np.array(basl[1]))*(np.cos(lat))

#Bx1=x2-x1
#By1=y2-y1
#Bz1=z2-z1


#u1=(1/lamda)*(((Bx1)*np.sin(h))+((By1)*np.cos(h)))
#v1=(1/lamda)*(((-1)*(Bx1)*np.sin(d)*np.cos(h))+(By1*np.sin(d)*np.sin(h))+(Bz1*np.cos(d)))
#w1=(1/lamda)*((Bx1*np.cos(d)*np.cos(h))+((-1)*(By1)*np.cos(d)*np.sin(h))+(Bz1*np.sin(d)))

#taug1=w1/c
#tau_E=taug1


#for i in range(len(perm_E)):
    
x1=(-1)*(N1)*np.sin(lat)+(U1)*np.sin(lat)
x2=(-1)*(N2)*np.sin(lat)+(U2)*np.sin(lat)

y1=E1
y2=E2

z1=(np.array(N1))*(np.cos(lat))
z2=(np.array(N2))*(np.cos(lat))

Bx=x2-x1
By=y2-y1
Bz=z2-z1


u=(1/lamda)*(((Bx)*np.sin((h)))+((By)*np.cos(h)))
v=(1/lamda)*(((-1)*(Bx)*np.sin(d)*np.cos(h))+(By*np.sin(d)*np.sin(h))+(Bz*np.cos(d)))
w=(1/lamda)*((Bx*np.cos(d)*np.cos(h))+((-1)*(By)*np.cos(d)*np.sin(h))+(Bz*np.sin(d)))
#print(((-1)*(By)*np.cos(d)*np.sin(h)))
taug=w/c
#print(w)

#tau=((By))*(1/c)*(np.sin(np.radians(h)))*(np.cos((d)))
#print(taugE)



fringecos=(np.cos((2*(np.pi)*(f))*taug))*gauss
fringesin=(np.sin((2*(np.pi)*(f))*taug))*gauss



plt.plot(h,fringecos,'b')
plt.plot(h,fringesin,'r')
#plt.plot(h,taugE,'r')
#plt.plot(h,tau,'b')

#print(np.shape(tau_E))
#print(np.shape(tau_N))
plt.show()








