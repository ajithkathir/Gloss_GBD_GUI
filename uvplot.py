#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 16:46:08 2019

@author: ajith
"""
import imageio
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
from itertools import permutations
from scipy import signal
import csv
import sys

zero=np.zeros(32)

gauss=signal.gaussian(641,std=128)

omega=0.000073   #frequncy of earth rotation

f=80*(10**6)
lat=np.radians(14)


vis_cos_2d=np.zeros([640,32])
vis_sin_2d=np.zeros([640,32])



array=np.loadtxt(open("arraycoordGRAPH.csv", "rb"), delimiter=",", skiprows=1)
#print(np.shape(array))


#a1=(int(input('Enter ant1: '))-1)
#a2=(int(input('Enter ant2: '))-1)

#fringecos=np.zeros([640,32])

#fringesin=np.zeros([640,32])
u=np.zeros([64,64])
v=np.zeros([64,64])

E=np.linspace(1,32,32)
W=np.linspace(33,64,32)
for i in range(1,64,1):
    a1=i
    for j in range(1,64,1):
        a2=j
        if a1 in E and a2 in E:
            continue
        if a1 in W and a2 in W:
            continue
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
        d=np.radians(40)
        #hd=np.linspace(-7,7,641)
        h=np.radians(0)
            
                
        x1=(-1)*(N1)*np.sin(lat)+(U1)*np.sin(lat)
        x2=(-1)*(N2)*np.sin(lat)+(U2)*np.sin(lat)
        
        y1=E1
        y2=E2
        
        z1=(np.array(N1))*(np.cos(lat))
        z2=(np.array(N2))*(np.cos(lat))
        
        Bx=x2-x1
        By=y2-y1
        Bz=z2-z1
        
        
        u[i,j]=(1/lamda)*(((Bx)*np.sin((h)))+((By)*np.cos(h)))
        
        v[i,j]=(1/lamda)*(((-1)*(Bx)*np.sin(d)*np.cos(h))+(By*np.sin(d)*np.sin(h))+(Bz*np.cos(d)))
        w=(1/lamda)*((Bx*np.cos(d)*np.cos(h))+((-1)*(By)*np.cos(d)*np.sin(h))+(Bz*np.sin(d)))
        
        
plt.subplot(121)        
plt.scatter(u,v,s=5)

#print(u)
#print(v)

cur_axes = plt.gca()
cur_axes.axes.get_xaxis().set_visible(False)
cur_axes.axes.get_yaxis().set_visible(False)
plt.savefig("my_image.png")




uv1 =imageio.imread('my_image.png')
uv2=np.array(uv1)
uv3=uv2[45:245,65:375,2]
uv4=uv3-255
print(uv4.shape)
np.set_printoptions(threshold=sys.maxsize)
print(uv4)
#plt.figure(2)
#plt.imshow(uv3)

f=np.fft.fft2(uv4)
fshift=np.fft.fftshift(f)
img=(np.abs(fshift))

plt.subplot(122)
plt.imshow(img)



#print(((-1)*(By)*np.cos(d)*np.sin(h)))
#taug=w/c


#fringecos=(np.cos((2*(np.pi)*(f))*taug))*gauss
#fringesin=(np.sin((2*(np.pi)*(f))*taug))*gauss

#vis_cos_2d[:,i-32]=fringecos  
#vis_sin_2d[:,i-32]=fringesin

#plt.subplot(121)
#plt.plot(fringecos,'b')
#plt.plot(fringesin,'r')



#phase=(np.arctan2(fringesin,fringecos))*57.3



#plt.subplot(122)
#plt.plot(h,phase)
#plt.plot(h,taugE,'r')
#plt.plot(h,tau,'b')
#print(np.shape(tau_E))
#print(np.shape(tau_N))
#plt.show()

#cen_indx=320
#print("cosine values around transit..")
#print(fringecos[cen_indx-2],fringecos[cen_indx-1],fringecos[cen_indx],fringecos[cen_indx+1],fringecos[cen_indx+2])

#print("sine values around transit..")
#print(fringesin[cen_indx-2],fringesin[cen_indx-1],fringesin[cen_indx],fringesin[cen_indx+1],fringesin[cen_indx+2])

#print("Phase values around transit..")
#print(phase[cen_indx-2],phase[cen_indx-1],phase[cen_indx],phase[cen_indx+1],phase[cen_indx+2])
