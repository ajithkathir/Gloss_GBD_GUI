#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 10:49:08 2019

@author: ajith
"""
import matplotlib
import astropy.io.fits as fit
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
import os
import tkinter.messagebox
from tkinter import filedialog


def UploadAction(event=None):
    filename = filedialog.askopenfilename()
    print('Selected:', filename)

gloss= tk.Tk()
gloss.geometry('500x300')
gloss.title('GLOSS plotting tool')
button = tk.Button(gloss, text='Browse', command=UploadAction)
button.pack()



def get_hours(time_str):
    h, m, s = time_str.split(':')
    return int(h) + int(m)/ 60.0 + int(s)/3600.0
#gloss_FITS_Filename="/home/gloss/GLOSS_fitsData/2018/1811/GLOSS_20181109T023020To20181109T114305.fits"
filename = filedialog.askopenfilename(title='Select the input fits file')
data =fit.open(filename)
Output_image_folder=os.getcwd()
hdr=data[0].header
time_freq=data[1].data


imdata=data[0].data
z=[]
t=[]

 #time axis data 
obser_start=hdr[11]
obser_stop=hdr[12]

 #convert hh:mm:ss to hours
starthr=get_hours(obser_start)
stophr=get_hours(obser_stop)
tt=np.array(time_freq.TIME)

#tt=np.linspace(starthr,stophr,hdr[3])


ff=np.array(time_freq.FREQ)
ff=ff[0:401]



for i in range(35,435):
   z.append(i)   
#fig=plt.figure(facecolor='black')  
fig=plt.figure(figsize=( 10.45,5.28,), dpi=100,facecolor='black')
fig.subplots_adjust(0.07, 0.1, 1.09, 0.90, 0.2,0.2)

ax = fig.add_subplot(111)
#ax.axes([0.08, 0.08, 0.94-0.08, 0.94-0.08])
x = plt.pcolormesh(tt, ff,imdata, cmap='jet', vmin=-50.1, vmax=-8)
#x=plt.imshow(imdata,aspect='auto',cmap='jet',origin='lower',extent=[2.5,12.5,35,435])
#=plt.imshow(imdata,aspect='auto',cmap='jet',origin='lower')
plt.axis([2.5, 11.5, ff[0], ff[400]])
plt.xticks(np.arange(round(starthr,1),round(stophr,1), 0.5))
#x=plt.imshow(imdata,aspect='auto',cmap='jet',origin='lower',extent=[2.5,11.5,35,435])
#plt.xticks(np.arange(0.0, 13.0, 1.0))

y=plt.colorbar(x,pad=0.05,aspect=50)
y.ax.tick_params(colors='white',labelsize='6',width=0.1)
y.set_label('Amplitude (dBm)',color='white')

#ax.subplot_adjust()
ax.spines['bottom'].set_color('white')
ax.spines['left'].set_color('white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
plt.xticks(fontsize='6')
plt.yticks(fontsize='6')
plt.yticks(np.arange(round(ff[0],1),round(ff[400]+0.1,1), 25))

plt.title('Gauribidanur LOw frequency Solar Spectrum'+' -- '+hdr[8],fontsize="15",color='white')

#plt.title('Gauribidanur Low frequency Solar Spectrum'+' -- '+hdr[8],fontsize="15",color='white')
r=plt.xlabel('Time(UT)',fontsize="10",color='white')

plt.ylabel('Frequency(MHz)',fontsize="10",color='white')
fig.savefig(Output_image_folder+'GBD_DSPEC_'+str(hdr[8]).replace('/',"") +'.jpeg', dpi=100, facecolor='black', edgecolor='w',orientation='landscape', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)  


#plt.show()



#if fnmatch.fnmatch(filename, '*.fits'):    
#    messagebox.showerror('Error','File not supported')


gloss.mainloop()





