#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 10:49:08 2019

@author: Ajith

@contributors: Indrajit, Kathiravan, Anshu

"""
from tkinter.ttk import *
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')

from astropy.io import fits
import numpy as np
import tkinter as tk
from tkinter import Tk,Frame,messagebox,Menu,Text,filedialog
from tkinter import *
import os
import datetime
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from sympy import *
import xlwt as xl


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def time():
    global now,nowdt
    no=str(datetime.datetime.now()) #.strftime('%H:%M:%S,%f')[:-7])
    now=no[0:19]
    nowdt=no[0:10]
    return now,nowdt



def OpenFile():
    global file

    file=filedialog.askopenfile(initialdir='pwd',title='select fits file',filetypes = (("fits files","*.fits"),("all files","*.*")))
    if '.fits' in file.name:
        text.config(state=NORMAL)
        time()
        text.insert(tk.END,now + '    You selected a fits file correctly'+'\n')
        text.see(tk.END)
        text.config(state=DISABLED)
        data =fits.open(file.name)
        hdr=data[0].header
        
        #Output_image_folder=os.getcwd()
        if 'GLOSS' in hdr['OBJECT']:
            enablingplot()
            text.config(state=NORMAL)
            time()
            text.insert(tk.END, now +'    You selected the correct GLOSS fits file'+'\n')
            text.see(tk.END)
            text.config(state=DISABLED)
            return file.name
        else:
            text.config(state=NORMAL)
            time()
            text.insert(tk.END,now +'    Error ######## Only GLOSS fits supported'+'\n','warning')
            text.see(tk.END)
            text.config(state=DISABLED)
            messagebox.showerror('Error','Only GLOSS fits supported')
            
              
    else:
        text.config(state=NORMAL)
        time()
        text.insert(tk.END,now +'    Extension Error ######## You chose a wrong file!! Please choose a fitsfile'+'\n','warning')
        text.see(tk.END)
        text.config(state=DISABLED)
        messagebox.showerror('Extension Error','You chose a wrong file!! Please choose a fitsfile')
        
                 
    return None
    

def enablingplot():
    menu.entryconfigure("Plot", state=tk.NORMAL)
def disablingplot():
    menu.entryconfigure("Plot", state=tk.DISABLED)

def enablinganalysis():
    menu.entryconfigure("Analysis", state=tk.NORMAL)
def disablinganalysis():
    menu.entryconfigure("Analysis", state=tk.DISABLED)
    
def enablingclear():
    menu.entryconfigure("Clear", state=tk.NORMAL)
def disablingclear():
    menu.entryconfigure("Clear", state=tk.DISABLED)    

def enablingsave():
    mfile.entryconfigure("Save",state=tk.NORMAL)
def disablingsave():
    mfile.entryconfigure("Save", state=tk.DISABLED)    


def enablingstart():
    manalys.entryconfigure("Start", state=tk.NORMAL)    
def disablingstart():
    manalys.entryconfigure("Start", state=tk.DISABLED)   

def enablinganalysn():
    manalys.entryconfigure("DiffFactor", state=tk.NORMAL)    
def disablinganalysn():
    manalys.entryconfigure("DiffFactor", state=tk.DISABLED)   
   
    
def enablingselect():
    manalys.entryconfigure("Select points", state=tk.NORMAL)    
def disablingselect():
    manalys.entryconfigure("Select points", state=tk.DISABLED)   
    

def enablingdelete():
    manalys.entryconfigure("Reset and Reselect", state=tk.NORMAL)    
def disablingdelete():
    manalys.entryconfigure("Reset and Reselect", state=tk.DISABLED)    



    
def readnplot():
    
    text.config(state=NORMAL)
    time()
    text.insert(tk.END,now +'    Plotting.........'+'\n')
    text.see(tk.END)
    text.config(state=DISABLED)
    
    global fp,t
    fp=[]
    t=[]

    def get_hours(time_str):
        h, m, s = time_str.split(':')
        return int(h) + int(m)/ 60.0 + int(s)/3600.0
   
        
    data =fits.open(file.name)
    hdr=data[0].header
    
    #Output_image_folder=os.getcwd()
    
    time_freq=data[1].data
    
    global imdata
    imdata=data[0].data
    z=[]
    t=[]
    
     #time axis data 
    obser_start=hdr[11]
    obser_stop=hdr[12]
    
    global tt,tt1,ff
    #convert hh:mm:ss to hours
    starthr=get_hours(obser_start)
    stophr=get_hours(obser_stop)
    tt=np.array(time_freq.TIME)
    tt=tt[0:imdata.shape[1]]
    
    
    ff=np.array(time_freq.FREQ)
    ff=ff[0:imdata.shape[0]]
    
    
    
    
    for i in range(35,435):
       z.append(i)   

    global fig
    fig=plt.figure(figsize=(15,6), dpi=100,facecolor='black')

    global ax
    ax = plt.subplot2grid((4, 4), (0,0), rowspan=3,colspan=3)
    

    x = plt.pcolormesh(tt, ff,imdata, cmap='jet', vmin=-50.1, vmax=-8)

    plt.axis([tt[0], tt[imdata.shape[1]-1], ff[0], ff[imdata.shape[0]-1]])
    
    def on_xlims_change(axes):
        global centretime
        centretime=(ax.get_xlim()[0]+ax.get_xlim()[1])/2
        global ax1, ax2 #,ax3
        ax1.cla()

        ax1.set_title("Freq (in MHz) profile "+str(round(centretime,4))+"UT",fontsize="12",color="white")
        ax1.spines['bottom'].set_color('white')
        ax1.spines['right'].set_color('white')
        ax1.tick_params(axis='x', colors='white')
        ax1.tick_params(axis='y', colors='white')
        ax1.plot(imdata[:,int(imdata.shape[1]/2)],ff)   
        ax1.invert_xaxis()
        ax1.set_xlabel('Amplitude (dBm)',color='white')

    
    
    def on_ylims_change(axes):
        global centrefreq
        centrefreq=(ax.get_ylim()[0]+ax.get_ylim()[1])/2
        ax2.cla()

        ax2.set_xlabel("Time (in UT) profile at "+str(round(centrefreq,4))+"MHz",fontsize="12",color="white")
        ax2.spines['bottom'].set_color('white')
        ax2.spines['left'].set_color('white')
        ax2.tick_params(axis='x', colors='white')
        ax2.tick_params(axis='y', colors='white')
        ax2.set_ylabel('Amplitude (dBm)',color='white')

        ax2.plot(tt,imdata[int(imdata.shape[0]/2),:])
        
        
    
    ax.callbacks.connect('xlim_changed', on_xlims_change)
    ax.callbacks.connect('ylim_changed', on_ylims_change)
    

    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')

    plt.xticks(fontsize='10')
    plt.yticks(fontsize='10')

    
    plt.title('Gauribidanur LOw frequency Solar Spectrum'+' -- '+hdr[8],fontsize="14",color='white')
   
    plt.ylabel('Frequency(MHz)',fontsize="12",color='white')
    plt.grid(True)

    global ax1, ax2 #,ax3
    ax1 = plt.subplot2grid((4, 4), (0, 3), rowspan=3,sharey=ax)
    ax1.set_title("Freq (in MHz) profile "+str(round(tt[int(tt.shape[0]/2)],5))+"UT",fontsize="12",color="white")
    ax1.spines['bottom'].set_color('white')
    ax1.spines['right'].set_color('white')
    ax1.tick_params(axis='x', colors='white')
    ax1.tick_params(axis='y', colors='white')
    ax1.plot(imdata[:,int(imdata.shape[1]/2)],ff)  
    ax1.invert_xaxis()
    ax1.set_xlabel('Amplitude (dBm)',color='white')
    ax1.yaxis.set_ticks_position('right')

    ax2 = plt.subplot2grid((4, 4), (3, 0), colspan=3,sharex=ax)
    ax2.spines['bottom'].set_color('white')
    ax2.spines['left'].set_color('white')
    ax2.tick_params(axis='x', colors='white')
    ax2.tick_params(axis='y', colors='white')
    ax2.set_ylabel('Amplitude (dBm)',color='white')
    

    ax2.plot(tt,imdata[int(imdata.shape[0]/2),:])
    ax2.set_xlabel("Time (in UT) profile at "+str(round(ff[int(ff.shape[0]/2)],5))+"MHz",fontsize="12",color="white")
    
    
    
    plt.gcf().canvas.draw()
    global canvas
    canvas = FigureCanvasTkAgg(fig, master=gloss)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    global toolbar
    toolbar = NavigationToolbar2Tk(canvas,gloss)
    toolbar.update()
    canvas.get_tk_widget().pack(anchor=tk.S,fill=tk.BOTH, expand=True)

    
    enablinganalysis()
    enablingclear()
    enablingsave()
    text.config(state=NORMAL)
    time()
    text.insert(tk.END,now +'    GLOSS spectrum plotted'+'\n')
    text.see(tk.END)
    text.config(state=DISABLED)

    
def clear():
    if messagebox.askyesno("Clear plot","Do you really want to clear existing plot?"):
        
        canvas.get_tk_widget().pack_forget()
        canvas._tkcanvas.pack_forget()       

        toolbar.destroy()
        text.config(state=NORMAL)
        time()
        text.insert(tk.END,now +'    cleared existing plot'+'\n') 
        text.see(tk.END)
        text.config(state=DISABLED)
        disablingplot()
        disablingclear()
        disablinganalysis()
        disablingsave()
    else:
        return None
    

    
    
def save():
    f=filedialog.asksaveasfilename(initialdir = 'pwd',title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    text2save = str(text.get(1.0, END)) # starts from `1.0`, not `0.0`
    f.write(text2save)
    text.config(state=NORMAL)
    time()
    text.insert(tk.END,now +'    saving the file in working directory'+'\n')
    text.see(tk.END)
    text.config(state=DISABLED)
    f.close() 
    
        
def close():
    if messagebox.askyesno("GLOSS plotting tool","Do you really want to exit?"):
        text.config(state=NORMAL)
        time()
        text.insert(tk.END,now +'    killing the window'+'\n')
        text.see(tk.END)
        text.config(state=DISABLED)
        txt=open("GLOSS--"+now+".log","w")
        txt.write(text.get(1.0,tk.END))
        txt.close()
        gloss.destroy()
    else:
        return None








global stop,evnts
stop=False
evnts=[]



def analysisclick():
     
    if stop==False:
        def onclick(event):
            

            enablingstart()
            xdata=event.xdata
            ydata=event.ydata
            
            if xdata==None:
                time()
                text.config(state=NORMAL)
                text.insert(tk.END,now+"    Warning ######### Select points only from the figure"+"\n","warning")
                text.see(tk.END)
                text.config(state=DISABLED)
            
            elif xdata>1 and ydata>1:
                time()
                text.config(state=NORMAL)
                text.insert(tk.END,now+'    Time selected : '+str(round(xdata,5))+'\n')
                text.insert(tk.END,now+'    Freq selected : '+str(round(ydata,5))+'\n')
                text.see(tk.END)
                text.config(state=DISABLED)
                global fp1,t1,t2
                fp1=(round(ydata,5))
                t1=(round(xdata,5))
                t2=None

                global x1,y1
                x1=event.x

                y1=event.y

                ax1.cla()
                ax2.cla()
                
                time_indx=np.where(tt==find_nearest(tt,xdata))

                freq_indx=np.where(ff==find_nearest(ff,ydata))
                
                cut1=np.linspace(0,(imdata.shape[0]-1),imdata.shape[0])
                cuttime=imdata[:,time_indx[0][0]]
                
                ax1.plot(cuttime,ff)
                ax1.set_title("Freq profile at "+str(round(xdata,5))+" hrs UTC",fontsize="12",color="white")
                ax1.set_xlabel('Amplitude (dBm)',color='white')
                ax1.spines['bottom'].set_color('white')
                ax1.spines['right'].set_color('white')
                ax1.tick_params(axis='x', colors='white')
                ax1.tick_params(axis='y', colors='white')
                ax1.yaxis.set_ticks_position('right')
                
                cut2=np.linspace(0,(imdata.shape[1]-1),imdata.shape[1])
                cutfreq=imdata[freq_indx[0][0],:]
                
                ax2.plot(tt,cutfreq)
                ax2.set_xlabel("Time (UT) profile at "+str(round(ydata,5))+" MHz",fontsize="12",color="white")
                ax2.set_ylabel('Amplitude (dBm)',color='white')        
                ax2.spines['bottom'].set_color('white')
                ax2.spines['left'].set_color('white')
                ax2.tick_params(axis='x', colors='white')
                ax2.tick_params(axis='y', colors='white')
                
                
              
                enablingdelete()
              
            elif xdata>1 and ydata<1:
                #print(reddot)
                global x2,y2
                x2=event.x
                y2=event.y
                global xdata2
                xdata2=event.xdata
                t2=round(xdata2,5)
              
                t.append(t2)
                fp.append(fp1)
                time()
                text.config(state=NORMAL)
                text.insert(tk.END,now+'    Time updated : '+str(round(xdata2,5))+'\n')
                text.insert(tk.END,now+'  New updated value will be taken for calculation'+'\n')
                #text.insert(tk.END,now+'    Freq selected : '+str(round(ydata,5))+'\n')
                text.see(tk.END)
                text.config(state=DISABLED)
               
                fig.canvas.get_tk_widget().create_rectangle(x2-3,canvas.get_tk_widget().winfo_height() - (y1-3),x2+3,canvas.get_tk_widget().winfo_height() - (y1+3),fill="black",tags=["dot","dot3"])

                
                evnts.append(x2)
                evnts.append(canvas.get_tk_widget().winfo_height() - (y1))
                #else:
                    
            
                
            else:
                #print(reddot)
                return None
                
                
                
                #return None
        def update():
            if t2==None:
                t.append(t1)
                fp.append(fp1)
                time()
                text.config(state=NORMAL)
                text.insert(tk.END,now+'    Time not updated'+'\n')
                text.insert(tk.END,now+'  Old value will be taken for calculation'+'\n')
                #text.insert(tk.END,now+'    Freq selected : '+str(round(ydata,5))+'\n')
                text.see(tk.END)
                text.config(state=DISABLED)
                
            else:
                t.append(t2)
                fp.append(fp1)
                time()
                text.config(state=NORMAL)
                text.insert(tk.END,now+'    Time updated : '+str(round(xdata2,5))+'\n')
                text.insert(tk.END,now+'  New updated value will be taken for calculation'+'\n')
                #text.insert(tk.END,now+'    Freq selected : '+str(round(ydata,5))+'\n')
                text.see(tk.END)
                text.config(state=DISABLED)
                
                fig.canvas.get_tk_widget().create_rectangle(x2-3,canvas.get_tk_widget().winfo_height() - (y1-3),x2+3,canvas.get_tk_widget().winfo_height() - (y1+3),fill="black",tags=["dot","dot3"])
                global reddot
                reddot=0
                
            
        fig.canvas.mpl_connect('button_press_event', onclick)
        text.config(state=NORMAL)
        time()
        text.insert(tk.END,now+"    Please select points from figure"+"\n")
        text.see(tk.END)
        text.config(state=DISABLED)
    else:
        return None
    
    

    
def deletealldots():
        global stop
        stop=True
               
        fig.canvas.get_tk_widget().delete("dot")
        fig.canvas.get_tk_widget().delete("line")
        
        global fp,t,evnts
        fp=[]
        t=[] 
        evnts=[]
        #print(fp)
        time()
        text.config(state=NORMAL)
        text.insert(tk.END,now+"    All selected points deleted"+"\n")
        text.insert(tk.END,now+"    Frequency and Time values are empty"+"\n")
        text.insert(tk.END,now+"    Select points from figure"+"\n")
        text.see(tk.END)
        text.config(state=DISABLED)
        
        messagebox.showinfo('Selected points deleted', "Frequency and Time values are empty. Continue to Select points from figure")

        disablingstart() 
        disablingdelete()
 

def mediansub():
    global imdata
    med=np.median(imdata)
    imdata=imdata-med
    canvas.get_tk_widget().pack_forget()
    canvas._tkcanvas.pack_forget()       
    #mainFrame.pack_forget()       
    toolbar.destroy()
    text.config(state=NORMAL)
    time()
    text.insert(tk.END,now +'    cleared existing plot'+'\n') 
    text.see(tk.END)
    text.config(state=DISABLED)
    disablingplot()
    disablingclear()
    disablinganalysis()
    disablingsave()
    readnplot() 
    
    
    
def analysis(Sc):
    #print(evnts)
    fig.canvas.get_tk_widget().create_line(evnts,width=2,tags=["line"])
    
    global wb,sheet1
    wb = xl.Workbook()
                        
    sheet1 = wb.add_sheet("Sheet1", cell_overwrite_ok=True)
            
    sheet1.write(0,1,"Fp Vs Time")
    sheet1.write(0,2,"Drift rate Vs Time")
    sheet1.write(0,3,"Height Vs Time")
    sheet1.write(0,4,"Velocity Vs Time")
    sheet1.write(1,0,"Newkirk Model")                
    sheet1.write(2,0,"BaumbachAllen_Model1")                
    sheet1.write(3,0,"BaumbachAllen_Model2")                
    sheet1.write(4,0," Saito Model")                
    sheet1.write(5,0,"Sittler & Guhathakurta Current sheet model")                
    sheet1.write(6,0,"Sittler & Guhathakurta Polar coronal hole model") 
    sheet1.write(7,0,"Hybrid model") 
    sheet1.write(8,0,"LeBlanc model")            
    
    
    
    if len(fp)>1:    
        def newton_raphson(f=None):
            
            e = 1e-5# setting the tolerance value
            dr = e + 1
        
            r0 = 1# initially assumed value of x
            count = 0# setting counter to know the no of itrations taken
            p = zeros(1,300)
            while (np.abs(dr) > e):    # initialising the itration and continue until the error is less than tolerance
                dr1=eval('f')
                drnum=dr1.subs(r,r0)            # calculating dx, diff is used for finding the differentiation of the fuction
                dr2=eval('diff(f)')
                drden=dr2.subs(r,r0)
                dr=drnum/drden
                #print(dr)
                r0 = r0 - dr    # updating the value of x
                count = count + 1    # incrimenting the counter
                p[count]=r 
                if (count > 300):
                    print('Error...! Solution not converging !!! \\n')       # printing the error message
                    break
                
            
                    
            return r0
    
    
    
        
        
            
        r=Symbol('r')
        
        Ne = (np.divide(fp,(9 * 10 ** (-3))))** 2# electron density from plasma frequency
        
        #fp = sqrt(Ne(i))*(9*10^(-3)) 
        
        
        #Sc = 1 #tk.simpledialog.askfloat('Enter Enhancement factor','Enter Enhancement Factor')
    
        if Sc != None:
                   
            Hsa=np.zeros(len(fp))
            Hnk=np.zeros(len(fp))
            Hba1=np.zeros(len(fp))
            Hlb=np.zeros(len(fp))
            Hhy=np.zeros(len(fp))
            Hsgc=np.zeros(len(fp))
            Hba2=np.zeros(len(fp))
            Hsgp=np.zeros(len(fp))
            
            
            
            
            for i in range(len(fp)):
                #Newkirk model
            
            
                Hnk[i]=4.321/np.log10(Ne[i]/(Sc*4.2*10**4))
               
               
                # Baumbach-Allen model 1
            
                BA1 = (1.55*(r** (-6)) + 2.99 * (r ** (-16)) + 0.036 * (r ** (-1.5))) - (Ne[i] / Sc * (10 ** (-8)))

            
                Hba1[i]=newton_raphson(BA1) #,itr[2]
            
                # Baumbach-Allen model 2
            
                Hba2[i] = (Sc * 4 * 10 ** 8 / Ne[i])** (1 / 9)
            
                # Saito model
            
                C1 = 1.36 * 10 ** 6
                d1 = 2.14; #print (d1)
            
                C2 = 1.68 * 10 ** 8
                d2 = 6.13; #print (d2)
            
            
                SA = C1 * r ** (-d1) + C2 * r ** (-d2) - Ne[i] / Sc
               
                
                Hsa[i]=newton_raphson(SA)

            
                # Hybrid model
            
                HYB = 15.45 * r ** (-16) + 3.16 * r ** (-6) + r ** (-4) + 0.0033 * r ** (-2) - Ne[i] / (Sc * (10 ** 8))
                Hhy[i]=newton_raphson(HYB)

            
                # Sittler & Guhathakurta model 
            
                zr = 1 / (1 + r)
            
                ap1 = 0.001272
                ap2 = 4.8039; #print (ap2)
                ap3 = 0.29696
                ap4 = -7.1743
                ap5 = 12.321
                #% polar coronal hole
            
                ac1 = 0.0032565
                ac2 = 3.6728; #print (ac2)
                ac3 = 4.8947
                ac4 = 7.6123
                ac5 = 5.9868
                #% current sheet model
            
                rho0 = 1.673 * 10 ** -14 / (2 * 10 ** -24)#% gm/cm^3 is converted to e/cm^3
            
            
                SGp = rho0 * (ap1 * (zr ** 2) * exp(ap2 * zr) * (1 + ap3 * zr + ap4 * zr ** 2 + ap5 * zr ** 3)) - Ne[i] / Sc#% polar coronal hole
            
                SGc = rho0 * (ac1 * (zr ** 2) * exp(ac2 * zr) * (1 + ac3 * zr + ac4 * zr ** 2 + ac5 * zr ** 3)) - Ne[i] / Sc#% current sheet model
            
                Hsgp[i]=newton_raphson(SGp)
              
                Hsgc[i]=newton_raphson(SGc)
              
                #Leblanc
            
                LB = 3.3e5 * r ** -2 + 4.1e6 * r ** -4 + 8e7 * r ** -6 - Ne[i] / Sc
            
                Hlb[i]=newton_raphson(LB)
                      
                
            drift=np.zeros(len(fp)-1)
            for i in range(len(fp)-1):
                drift[i]=((fp[i+1]-fp[i])/(t[i+1]-t[i]))
            heights=(Hnk, Hba1,Hba2, Hsa,Hsgc,Hsgp,Hhy, Hlb)
            #print(drift)
            
            t_av=np.zeros(len(fp)-1)
            for i in range(len(fp)-1):
                t_av[i]=((t[i+1]+t[i])/2)
            
            Vsa=np.zeros(len(fp)-1)
            Vnk=np.zeros(len(fp)-1)
            Vba1=np.zeros(len(fp)-1)
            Vlb=np.zeros(len(fp)-1)
            Vhy=np.zeros(len(fp)-1)
            Vsgc=np.zeros(len(fp)-1)
            Vba2=np.zeros(len(fp)-1)
            Vsgp=np.zeros(len(fp)-1)
            for i in range(len(fp)-1):
                Vnk[i]=(Hnk[i+1]-Hnk[i])/(t[i+1]-t[i])
                Vba1[i]=(Hba1[i+1]-Hba1[i])/(t[i+1]-t[i])
                Vba2[i]=(Hba2[i+1]-Hba2[i])/(t[i+1]-t[i])
                Vsa[i]=(Hsa[i+1]-Hsa[i])/(t[i+1]-t[i])
                Vsgc[i]=(Hsgc[i+1]-Hsgc[i])/(t[i+1]-t[i])
                Vsgp[i]=(Hsgp[i+1]-Hsgp[i])/(t[i+1]-t[i])
                Vhy[i]=(Hhy[i+1]-Hhy[i])/(t[i+1]-t[i])
                Vlb[i]=(Hlb[i+1]-Hlb[i])/(t[i+1]-t[i])
                
                
            velocities=(Vnk,Vba1,Vba2,Vsa,Vsgc,Vsgp,Vhy,Vlb)
            
            las=filedialog.askopenfile(parent=mainwin,initialdir='pwd',title='select lasco file for overplotting',filetypes = (("yht files","*.yht"),("all files","*.*")))            
            #print(las.name)
            lasco=np.loadtxt((las.name),comments="#",delimiter="\t",dtype=str)
            
            ht=np.zeros(lasco.shape[0])
            
            for i in range(lasco.shape[0]):                 
                test=(lasco[i].rsplit(" "))
                if test[1]=="":
                    ht[i]=float(test[2])
                else:
                    ht[i]=float(test[1])
                
            tim=np.zeros(lasco.shape[0])
            for i in range(lasco.shape[0]):
                test=(lasco[i].rsplit(" "))        
                if test[1]=="":
                    tim[i]=hms_to_hr(test[4])
                else:
                    tim[i]=hms_to_hr(test[3])
            
            def NewkirkModel():
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    NewkirkModel height values="+str(heights[0])+"\n")
                text.insert(tk.END,now+"    NewkirkModel velocity values="+str(velocities[0])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f=plt.figure(figsize=(10,10),num="Newkirk Model")
                
                plt.subplot(321)
                plt.scatter(t,fp)
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(1,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(1,2,str(drpoly))
                plt.plot(t_av,drplt)
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[0])
                htfit=np.polyfit(t, heights[0], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(1,3,str(htpoly))
                plt.plot(t,htplt)
                
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[0])
                         
                vfit=np.polyfit(t_av, velocities[0], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(1,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
            
                plt.subplot(313)
                plt.plot(t,heights[0],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global newfit
                
                newfit=[fppoly,drpoly,htpoly,vpoly]
                
                return newfit
                
            def BaumbachAllen_Model1():
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Baumbach-Allen Model1 height values="+str(heights[1])+"\n")
                text.insert(tk.END,now+"    Baumbach-Allen Model1 velocity values="+str(velocities[1])+"\n")
                
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="Baumbach-Allen Model - 1")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(2,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(2,2,str(drpoly))
                plt.plot(t_av,drplt)                
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[1])
                htfit=np.polyfit(t, heights[1], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(2,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[1])
                vfit=np.polyfit(t_av, velocities[1], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(2,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[1],marker='o', linestyle='--')
                text.config(state=NORMAL)
                
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global ba1fit
                
                ba1fit=[fppoly,drpoly,htpoly,vpoly]
                
                return ba1fit
                
            def BaumbachAllen_Model2():
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Baumbach-Allen Model2 height values="+str(heights[2])+"\n")
                text.insert(tk.END,now+"    Baumbach-Allen Model2 velocity values="+str(velocities[2])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="Baumbach-Allen Model - 2")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(3,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(3,2,str(drpoly))
                plt.plot(t_av,drplt)
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[2])
                htfit=np.polyfit(t, heights[2], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(3,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[2])
                vfit=np.polyfit(t_av, velocities[2], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(3,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[2],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                
                global ba2fit
                ba2fit=[fppoly,drpoly,htpoly,vpoly]
                
                return ba2fit
                
                
                
            def SaitoModel():    
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Saito Model height values="+str(heights[3])+"\n")
                text.insert(tk.END,now+"    Saito Model velocity values="+str(velocities[3])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="Saito Model")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(4,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(4,2,str(drpoly))
                plt.plot(t_av,drplt)                
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[3])
                htfit=np.polyfit(t, heights[3], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(4,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[3])
                vfit=np.polyfit(t_av, velocities[3], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(4,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[3],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global saitofit
                saitofit=[fppoly,drpoly,htpoly,vpoly]
                
                return saitofit
                
            def SittleGuhathakurtaCurrent_model():
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Sittler & Guhathakurta Current Sheet model height values="+str(heights[3])+"\n")
                text.insert(tk.END,now+"    Sittler & Guhathakurta Current Sheet model velocity values="+str(velocities[3])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="Sittler & Guhathakurta Current Sheet Model")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(5,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(5,2,str(drpoly))
                plt.plot(t_av,drplt)
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[4])
                htfit=np.polyfit(t, heights[4], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(5,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[4])
                vfit=np.polyfit(t_av, velocities[4], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(5,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[4],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global sgcfit
                sgcfit=[fppoly,drpoly,htpoly,vpoly]
                
                return sgcfit
                
                
            def SittleGuhathakurtaPolar_model():
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Sittler & Guhathakurta Polar coronal hole model  height values="+str(heights[5])+"\n")
                text.insert(tk.END,now+"    Sittler & Guhathakurta Polar coronal hole model  velocity values="+str(velocities[5])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="Sittler & Guhathakurta Polar coronal hole Model")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(6,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(6,2,str(drpoly))
                plt.plot(t_av,drplt)
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[5])
                htfit=np.polyfit(t, heights[5], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(6,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[5])
                vfit=np.polyfit(t_av, velocities[5], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(6,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[5],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global sgpfit
                sgpfit=[fppoly,drpoly,htpoly,vpoly]
                
                return sgpfit
                
                
            def Hybridmodel():
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Hybrid model height values="+str(heights[6])+"\n")
                text.insert(tk.END,now+"    Hybrid model velocity values="+str(velocities[6])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="Hybrid Model")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(7,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(7,2,str(drpoly))
                plt.plot(t_av,drplt)
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[6])
                htfit=np.polyfit(t, heights[6], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(7,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[6])
                vfit=np.polyfit(t_av, velocities[6], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(7,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[6],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global hybridfit
                hybridfit=[fppoly,drpoly,htpoly,vpoly]
                
                return hybridfit
                
                
            def LeblancModel():
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Leblanc Model height values="+str(heights[7])+"\n")
                text.insert(tk.END,now+"    Leblanc Model velocity values="+str(velocities[7])+"\n")
                text.see(tk.END)
                text.config(state=DISABLED)
                f=plt.figure(figsize=(10,10),num="LeBlanc Model")
                
                plt.subplot(321)
                plt.scatter(t,fp)           #,marker='o',linestyle="--")
                fpfit=np.polyfit(t, fp, n)
                fppoly=np.poly1d(fpfit)
                ys=fppoly(t)
                sheet1.write(8,1,str(fppoly))
                plt.plot(t,ys)
                plt.title("Freq Vs Time")
                
                plt.subplot(322)
                plt.scatter(t_av,drift)           
                drfit=np.polyfit(t_av, drift, n-1)
                drpoly=np.poly1d(drfit)
                drplt=drpoly(t_av)
                sheet1.write(8,2,str(drpoly))
                plt.plot(t_av,drplt)
                plt.title("Drift rate")
                
                plt.subplot(323)
                plt.scatter(t,heights[7])
                htfit=np.polyfit(t, heights[7], n)
                htpoly=np.poly1d(htfit)
                htplt=htpoly(t)
                sheet1.write(8,3,str(htpoly))
                plt.plot(t,htplt)
                plt.title("Height VS Time plot")
                
                plt.subplot(324)
                plt.scatter(t_av,velocities[7])
                vfit=np.polyfit(t_av, velocities[7], n-1)
                vpoly=np.poly1d(vfit)
                vplt=vpoly(t_av)
                sheet1.write(8,4,str(vpoly))
                plt.plot(t_av,vplt)
                plt.title("Velocity plot")
                
                plt.subplot(313)
                plt.plot(t,heights[7],marker='o', linestyle='--')
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
                text.see(tk.END)
                text.config(state=DISABLED)                
                
                plt.plot(tim,ht,marker="o",linestyle="--")    
                plt.title("LASCO Height Vs Time plot overplotted with GLOSS")
                plt.xlabel("Time in hrs")
                plt.ylabel("Height")
                
                text.config(state=NORMAL)
                time()
                text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
                text.see(tk.END)
                text.config(state=DISABLED)
                
                f.show()
                global leblancfit
                leblancfit=[fppoly,drpoly,htpoly,vpoly]
                
                return leblancfit
                
            
            n=len(fp)-2
            
            NewkirkModel()
            BaumbachAllen_Model1()
            BaumbachAllen_Model2()
            SaitoModel()
            SittleGuhathakurtaCurrent_model()   
            SittleGuhathakurtaPolar_model()
            Hybridmodel()
            LeblancModel()
            
            time()
            wb.save("GLOSS--"+nowdt+"--fit_table.xlsx")
            
             

    else:    
        text.config(state=NORMAL)
        time()
        text.insert(tk.END,now+"    Inadequate points for calculation")
        text.insert(tk.END,now+"    Select atleast two points")
        text.see(tk.END)
        text.config(state=DISABLED)
        messagebox.showerror('Inadequate points','Select atleast two points')
        





def analysis1():
    global mainwin
    mainwin=gloss
    analysis(1)
    enablinganalysn()


def analysis_n():
        
    def clearval():
        e1.delete(0,'end')
    
    def calc():
        
        Sc1=e1.get()
        if Sc1.isnumeric() == True:
            Sc1=float(Sc1)
            global mainwin
            mainwin=master
            analysis(Sc1)
            
            tree.insert('','end',text="Newkirk Model",value=(str(newfit[0]),str(newfit[1]),str(newfit[2]),str(newfit[3])),tags='T')
            tree.insert('','end',text="Baumbach-Allen Model 1",value=(str(ba1fit[0]),str(ba1fit[1]),str(ba1fit[2]),str(ba1fit[3])),tags='T')
            tree.insert('','end',text="Baumbach-Allen Model 2",value=(str(ba2fit[0]),str(ba2fit[1]),str(ba2fit[2]),str(ba2fit[3])),tags='T')
            tree.insert('','end',text="Saito Model",value=(str(saitofit[0]),str(saitofit[1]),str(saitofit[2]),str(saitofit[3])),tags='T')
            tree.insert('','end',text="Sittler & Guhathakurta Current Sheet model",value=(str(sgcfit[0]),str(sgcfit[1]),str(sgcfit[2]),str(sgcfit[3])),tags='T')
            tree.insert('','end',text="Sittler & Guhathakurta Polar coronal hole",value=(str(sgpfit[0]),str(sgpfit[1]),str(sgpfit[2]),str(sgpfit[3])),tags='T')
            tree.insert('','end',text="Hybrid model",value=(str(hybridfit[0]),str(hybridfit[1]),str(hybridfit[2]),str(hybridfit[3])),tags='T')
            tree.insert('','end',text="LeBlanc model",value=(str(leblancfit[0]),str(leblancfit[1]),str(leblancfit[2]),str(leblancfit[3])),tags='T')
            
            tree.tag_configure('T', font='Times 15')
            
            
    
        else:
            text.config(state=NORMAL)
            time()
            text.insert(tk.END,now +'    Value Error ######## You should give a number!! Please enter a number'+'\n','warning')
            text.see(tk.END)
            text.config(state=DISABLED)
            messagebox.showerror('Value Error','You should give a number!! Please enter a number')
            
            
    master=tk.Tk()
    master.geometry("%dx%d+%d+%d" % (2000, 800, 100, 100))
    master.title("Fit equation table")
    #master.resizable(0,0)
    treeframe = tk.Frame(master)
    treeframe.pack(fill=BOTH,expand=True)
    
    entry=Label(treeframe,text="Enter Enhancement Factor",font=("Times",15))
    entry.pack(side=TOP,anchor=W)
    e1=Entry(treeframe,width=15,font=("Times",15))
    e1.pack(side=TOP,anchor=NW)
    
    
    b1=tk.Button(treeframe, text='Calculate fit',font=("Times",15),command=calc)
    b1.pack(side=TOP,anchor=W)
    b2=tk.Button(treeframe, text='Clear', font=("Times",15), command=clearval)
    b2.pack(side=TOP,anchor=W)
    
    style = Style(treeframe)
    style.configure("mystyle.Treeview.Heading", font=("Times", 15))
    style.configure("mystyle.Treeview", rowheight=60)    
    tree=Treeview(treeframe,style="mystyle.Treeview")
    
    
    
    tree["columns"]=("one","two","three","four")
    tree.column("#0", width=350, minwidth=200, stretch=tk.YES)
    tree.column("one", width=300, minwidth=200, stretch=tk.YES)
    tree.column("two", width=350, minwidth=200,stretch=tk.YES)
    tree.column("three", width=300, minwidth=200, stretch=tk.YES)
    tree.column("four", width=350, minwidth=200, stretch=tk.YES)
        
    tree.rowconfigure(1,weight=1)
    tree.heading("#0",text="Models",anchor=tk.W)
    tree.heading("one", text="Freq vs Time",anchor=tk.W)
    tree.heading("two", text="Drift rate vs Time",anchor=tk.W)
    tree.heading("three", text="Height vs Time",anchor=tk.W)
    tree.heading("four", text="Velocity vs Time",anchor=tk.W)
    
    
    tree.pack(side=tk.TOP,fill=tk.Y,expand=True)
    
    master.mainloop()
        
    
    
    
    
    
def hms_to_hr(s):
    t = 0
    for u in s.split(':'):
        t = 60 * t + int(u)
    return t/3600
    

def lasco():
    las=filedialog.askopenfile(initialdir='pwd',title='select lasco file',filetypes = (("yht files","*.yht"),("all files","*.*")))            
    #print(las.name)
    lasco=np.loadtxt((las.name),comments="#",delimiter="\t",dtype=str)
    
    ht=np.zeros(lasco.shape[0])
    
    for i in range(lasco.shape[0]):                 
        test=(lasco[i].rsplit(" "))
        if test[1]=="":
            ht[i]=float(test[2])
        else:
            ht[i]=float(test[1])
        
    tim=np.zeros(lasco.shape[0])
    for i in range(lasco.shape[0]):
        test=(lasco[i].rsplit(" "))        
        if test[1]=="":
            tim[i]=hms_to_hr(test[4])
        else:
            tim[i]=hms_to_hr(test[3])
    
    text.config(state=NORMAL)
    time()
    text.insert(tk.END,now+"    Plotting LASCO Height-Time ..... \n")
    text.see(tk.END)
    text.config(state=DISABLED)
    
    f=plt.figure()
    plt.plot(tim,ht,marker="o",linestyle="--")    
    plt.suptitle("LASCO Height Vs Time plot overplotted with GLOSS")
    plt.xlabel("Time in hrs")
    plt.ylabel("Height")
    text.config(state=NORMAL)
    time()
    text.insert(tk.END,now+"    LASCO Height-Time plotted \n")
    text.see(tk.END)
    text.config(state=DISABLED)
    f.show()
    
    
    
    
    
gloss= tk.Tk()
gloss.style = Style()
gloss.style.theme_use("clam")
gloss.minsize(1600,600)
gloss.title('GLOSS plotting tool')
scrollbar = Scrollbar(gloss)
scrollbar.pack(side=RIGHT, fill=Y)

mainFrame = tk.Frame(gloss)
mainFrame.pack()


text=Text(gloss,fg="black",width=65,height=100,yscrollcommand=scrollbar.set,font=('Times', 12))
text.pack(side=tk.RIGHT,anchor=tk.SE, expand=False)
text.tag_config('warning', foreground="red")
text.tag_config('title', foreground="black",font=('Times',20,"bold"))
time()
text.insert(tk.END,"                            LOGGER \n","title")
text.insert(tk.END,"Plotting tool opened at:"+now+"\n","title")
text.insert(tk.END,"--------------------------------------------------------- \n","title")
text.config(state=DISABLED)


menu = tk.Menu(gloss)

mfile = tk.Menu(menu)
mfile.add_command(label="Open", command=OpenFile)
mfile.add_command(label="Save", command=save,state=tk.DISABLED)
mfile.configure(font=("Times",18))

mplot = tk.Menu(menu)
mplot.add_command(label="Plot spectrum", command=readnplot)
mplot.configure(font=("Times",18))


manalys=tk.Menu(menu)
manalys.add_command(label="Median substraction",command=mediansub)

manalys.add_command(label="Select points",command=analysisclick)

manalys.add_command(label="Reset and Reselect",command=deletealldots,state=tk.DISABLED)

manalys.add_command(label="Start",command=analysis1,state=tk.DISABLED)
manalys.add_command(label="DiffFactor",command=analysis_n,state=tk.DISABLED)
manalys.configure(font=("Times",18))

mlasco=tk.Menu(menu)
mlasco.add_command(label="Plot LASCO Height Vs Time",command=lasco)
mlasco.configure(font=("Times",18))

menu.add_cascade(label='File', menu=mfile)
menu.add_cascade(label='Plot', menu=mplot, state=tk.DISABLED)
menu.add_cascade(label='Analysis',menu=manalys,state=tk.DISABLED)
menu.add_cascade(label='LASCO',menu=mlasco)
menu.add_command(label="Clear", command=clear,state=tk.DISABLED)
menu.add_command(label="Close", command=close)
menu.configure(background="black",foreground='white',font=("Times",18))    
gloss.config(menu=menu)
gloss.option_add('*Dialog.msg.font', 'Times 15')
messagebox.showinfo("Welcome to GLOSS plotting tool", "This tool is designed and developed for plotting GLOSS spectrograph data and analysis of typeII solar event")


gloss.mainloop()





