#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 10:49:08 2019

@author: ajith
"""
from tkinter.ttk import *
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import tkinter as tk
from tkinter import * #Tk,Frame,messagebox,Menu,Text
import os
import datetime
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from sympy import *







def time():
    global now
    now=str(datetime.datetime.now())
    return now



def OpenFile():
    global file

    file=filedialog.askopenfile(initialdir='pwd',title='select fits file')
    if '.fits' in file.name:
        time()
        text.insert(tk.END,now + '...You selected a fits file correctly'+'\n')
        data =fits.open(file.name)
        hdr=data[0].header
        
        #Output_image_folder=os.getcwd()
        if 'GLOSS' in hdr['OBJECT']:
            enabling()
            time()
            text.insert(tk.END, now +'...You selected the correct GLOSS fits file'+'\n')
            return file.name
        else:
            messagebox.showerror('Error','Only GLOSS fits supported')
            time()
            text.insert(tk.END,now +'...Error ######## Only GLOSS fits supported'+'\n','warning')
              
    else:
        messagebox.showerror('Extension Error','You chose a wrong file!! Please choose a fitsfile')
        time()
        text.insert(tk.END,now +'...Extension Error ######## You chose a wrong file!! Please choose a fitsfile'+'\n','warning')
                 
    return None
    #if enable==1:
    #    menu.entryconfigure(2, state=tk.NORMAL)
    #else:
     #   menu.entryconfigure(2, state=tk.DISABLED)


def enabling():
    #print('enabling')
    menu.entryconfigure("Plot", state=tk.NORMAL)




def readnplot():
    
    def get_hours(time_str):
        h, m, s = time_str.split(':')
        return int(h) + int(m)/ 60.0 + int(s)/3600.0
    #gloss_FITS_Filename="/home/gloss/GLOSS_fitsData/2018/1811/GLOSS_20181109T023020To20181109T114305.fits"
    #filename = filedialog.askopenfilename(title='Select the input fits file')
    #filename=fitsname
    #filename=
    
        
    data =fits.open(file.name)
    hdr=data[0].header
    
    #Output_image_folder=os.getcwd()
    
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
    global fig
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
    
    
    r=plt.xlabel('Time(UT)',fontsize="10",color='white')
    
    plt.ylabel('Frequency(MHz)',fontsize="10",color='white')
    #fig.savefig(Output_image_folder+'GBD_DSPEC_'+str(hdr[8]).replace('/',"") +'.jpeg', dpi=100, facecolor='black', edgecolor='w',orientation='landscape', papertype=None, format=None,transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None)  
    #plt.show()
    plt.gcf().canvas.draw()
    global canvas
    canvas = FigureCanvasTkAgg(fig, master=gloss)  # A tk.DrawingArea.
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.X, expand=True)
    global toolbar
    toolbar = NavigationToolbar2Tk(canvas, mainFrame)
    toolbar.update()
    mainFrame.pack()
    #canvas.get_tk_widget().pack(side=tk.BOTTOM,anchor=tk.W,fill=tk.BOTH, expand=True)
    time()
    text.insert(tk.END,now +'...GLOSS spectrum plotted'+'\n')


def clear():

    canvas.get_tk_widget().pack_forget()
    mainFrame.pack_forget()
    time()
    text.insert(tk.END,now +'...cleared existing plot'+'\n')

    
    
def save():
    f=filedialog.asksaveasfilename(initialdir = 'pwd',title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    text2save = str(text.get(1.0, END)) # starts from `1.0`, not `0.0`
    f.write(text2save)
    time()
    text.insert(tk.END,now +'...saving the file in working directory'+'\n')
    f.close() 
    
        
def close():
        if messagebox.askyesno("GLOSS plotting tool","Do you really want to exit?"):
            time()
            text.insert(tk.END,now +'...killing the window'+'\n')
            gloss.destroy()
        else:
            return None


def analysis():
    
    
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
            # drawnow();
            # plot(abs(p),'r','linewidth',3);
            # grid;
            if (count > 300):
                print('Error...! Solution not converging !!! \\n')       # printing the error message
                break
            
        
                # plot(abs(p));
                # if (count < 300)
                #     fprintf('The solution = ');  %printing the result
                #     x
                #     fprintf('\nNumber of itration taken = %d\n',count);
                # end
                
        return r0



    
    
    #clc()
    #clear(mstring('all'))
    #format(mstring('long'))
    #syms(mstring('r'))
    #fp=input('Enter the plasma frequency : ')%fp=70;                        % plasma frequency
        
    r=Symbol('r')
    
    fp = [80, 55]
    #fp=linspace(50,120,100);
    #fp=[51,54,58,61,64,68,71,74,78,81,84,88,113,140,143,147,150,153,160,164,167,170,174,177,183,189,191,194,197,201,204,208,211,214,218,221,224,228,231,241,244,248,251,254,258,260,264,267,270,277,290,300,304,307,310,314,324,327,330,334,337,340,344,347,357,367,370,374,380,384,387,390,394,397,400,404]
    Ne = (np.divide(fp,(9 * 10 ** (-3))))** 2# electron density from plasma frequency
    
    #fp = sqrt(Ne(i))*(9*10^(-3)) 
    
    Sc = tk.simpledialog.askinteger('Enter factor','Enter Enhancement Factor')
    
    Hsa=np.zeros(len(fp))
    Hnk=np.zeros(len(fp))
    Hba1=np.zeros(len(fp))
    Hlb=np.zeros(len(fp))
    Hhy=np.zeros(len(fp))
    Hsgc=np.zeros(len(fp))
    Hba2=np.zeros(len(fp))
    itr=np.zeros(9)
    Hsgp=np.zeros(len(fp))
    
    
    
    
    for i in range(len(fp)):
        #Newkirk model
    
    
        Hnk[i]=4.321/np.log10(Ne[i]/(Sc*4.2*10**4))
       
       
        # Baumbach-Allen model 1
    
        BA1 = (1.55*(r** (-6)) + 2.99 * (r ** (-16)) + 0.036 * (r ** (-1.5))) - (Ne[i] / Sc * (10 ** (-8)))
        #BA1=(10^14)*(1.55*r^(-6)+2.99*r^(-16))-(Ne(i)*10^6) ;
    
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
        #newton_raphson(SA)[1]= itr[3]
    
        # Hybrid model
    
        HYB = 15.45 * r ** (-16) + 3.16 * r ** (-6) + r ** (-4) + 0.0033 * r ** (-2) - Ne[i] / (Sc * (10 ** 8))
        Hhy[i]=newton_raphson(HYB)
        #newton_raphson(HYB)[1]= itr[4] 
    
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
        #,newton_raphson(SGp)[1]=, itr[5]
    
        Hsgc[i]=newton_raphson(SGc)
        #,newton_raphson(SGc)[1]= , itr[6] 
    
        #Leblanc
    
        LB = 3.3e5 * r ** -2 + 4.1e6 * r ** -4 + 8e7 * r ** -6 - Ne[i] / Sc
    
        Hlb[i]=newton_raphson(LB)
        #,newton_raphson(BA1)[1]=, itr[7] 
    
    
        # A1=2.855e8;
        # P1=7.142;
        # anshu=A1*r^(-P1) - Ne(i)/Sc ;
        # [anshu_h(i),itr7]=newton_raphson(anshu) ;
    
    
        # A1 = 1.521*10^8 ; a1 = 7.279 ;
        # A2 = 1.84*10^8 ; a2 = 7.938 ;
        # A3 = 2.07*10^7 ; a3 = 4.852 ;
        # A4 = 7.52*10^5 ; a4 = 2.024 ;
        # 
        # 
        # WANG = A1*r^(-a1) + A2*r^(-a2) + A3*r^(-a3) + A4*r^(-a4)- Ne(i)/Sc ;
        # 
        # [Hwang(i),itr8]=newton_raphson(WANG) ;
    
    
        #heights=vertcat(Hnk,Hba1,Hba2,Hsa,Hsgp,Hsgc,Hhy,Hlb)';
        #heights = 
        
        
        
    heights=(Hnk, Hba1,Hba2, Hsa, Hsgc, Hhy, Hlb)    
    text.insert(tk.END,"NewkirkModel height values="+str(heights[0])+"\n")
    text.insert(tk.END,"Baumbach-Allen Model1 height values="+str(heights[1])+"\n")
    text.insert(tk.END,"Baumbach-Allen Model2 height values="+str(heights[2])+"\n")

    text.insert(tk.END,"Saito Model height values="+str(heights[3])+"\n")
    text.insert(tk.END,"Sittler & Guhathakurta model height values="+str(heights[4])+"\n")
    text.insert(tk.END,"Hybrid model height values="+str(heights[5])+"\n")
    text.insert(tk.END,"Leblanc Model height values="+str(heights[6])+"\n")
    
    












gloss= tk.Tk()
gloss.style = Style()
gloss.style.theme_use("clam")
gloss.minsize(1200,600)
gloss.title('GLOSS plotting tool')
scrollbar = Scrollbar(gloss)
scrollbar.pack(side=RIGHT, fill=Y)
#gloss.configure(background="white")
mainFrame = tk.Frame(gloss)
mainFrame.pack()
global text
text=Text(gloss,width=55,height=60,fg="black",yscrollcommand=scrollbar.set)
text.pack(side=tk.RIGHT,anchor=tk.SE, expand=False)
text.tag_config('warning', foreground="red")



menu = tk.Menu(gloss)

mfile = tk.Menu(menu)
mfile.add_command(label="Open", command=OpenFile)
mfile.add_command(label="Save", command=save)

mplot = tk.Menu(menu)
mplot.add_command(label="Plot spectrum", command=readnplot)

manalys=tk.Menu(menu)
manalys.add_command(label="start",command=analysis)

menu.add_cascade(label='File', menu=mfile)
menu.add_cascade(label='Plot', menu=mplot, state=tk.DISABLED)
menu.add_cascade(label='Analysis',menu=manalys)
menu.add_command(label="Clear", command=clear)
menu.add_command(label="Close", command=close)
menu.configure(background="black",foreground='white')    
gloss.config(menu=menu)
gloss.mainloop()





