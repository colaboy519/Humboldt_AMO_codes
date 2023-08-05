#   cd /users/amo/walterni/MyGeneralScripts/ && python LaserExcitation.py


import math
import numpy as np
from scipy.fft import fft, ifft, ifftshift
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import sys 


from UnitConversion import *
Npoints=20000
sf=10

plot_pulse=0
if plot_pulse==1:
    fig, axs = plt.subplots(4,1)



def LaserFT(omega_,Tau_,intensity_):  # all inputs in a.u.

    E0=intensity_**0.5
    #print('EO peak in a.u.:', E0)
    
    timedomain_x=np.linspace(-Tau_,Tau_,Npoints)
    timedomain_y=np.linspace(-Tau_,Tau_,Npoints)

   

    for t in range(len(timedomain_x)):
        oscill=math.cos(omega_*timedomain_x[t])*E0
        ampl=math.cos((math.pi/(2*Tau_))*timedomain_x[t])
        timedomain_y[t]=oscill


        timedomain_y[t]=ampl**2*oscill
        #timedomain_y[t]=ampl**2*oscill

    if plot_pulse==1:
        axs[0].plot(timedomain_x,timedomain_y)



    timedomain_y_ext=np.concatenate(([0]*Npoints*sf,timedomain_y,[0]*Npoints*sf))

    #timedomain_y_ext=timedomain_y
    freqdomain_y=fft(ifftshift(timedomain_y_ext))

    N_true=round(len(freqdomain_y)/2)
    freqdomain_y_1=freqdomain_y[:N_true]


    freq_max=1/(2*Tau_)*math.pi*Npoints
    freqdomain_x_1=np.linspace(0,1,N_true)*freq_max


    if plot_pulse==1:
        axs[1].plot(np.real(freqdomain_y),color='r')
        axs[1].plot(np.imag(freqdomain_y),color='b')
        axs[1].plot(np.abs(freqdomain_y),color='k')
  


        axs[2].plot(freqdomain_x_1,np.real(freqdomain_y_1),color='r')
        axs[2].plot(freqdomain_x_1,np.imag(freqdomain_y_1),color='b')
        axs[2].plot(freqdomain_x_1,np.abs(freqdomain_y_1),color='k')
        axs[2].set_xlim([0,3*omega_])  


    freqdomain_x_2=np.concatenate(([-100,-10**(-10)],freqdomain_x_1,[100]))
    freqdomain_y_2=np.concatenate(([0,0],freqdomain_y_1,[0]))


 

    

    #freqdomain_y_1=np.abs(np.fft.fftshift( np.fft.fft(timedomain_y_ext) )) 

   

    intfct_r=interp1d(freqdomain_x_2,np.real(freqdomain_y_2))
    intfct_i=interp1d(freqdomain_x_2,np.imag(freqdomain_y_2))
    intfct_a=interp1d(freqdomain_x_2,np.abs(freqdomain_y_2))

    if plot_pulse==1:
        xx=np.linspace(0,3*omega_,10000)
        axs[3].plot(xx,intfct_r(xx),color='r')
        axs[3].plot(xx,intfct_i(xx),color='b')
        axs[3].plot(xx,intfct_a(xx),color='k')

        plt.show()
       

   

    return intfct_r

if plot_pulse==1:
    intz=LaserFT(0.5,200,Wqcm2au(10**11))




  ## fourier Transform
    #


    #freqdomain_x_1 =np.fft.fftshift( np.fft.fftfreq(timedomain_y_ext.shape[0],timedomain_x[1]-timedomain_x[0]) )

    
    #fmax=1/(2*Tau_)*math.pi*Npoints
    #freqdomain_x=np.linspace(0,1,round(len(freqdomain_y_1)/2))*fmax
    #freqdomain_y=freqdomain_y_1[round(len(freqdomain_y_1)/2):]
    

    #freqdomain_x_ext=np.concatenate(([-100],freqdomain_x,[+100]))
    #freqdomain_y_ext=np.concatenate(([0],freqdomain_y,[0]))

    #intfct=interp1d(freqdomain_x_ext,freqdomain_y_ext)
    
    
    #aa=freqdomain_y[0]/max(freqdomain_y)
    #bb=freqdomain_y[-1]/max(freqdomain_y)

    #print(freqdomain_x[0],freqdomain_x[-1])
    #print(freqdomain_y[0],freqdomain_y[-1])
    #print(aa,bb)

    
        
        #axs[1].plot(freqdomain_x_1,freqdomain_y_1)
        #axs[2].plot(freqdomain_x,freqdomain_y)







# timexy,timeyy,freqxy,freqyy,inty=LaserFT(1,200,1,'y')
# plt.plot(freqxy,freqyy,color='red')

# plt.show()














