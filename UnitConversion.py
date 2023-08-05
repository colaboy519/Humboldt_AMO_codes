import math
import numpy as np
from sys import argv
import sys
#from varname import varname

hplanck=6.62607015*10**(-34)
eelectron=1.60217663*10**(-19)
clight=299792458
pi=math.pi


def au2nm(e_au):
    e_ev=27.211396*e_au
    e_J=e_ev*eelectron
    e_Hz=e_J/hplanck
    e_nm=clight/e_Hz*10**9
    return e_nm

def nm2au(e_nm):
    e_Hz=clight/e_nm*10**9
    e_J=e_Hz*hplanck
    e_ev=e_J/eelectron
    e_au=e_ev/27.211396
    return e_au

def fs2au(t_fs):
   t_au=t_fs/0.02418884254
   return t_au

def au2fs(t_au):
   t_fs=t_au*0.02418884254
   return t_fs

def sign(n): 
    if n<0: return -1 
    elif n>0: return 1 
    else: 0 

def str2val(teststr):
    posE=teststr.find('E')
    posD=teststr.find('D')
    posED=max([posE,posD])
    if posED !=-1:
        base=float(teststr[0:posED])
        exp=int(teststr[posED+1:])
        val=base*10**exp
    if posED ==-1:
        val=float(teststr)

    return val




def cart2sphr(z_,r_,R_int_):
        x=(math.sqrt(r_**2+(z_+R_int_/2)**2)+math.sqrt(r_**2+(z_-R_int_/2)**2))/R_int_
        y=(math.sqrt(r_**2+(z_+R_int_/2)**2)-math.sqrt(r_**2+(z_-R_int_/2)**2))/R_int_
        return x,y



def MatrixListTranspose(M_):
    EnL=[]
    for column in range(len(M_[0])):
        EnC=[]
        for line in range(len(M_)):
            EnC.append(M_[line][column])
        EnL.append(EnC)
    return EnL
#print (cart2sphr(-1,-1,1))

def p_str(R_int_,Ndig_):
    '''
    '''
    R_nbr=int(R_int_)
    R_dec=int(round((R_int_-R_nbr)*(10**Ndig_),Ndig_))
    R_nbr_str="% s" % R_nbr
    R_dec_str = "% s" % R_dec
    if len(R_dec_str)>Ndig_:
        R_dec_str=R_dec_str[:Ndig_]
    strfill=''.join(['0']*(Ndig_-len(R_dec_str)))
    return ''.join((R_nbr_str,'p',strfill,R_dec_str))

   
def sum_laser(w,t,int1,int2):

    if w=='scan':summary='freqscan'+'_t'+str(t)+'_i'+p_str(int1,3)+'E'+str(int2)
    elif t=='scan': summary='durationscan'+'_w'+p_str(w,3)+'_i'+p_str(int1,3)+'E'+str(int2)
    elif int1=='scan':   summary='intscan'+'_w'+p_str(w,3)+'_t'+str(t)+'_i'+'X'+'E'+str(int2)
    elif int2=='scan': summary='intscan'+'_w'+p_str(w,3)+'_t'+str(t)+'_i'+p_str(int1,3)+'EX'
    else:  summary='w'+p_str(w,3)+'_t'+str(t)+'_i'+p_str(int1,3)+'E'+str(int2)


    

    return summary




colorsN=['cadetblue','yellowgreen','darkmagenta','silver','gold','darkblue','yellowgreen','silver','darkblue','yellowgreen','silver']

#https://matplotlib.org/2.0.2/examples/color/named_colors.html
# https://matplotlib.org/stable/tutorials/colors/colormaps.html
#print(fs2au(28.26))


def namestr(obj):
    return print( [name for name in globals() if globals()[name] is obj],obj)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def closest(lst, K):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

def Wqcm2au(Wqcm):
    au=Wqcm/(3.509*10**(16))
    return au



def bq(compl_):
    return abs(compl_*np.conj(compl_))


def val2zerostr(n,d):
    str0=''.join(['0']*d)
    if n==0:
        strnbr=''.join(['0']*d) 
    else:
        n_ziffern=math.floor(math.log10(n)+1)
        if n_ziffern>d:
            print('ERROR: Problem in fct val2zerostr')
        str0=''.join(['0']*(d-n_ziffern))   
        strnbr=str0+str(n)
    return strnbr


def dashedstring(str_):
    p1=str_.find('-')
    if p1==-1:
        res=[float(str_)]
    else:
        p2=p1+str_[p1+1:].find('-')+1
        val1=float(str_[:p1])
        val2=float(str_[p1+1:p2])
        val3=float(str_[p2+1:])
        if val3<=val2:
            print('ERROR IN SCAN INPUT')
            sys.exit()
        res=np.arange(val1, val3+val2, val2)
    return res


