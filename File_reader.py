from scipy.io import FortranFile
import numpy as np
from UnitConversion import *




def dip_reader(filename_):
    path_dip_out='/users/amo/walterni/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/DIPOLE_2e/out/H_H/H_H_'
    f_dip= FortranFile(path_dip_out+filename_+'/SpS_s_gu.'+Gauge,'r')
    N1=f_dip.read_ints(np.int32)[0]
    N2=f_dip.read_ints(np.int32)[0]
    E1=f_dip.read_reals(float)
    E2=f_dip.read_reals(float)
    M=[0]*N1
    for line in range (N1):
        cont=f_dip.read_reals(float)
        M[line]=cont
    return E1,E2,M


def Nucfix_base_reader(filename_):
    path='/users/amo/walterni/nucfix/input/basis/h2b'+filename_+'.inp'
    f=open(path, 'r')
    kk=f.read(9999)
    f.close()
    pos=kk.find('box radius')
    R_max=float(str2val(kk[pos+20:pos+27]))
    pos=kk.find('Order')
    B_k=int(kk[pos+30:pos+32])
    pos=kk.find('Number')
    B_N=int(kk[pos+30:pos+34])
    
    B_s=B_N-B_k+1         # number intervalls
    r_intlen=R_max/(B_s)  # length intervall
    startknots=[0]*B_k
    endknots=[R_max]*B_k
    middleknots=np.arange(r_intlen,R_max,r_intlen)
    B_knots=np.concatenate((startknots,middleknots, endknots))
    
    return B_N, B_k,R_max,B_knots


 # freq, length, Int (rest see template)

def lines_reader(path_,line_numbers_):  #first line=1
    for i in range (len(line_numbers_)):
        line_numbers_[i]=line_numbers_[i]-1
    lines = []
    with open(path_, 'r') as fp:
        for i, line in enumerate(fp):
            # read line 4 and 7
            if i in line_numbers_:
                lines.append(line.strip())
            elif i > max(line_numbers_):
                # don't read after line 7 to save time
                break
    return lines

def valselect(str_):
    return str_[str_.find(':')+1:]
  

def pulse_reader(pulsename_):
    path='/users/amo/walterni/AMO_TOOLS_USER/TIME_PROP/PULSE_INPUTS/'+pulsename_+'.pls'
    res=lines_reader(path,[25,47,76])
    freq=float(valselect(res[0]))
    length=int(valselect(res[1]))
    int1=float(res[2][res[2].find(':')+1:res[2].find('D+')])
    int2=int(res[2][res[2].find('D+')+1:])
    #print(pulsename_,'w =',freq, ', t =',length,', int =',int1,'D+',int2)
    return freq, length, int1, int2
   

def laserfilename_creator(freq,length,int1,int2):
    gb='w'+p_str(freq,3)+'_t'+str(length)+'_i'+p_str(int1,3)+'E'+str(int2)
    return gb







def Nucfix_reader(filename1_,filename2_,filename3_):
    for filename in (filename1_,filename2_):
        f_ew= FortranFile('/users/amo/walterni/nucfix/scr/'+filename+'.ewb', 'r')
        N_vib=f_ew.read_ints(np.int32)[0]
        EnT=[0]*N_vib
        for line in range(N_vib):
            cont=f_ew.read_reals(float)[0]
            EnT[line]=cont
        f_ew.close()
        f_ev= FortranFile('/users/amo/walterni/nucfix/scr/'+filename+'.evb', 'r')
        k1,k2=f_ev.read_ints(np.int32)
        Ev=f_ev.read_reals(float)
        EvM=np.reshape(Ev, (N_vib,N_vib))  # was B_N-2 instead of N_vib before
        f_ev.close

        if filename==filename1_:
            E1,EV1=EnT,EvM
        if filename==filename2_:
            E2,EV2=EnT,EvM


    f_M= FortranFile('/users/amo/walterni/nucfix/out/'+filename3_+'.cnm.out', 'r')
    M=[0]*len(E2)
    for line in range(len(E2)):
        cont=f_M.read_reals(float)[0]
        M[line]=cont
    f_M.close

    
    
    
    ### M is Matrix element from first state in ground state in file1 to all other states in file 2
    return E1,E2,EV1,EV2,M

def NF_pot_reader(statename_):
    path_pot='/users/amo/walterni/nucfix/input/potentials/h2st'+statename_+'p001.pot'
    f_pot= open(path_pot, 'r')
    kk=f_pot.read(9999)
    f_pot.close
    pos=kk.find('Number of grid points')
    N_grid_pot=int(kk[pos+51:pos+55])
    grT=[]
    pointer=pos+385
    symbs=(' ','\n','\t','*')
    teststr=20
    counter=0
    while counter <2*N_grid_pot: 
        posref=teststr
        for symb in symbs:
            pos=kk[pointer:pointer+teststr].find(symb)
            if 0<=pos and pos<posref:
                posref=pos
        if posref==0:
                pointer=pointer+1
        if posref>=1:
                cont=kk[pointer:pointer+posref]
                cont_z=str2val(cont)
                pointer=pointer+posref+1
                counter=counter+1
                grT.append(cont_z)
    grx=grT[::2]
    gry=grT[1::2]
    return grx,gry