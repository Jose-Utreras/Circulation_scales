import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.nddata.utils import block_reduce
import os , sys, glob
from common_functions import colorplot
from scipy.special import erf
import subprocess
from scipy.optimize import curve_fit
from subprocess import check_output
import time
from scipy.interpolate import interp1d
from astropy.table import Table , Column ,vstack,hstack
import linecache

def read_line(name,number):
    linea = linecache.getline(name, number)
    linea = linea.split('\t')
    for i, ele in enumerate(linea):
        linea[i]=float(ele)
    return linea[0],linea[1],linea[2],np.array(linea[3:])

def sigma_array(name,n1,n2,pc,vo,steps, n1_min, n1_max, n2_min, n2_max, p_min, p_max):
    lista = sorted(glob.glob('Samplers/'+name_file+'_sampler*'))
    n_procs=len(lista)
    nsamp = steps/n_procs

    aux= int(1.0*(steps-1)*1.0*(n1-n1_min)/(n1_max-n1_min))
    number_file = int(aux/nsamp)

    i1 = int((aux%nsamp)*steps*steps)
    i2 = int((steps-1)*1.0*(n2-n2_min)/(n2_max-n2_min))*steps
    i3 = int((steps-1)*1.0*(np.log(1.0*pc/p_min)/np.log(1.0*p_max/p_min)))

    n11, n21, pc1, s1 = read_line(lista[number_file],i1+i2+i3)
    n12, n22, pc2, s2 = read_line(lista[number_file],i1+i2+i3+1)
    n13, n23, pc3, s3 = read_line(lista[number_file],i1+i2+i3+steps)
    n14, n24, pc4, s4 = read_line(lista[number_file],i1+i2+i3+steps+1)

    number_file = int((aux+1)/nsamp)
    i1 = int(((aux+1)%nsamp)*steps*steps)

    n15, n25, pc5, s5 = read_line(lista[number_file],i1+i2+i3)
    n16, n26, pc6, s6 = read_line(lista[number_file],i1+i2+i3+1)
    n17, n27, pc7, s7 = read_line(lista[number_file],i1+i2+i3+steps)
    n18, n28, pc8, s8 = read_line(lista[number_file],i1+i2+i3+steps+1)

    v1=1.0/np.abs((n1-n11)*(n2-n21)*(pc-pc1))
    v2=1.0/np.abs((n1-n12)*(n2-n22)*(pc-pc2))
    v3=1.0/np.abs((n1-n13)*(n2-n23)*(pc-pc3))
    v4=1.0/np.abs((n1-n14)*(n2-n24)*(pc-pc4))
    v5=1.0/np.abs((n1-n15)*(n2-n25)*(pc-pc5))
    v6=1.0/np.abs((n1-n16)*(n2-n26)*(pc-pc6))
    v7=1.0/np.abs((n1-n17)*(n2-n27)*(pc-pc7))
    v8=1.0/np.abs((n1-n18)*(n2-n28)*(pc-pc8))
    vtot=v1+v2+v3+v4+v5+v6+v7+v8
    v1/=vtot
    v2/=vtot
    v3/=vtot
    v4/=vtot
    v5/=vtot
    v6/=vtot
    v7/=vtot
    v8/=vtot

    sigma = v1*s1+v2*s2+v3*s3+v4*s4+v5*s5+v6*s6+v7*s7+v8*s8
    return vo*sigma

###############################################################################

name_file=str(np.loadtxt('Input.txt',usecols=(1,),dtype=bytes).astype(str))

L , N, n1_min, n1_max, n2_min, n2_max, dv_min, dv_max, nsteps  = np.loadtxt('Input.txt',usecols=(2,3,6,7,8,9,10,11,12),unpack=True)

scales      = np.loadtxt('Files/'+name_file+'_scales.txt', dtype=int)
N_scales    = len(scales)
N           = int(N)
DX          = L/N
DA          = DX**2
p_min       = 4.0/L
p_max       = 1.0/(4.0*DX)


lista = sorted(glob.glob('Chains/'+name_file+'/*chain'))


rmax_array  = []
rmin_array  = []

for li in lista:
    elements=li.split('_')

    rmax_array.append(float(elements[-2]))
    rmin_array.append(float(elements[-3]))

Profile     = []
Profile_up  = []
Profile_lw  = []
Scales      = []



for k, rmin , rmax, li in zip(range(len(rmax_array)),rmin_array,rmax_array,lista):
    print(rmin,rmax)
    scales_str = []
    scales_flt = []
    resol      = ''

    tab=Table.read(li,path='data')

    n1_arr=tab['n1']
    n2_arr=tab['n2']
    kc_arr=tab['pc']
    vo_arr=tab['vo']

    walker_output=[]
    counter=0

    n1_16,n1_50,n1_84=np.percentile(n1_arr,[16,50,84])
    n2_16,n2_50,n2_84=np.percentile(n2_arr,[16,50,84])
    kc_16,kc_50,kc_84=np.percentile(kc_arr,[16,50,84])
    vo_16,vo_50,vo_84=np.percentile(vo_arr,[16,50,84])

    args=[int(nsteps),n1_min, n1_max, n2_min, n2_max, p_min, p_max]

    print('Computing sigma')
    sigma = sigma_array(name_file,n1_50,n2_50,kc_50,vo_50,*args)
    print('Computing error')
    dn1   = 0.5*(sigma_array(name_file,n1_84,n2_50,kc_50,vo_50,*args)-sigma_array(name_file,n1_16,n2_50,kc_50,vo_50,*args))
    dn2   = 0.5*(sigma_array(name_file,n1_50,n2_84,kc_50,vo_50,*args)-sigma_array(name_file,n1_50,n2_16,kc_50,vo_50,*args))
    dkc   = 0.5*(sigma_array(name_file,n1_50,n2_50,kc_84,vo_50,*args)-sigma_array(name_file,n1_50,n2_50,kc_16,vo_50,*args))
    dvo   = 0.5*(vo_84-vo_16)*(sigma/vo_50)

    error = np.sqrt(dn1**2+dn2**2+dkc**2+dvo**2)

    res            = scales*DX

    p10            = np.zeros(N_scales,dtype=float)
    p20            = np.zeros(N_scales,dtype=float)
    p30            = np.zeros(N_scales,dtype=float)
    p40            = np.zeros(N_scales,dtype=float)
    p50            = np.zeros(N_scales,dtype=float)
    p60            = np.zeros(N_scales,dtype=float)
    p70            = np.zeros(N_scales,dtype=float)
    p80            = np.zeros(N_scales,dtype=float)
    p90            = np.zeros(N_scales,dtype=float)

    q10            = np.zeros(N_scales,dtype=float)
    q20            = np.zeros(N_scales,dtype=float)
    q30            = np.zeros(N_scales,dtype=float)
    q40            = np.zeros(N_scales,dtype=float)
    q50            = np.zeros(N_scales,dtype=float)
    q60            = np.zeros(N_scales,dtype=float)
    q70            = np.zeros(N_scales,dtype=float)
    q80            = np.zeros(N_scales,dtype=float)
    q90            = np.zeros(N_scales,dtype=float)

    radius  = 0.5*(rmin+rmax)

    print('Loading Circulation')
    Tabla       = np.loadtxt('Files/'+name_file+'_vort.txt')
    Radius      = Tabla[:,1]
    ring        = (Radius>=rmin)&(Radius<=rmax)
    Tabla       = Tabla[ring]

    for i in range(N_scales):
        X    = Tabla[:,2+2*i]/DA
        Y    = Tabla[:,3+2*i]/DA
        Z    = Y + np.random.normal(0,sigma[i],size=len(Y))/DA

        p10[i], p20[i], p30[i], p40[i], p50[i], p60[i], p70[i], p80[i], p90[i] = np.percentile(X,[10,20,30,40,50,60,70,80,90])
        q10[i], q20[i], q30[i], q40[i], q50[i], q60[i], q70[i], q80[i], q90[i] = np.percentile(Z,[10,20,30,40,50,60,70,80,90])

    psum=p10+p20+p30+p40+p50+p60+p70+p80+p90
    qsum=q10+q20+q30+q40+q50+q60+q70+q80+q90

    delta=(psum-qsum)*DA/9

    Gamma_turb     = np.zeros(N_scales,dtype=float)
    Gamma_up       = np.zeros(N_scales,dtype=float)
    Gamma_lw       = np.zeros(N_scales,dtype=float)

    for i in range(N_scales):
        X    = Tabla[:,2+2*i]/DA
        Y    = Tabla[:,3+2*i]/DA
        Z    = np.random.normal(delta[i],sigma[i],size=len(Y))/DA

        Up      = np.random.normal(delta[i],sigma[i]+error[i],size=len(Y))/DA
        Lw      = np.random.normal(delta[i],max(sigma[i]-error[i],0),size=len(Y))/DA

        Gamma_turb[i]  = (np.abs(Z)**2).sum()/(np.abs(Y)**2).sum()
        Gamma_up[i]    = (np.abs(Up)**2).sum()/(np.abs(Y)**2).sum()
        Gamma_lw[i]    = (np.abs(Lw)**2).sum()/(np.abs(Y)**2).sum()

        p10[i], p20[i], p30[i], p40[i], p50[i], p60[i], p70[i], p80[i], p90[i] = np.percentile(X,[10,20,30,40,50,60,70,80,90])
        q10[i], q20[i], q30[i], q40[i], q50[i], q60[i], q70[i], q80[i], q90[i] = np.percentile(Y+Z,[10,20,30,40,50,60,70,80,90])

    fun=interp1d(Gamma_turb,res)
    Scales.append(radius)
    try:
        Profile.append(fun(1.0))
    except:
        Profile.append(res.min())

    fun=interp1d(Gamma_up,res)
    try:
        Profile_up.append(fun(1.0))
    except:
        Profile_up.append(2*res.min())

    fun=interp1d(Gamma_lw,res)
    try:
        Profile_lw.append(fun(1.0))
    except:
        Profile_lw.append(0.0)

    print(Profile)

    f1         = plt.figure(figsize=(10,8))
    ax1        = f1.add_subplot(111)

    ax1.plot(res, Gamma_turb , color=colorplot(11,1), label='model',linestyle='-.')
    ax1.fill_between(res, np.minimum(Gamma_lw,Gamma_up),np.maximum(Gamma_lw,Gamma_up), color=colorplot(11,1), label='',alpha=0.2)
    ax1.set_xlabel(r'Scale [pc]',size=20)
    ax1.set_xscale('log')
    ax1.legend(loc='upper right', prop={'size':15}, ncol = 2)
    f1.savefig('Plots/'+name_file+'_%05d_' %int(radius)+'fraction.png')
    plt.close(f1)


    tabla=Table()
    tabla['resolution']   = Column(np.array(res)       , description='resolution in pc')
    tabla['fraction']     = Column(np.array(Gamma_turb), description='fraction')
    tabla['fraction_low'] = Column(np.array(Gamma_lw)  , description='lower limit')
    tabla['fraction_upp'] = Column(np.array(Gamma_up)  , description='upper limit')
    tabla.write('Tables/'+name_file+'_%05d_'%int(radius)+'fraction',path='data',format='hdf5',overwrite=True)
    del tabla

    plt.plot(res,p10,linestyle='-',color='black')
    plt.plot(res,p20,linestyle='-',color='black')
    plt.plot(res,p30,linestyle='-',color='black')
    plt.plot(res,p40,linestyle='-',color='black')
    plt.plot(res,p50,linestyle='-',color='black')
    plt.plot(res,p60,linestyle='-',color='black')
    plt.plot(res,p70,linestyle='-',color='black')
    plt.plot(res,p80,linestyle='-',color='black')
    plt.plot(res,p90,linestyle='-',color='black')

    plt.plot(res,q10,linestyle='--',color='black')
    plt.plot(res,q20,linestyle='--',color='black')
    plt.plot(res,q30,linestyle='--',color='black')
    plt.plot(res,q40,linestyle='--',color='black')
    plt.plot(res,q50,linestyle='--',color='black')
    plt.plot(res,q60,linestyle='--',color='black')
    plt.plot(res,q70,linestyle='--',color='black')
    plt.plot(res,q80,linestyle='--',color='black')
    plt.plot(res,q90,linestyle='--',color='black')

    plt.savefig('Plots/'+name_file+'_%05d_' %int(radius)+'_percentiles.png')
    plt.close()



Scales      = np.array(Scales)
Profile     = np.array(Profile)
Profile_lw  = np.array(Profile_lw)
Profile_up  = np.array(Profile_up)

plt.figure(figsize=(10,8))
plt.plot(Scales,Profile)
plt.fill_between(Scales,Profile_lw,Profile_up,alpha=0.2)
plt.savefig('Plots/'+name_file+'_profile.png')

tabla=Table()
tabla['radius']    = Column(np.array(Scales)      , description='radius in pc')
tabla['scale']     = Column(np.array(Profile)     , description='transition scale')
tabla['scale_low'] = Column(np.array(Profile_lw)  , description='lower limit')
tabla['scale_upp'] = Column(np.array(Profile_up)  , description='upper limit')
tabla.write('Tables/'+name_file+'_profile' ,path='data',format='hdf5',overwrite=True)
