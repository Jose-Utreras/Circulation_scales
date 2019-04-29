import numpy as np
import os , sys, glob
from common_functions import  gaussian, Smooth
import subprocess
from scipy.optimize import curve_fit
from subprocess import check_output
import time
from scipy.interpolate import interp1d

eta=2*0.6744897501
kms_to_pcmyr=1.022

cmdargs = sys.argv

av_file     =   cmdargs[-1]
Rmax        =   float(cmdargs[-2])
Rmin        =   float(cmdargs[-3])

name_file=str(np.loadtxt('Input.txt',usecols=(1),dtype=bytes).astype(str))

L , N, n1_min, n1_max, n2_min, n2_max, dv_min, dv_max  = np.loadtxt('Input.txt',usecols=(2,3,6,7,8,9,10,11),unpack=True)

scales      = np.loadtxt('Files/'+name_file+'_scales.txt', dtype=int)
N_scales    = len(scales)
N           = int(N)
DX          = L/N
DA          = DX**2
p_min       = 4.0/L
p_max       = 1.0/(4.0*DX)

resolutions = scales*DX
sigma_t     = 1.0*kms_to_pcmyr*2*DX # total error in circulation

stds=np.zeros(N_scales, dtype=float)
serr=np.zeros(N_scales, dtype=float)

aux_res     = np.insert(scales, 0, scales[0]/2.0)
aux_res     = np.insert(aux_res, len(aux_res), scales[-1]*2)
res_temp    = np.exp(np.linspace(np.log(0.75*scales.min()),np.log(1.5*scales.max()),1000))

ratio=[]

Tabla       = np.loadtxt('Files/'+name_file+'_vort.txt')
Radius      = Tabla[:,1]
ring        = (Radius>=Rmin)&(Radius<=Rmax)
Tabla       = Tabla[ring]

if (av_file=='model')|(av_file=='exp')|(av_file=='expansion'):
    Tabla2      = np.loadtxt('Files/'+name_file+'_other.txt')
    Radius      = Tabla2[:,1]
    ring        = (Radius>=Rmin)&(Radius<=Rmax)
    Tabla2      = Tabla2[ring]

    del Radius

    ID          = Tabla[:,0]
    Tabla=Tabla[ID.argsort()]
    del ID

    ID          = Tabla2[:,0]
    Tabla2=Tabla2[ID.argsort()]
    del ID


for i, R in enumerate(scales):
    print(R)
    c_sig               = sigma_t/R**1.5

    vort                = Tabla[:,2+2*i]

    if (av_file=='model'):
        omeg            = Tabla2[:,2+2*i]
        print('MODEL')
    elif(av_file=='exp')|(av_file=='expansion'):
        omeg            = Tabla2[:,3+2*i]
        print('EXPANSION')
    else:
        omeg            = Tabla[:,3+2*i]
        print('ROT-CURVE')

    mapa    = vort-omeg
    mapa   -= np.mean(mapa)
    try:
        c_len=len(mapa)
        pu,pl=np.percentile(mapa,[75,25])
        c_width=16*(pu-pl)*c_len**(-0.333)
        c_bins=1.0/c_width

        c_max=mapa.max()+0.5*c_width
        c_min=mapa.min()-0.5*c_width
        c_bins=max(int((c_max-c_min)*c_bins),4)
        c_width=1.0*(c_max-c_min)/(1.0*c_bins)

    except:
        c_len=len(mapa)
        c_max=mapa[0]+0.5*c_width
        c_min=mapa[0]-0.5*c_width
        c_bins=4
        c_width=1.0*(c_max-c_min)/(1.0*c_bins)
    if c_len<2:
        c_max=mapa[0]+2*c_sig
        c_min=mapa[0]-2*c_sig
        c_bins=4
        c_width=1.0*(c_max-c_min)/(1.0*c_bins)
    np.savetxt('Files/'+name_file+'_circ_vector.txt',mapa)
    parameters = name_file +' '+str(c_min)+' '+str(c_max)+' '+str(c_bins)+' '
    parameters+= str(c_width)+' '+str(c_sig)+' '+str(c_len)
    subprocess.call(['./histogram '+parameters], stdin = sys.stdin, shell=True)

    y,errors=np.loadtxt('Files/'+name_file+'_hist_data.txt',unpack=True)
    x=np.linspace(c_min,c_max,c_bins+1)
    x=0.5*(x[1:]+x[:-1])
    errors=np.sqrt(y+errors**2)
    errors[errors<1.0]=1.0

    std_lim=np.average((x-np.average(x,weights=y))**2,weights=y)**0.5
    a_min,a_max=np.average(y),y.max()*2
    xaux=np.insert(x,0,2*x[0]-x[1])
    yaux=np.insert(y,0,0)
    ysum=np.cumsum(yaux)/yaux.sum()
    fsum=interp1d(ysum,xaux)
    mean_lim=max(np.abs(fsum(0.6)),np.abs(fsum(0.4)))
    mu_min,mu_max=-mean_lim,mean_lim
    sig_lim=0.5*(fsum(0.84)-fsum(0.16))
    sig_min,sig_max=min(sig_lim,std_lim)/2.0,1.5*max(sig_lim,std_lim)

    try:
        popt, pcov   = curve_fit(gaussian,x,y,sigma=errors,bounds=([a_min,mu_min,sig_min],[a_max,mu_max,sig_max]))
        popt1,pcov1  = curve_fit(gaussian,x,y,sigma=errors,method='trf',bounds=([a_min,mu_min,sig_min],[a_max,mu_max,sig_max]))
        chi_0        = (y-gaussian(x,*popt))/errors
        chi_1        = (y-gaussian(x,*popt1))/errors
        chi_0        = chi_0**2
        chi_1        = chi_1**2
        chi_0=chi_0.sum()
        chi_1=chi_1.sum()
        if chi_0>chi_1:
            popt=popt1
            pcov=pcov1
    except:
        pass
    try:
        stds[i]=np.abs(popt[-1])
        serr[i]=np.sqrt(pcov.diagonal()[-1])
    except:
        stds[i]=c_sig
        serr[i]=c_sig


    outlier=False
    if i>2:
        aux_ratio = np.array(ratio)
        mu        = np.mean(aux_ratio)
        dev       = np.std(aux_ratio)
        QQ        = 0.5*(c_max-c_min)/stds[i]
        if (QQ<mu-5*dev)|(mu+5*dev<QQ):
            outlier=True
        del aux_ratio,QQ

    ratio.append(0.5*(c_max-c_min)/stds[i])
    if outlier:
        print('HAY UN OUTLIER!!!!')
        stds[i]=0.5*(c_max-c_min)/mu
        serr[i]=dev



if (serr[-1]>100.0)|(np.isnan(serr[-1])):
    serr[-1]=stds[-1]
infinitos=(serr>100.0)|(serr==np.inf)|(np.isnan(serr))
finitos=~infinitos
function=interp1d(scales[finitos],serr[finitos])
serr[infinitos]=function(scales[infinitos])

aux_stds=np.insert(stds,0,stds[0])
aux_stds=np.insert(aux_stds,len(aux_stds),stds[-1])

fun=interp1d(aux_res,aux_stds)
stds_temp=fun(res_temp)

smooth_stds=Smooth(np.log(res_temp),stds_temp,n=40,b=1)

fun_smooth=interp1d(res_temp,smooth_stds)

diff=np.abs(stds_temp-smooth_stds)
dev=np.std(diff[diff>0])
diff=np.abs(stds-fun_smooth(scales))
diff=np.sqrt(dev**2+diff**2)
stds=fun_smooth(scales)
stot=np.sqrt(serr**2+diff**2)
stot=np.maximum(stot,0.05*stds)

stds_string='stds'
stot_string='stot'

for ele1,ele2 in zip(stds,stot):
    stds_string+=','+str(ele1)
    stot_string+=','+str(ele2)
stds_string+='\n'
stot_string+='\n'

Resolution=''
sigma_array=''
error_array=''

print(stds)
print(stot)

for sa,re,se in zip(stds,resolutions,stot):
    sigma_array+=' %.03e' %sa
    Resolution +=' %.03e' %re
    error_array+=' %.03e' %se

print('Writing parameters')

array_file=name_file+'_%05d' %int(Rmin)+'_%05d' %int(Rmax)


linea= Resolution + sigma_array + error_array
variables=open('Files/'+name_file+'_%05d' %int(Rmin) +'_%05d' %int(Rmax)+'.txt','w')
variables.write(linea)
variables.close()
