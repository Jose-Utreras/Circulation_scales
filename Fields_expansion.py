import numpy as np
import math
from common_functions import *
import sys
import yt
from yt.units import kpc,pc,km,second,yr,Myr,Msun

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from yt.fields.api import ValidateParameter

def Smooth(x,y,n=1,b=1):
        NN=len(y[:])
        z=np.zeros(NN)
        d=np.zeros(NN)
        X=np.zeros(NN+2*n)
        Y=np.zeros(NN+2*n)


        for i in range(n):
                X[n-i-1]=x[0]-(i+1)*(x[1]-x[0])
                Y[n-i-1]=y[0]
        for i in np.arange(NN,NN+2*n,1):
                X[i]=x[NN-1]+(i+1-NN)*(x[1]-x[0])
                Y[i]=y[NN-1]
        count = n
        for xa,xb in zip(x,y):
                X[count]=xa
                Y[count]=xb
                count+=1
        for i in range(len(x)):
                for j in range(2*n+1):
                        z[i]=z[i]+np.exp( -(X[i+n]-X[i+j])**2 /(2*b**2) )*Y[i+j]
                        d[i]=d[i]+np.exp( -(X[i+n]-X[i+j])**2 /(2*b**2) )

        return z/d
def _vc(field,data):
        if data.has_field_parameter("bulk_velocity"):
                bv = data.get_field_parameter("bulk_velocity").in_units("cm/s")
        else:
                bv = data.ds.arr(np.zeros(3), "cm/s")
        xv = data["gas","velocity_x"] - bv[0]
        yv = data["gas","velocity_y"] - bv[1]
        center = data.get_field_parameter('center')
        x_hat = data["x"] - center[0]
        y_hat = data["y"] - center[1]
        r = np.sqrt(x_hat*x_hat+y_hat*y_hat)
        x_hat /= r
        y_hat /= r

        return (yv*x_hat-xv*y_hat)


yt.add_field("vc", function=_vc,take_log=False, units=r"km/s",validators=[ValidateParameter('bulk_velocity')])


def _radial_velocity(field,data):
    center = data.get_field_parameter('center')
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity").in_units("cm/s")
    else:
        bv = data.ds.arr(np.zeros(3), "cm/s")

    x  = data["x"] - center[0]
    y  = data["y"] - center[1]
    r  = np.sqrt(x**2+y**2)
    vx = data["gas","velocity_x"] - bv[0]
    vy = data["gas","velocity_y"] - bv[1]

    vr=vx*x+vy*y
    vr/=r

    return vr

yt.add_field("radial_velocity", function=_radial_velocity,validators=[ValidateParameter('center')], units=r"km/s")

def _Disk_H(field, data):
    center = data.get_field_parameter('center')
    z = data["z"] - center[2]
    return np.abs(z)

yt.add_field("Disk_H",
             function=_Disk_H,
             units="pc",
             take_log=False,
             validators=[ValidateParameter('center')])

def _Z(field, data):
    center = data.get_field_parameter('center')
    z = data["z"] - center[2]
    return z

yt.add_field("Z",
             function=_Z,
             units="pc",
             take_log=False,
             validators=[ValidateParameter('center')])

def _Disk_Radius(field, data):
    center = data.get_field_parameter('center')
    x = data["x"] - center[0]
    y = data["y"] - center[1]
    r = np.sqrt(x*x+y*y)
    return r
yt.add_field("Disk_Radius",
             function=_Disk_Radius,
             units="cm",
             take_log=False,
             validators=[ValidateParameter('center')])


def _Disk_Angle(field, data):
    center = data.get_field_parameter('center')
    x = data["x"] - center[0]
    y = data["y"] - center[1]
    r = np.arctan2(y,x)
    return r
yt.add_field("Disk_Angle",
             function=_Disk_Angle,
             units="dimensionless",
             take_log=False,
             validators=[ValidateParameter('center')])



def _Omega(field,data):
        if data.has_field_parameter("bulk_velocity"):
                bv = data.get_field_parameter("bulk_velocity").in_units("cm/s")
        else:
                bv = data.ds.arr(np.zeros(3), "cm/s")
        xv = data["gas","velocity_x"] - bv[0]
        yv = data["gas","velocity_y"] - bv[1]
        center = data.get_field_parameter('center')
        x_hat = data["x"] - center[0]
        y_hat = data["y"] - center[1]
        r = np.sqrt(x_hat*x_hat+y_hat*y_hat)
        x_hat /= r**2
        y_hat /= r**2

        return yv*x_hat-xv*y_hat


yt.add_field("Omega", function=_Omega,take_log=False, units=r"1/yr",validators=[ValidateParameter('bulk_velocity')])


def Sigma_profile(M,R,DR,Rmin,Rmax):

    Redges=np.arange(Rmin*pc,Rmax*pc,DR*pc)
    Rcen=0.5*(Redges[1:]+Redges[0:-1])
    N=len(Rcen)
    Scen=np.zeros(N)
    for k in range(N):
        ring=(Redges[k]<R)&(Redges[k+1]>R)
        if len(ring[ring])<4:
            Scen[k]=np.nan
        else:
            auxmas=M[ring]
            Scen[k]=np.sum(auxmas)

    correct=~np.isnan(Scen)
    Rcen=Rcen[correct]
    Scen=Scen[correct]
    Scen/=2*np.pi*Rcen*DR
    Rcen=np.insert(Rcen,0,0)
    Scen=np.insert(Scen,0,2*Scen[0]-Scen[1])

    return Rcen, Scen

def Velocity_profile(V,R,DR,Rmin,Rmax,function):

    Redges=np.arange(Rmin*pc,Rmax*pc,DR*pc)
    Rcen=0.5*(Redges[1:]+Redges[0:-1])
    N=len(Rcen)
    Vcen=np.zeros(N)
    for k in range(N):
        ring=(Redges[k]<R)&(Redges[k+1]>R)
        if len(ring[ring])<4:
            Vcen[k]=np.nan
        else:
            auxome=V[ring]
            Vcen[k]=function(auxome)

    correct=~np.isnan(Vcen)
    Rcen=Rcen[correct]
    Vcen=Vcen[correct]

    return Rcen,Vcen


def Vorticity_profile(V,R,DR,Rmin,Rmax,function):
    x,y = Velocity_profile(V,R,DR,Rmin,Rmax,function)
    U=np.copy(y)
    du=x[1:-1]*(U[2:]-U[:-2])/((x[2:]-x[:-2])*U[1:-1])
    du=np.insert(du,-1,du[-1])
    du=np.insert(du,0,du[0])
    y=y/x
    x=np.insert(x,0,0)
    y=np.insert(y,0,y[0]*2)
    y[1:]*=(1+du)

    return x,y

def add_extremes(X):
    X=np.insert(X,0,0)
    X=np.insert(X,len(X),2*X[-1]-X[-2])
    return X

def derivative(X,Y):
    dY=np.zeros_like(Y)
    dY[0]    = Y[1]-Y[0]
    dY[-1]   = Y[-1]-Y[-2]
    dY[1:-1] = 0.5*(Y[2:]-Y[:-2])

    dX=np.zeros_like(X)
    dX[0]    = X[1]-X[0]
    dX[-1]   = X[-1]-X[-2]
    dX[1:-1] = 0.5*(X[2:]-X[:-2])

    dY      /= dX
    return dY


def XY_field(N,L):

    X=np.array(list(np.reshape(range(N),(1,N)))*N)
    X=np.reshape(X,(N,N))-(N-1)/2
    X*=L/N

    return X,X.T

def polar_field(N,L):
    X,Y=XY_field(N,L)
    Radius=np.sqrt(X**2+Y**2)
    Theta=np.arctan2(Y,X)

    return Radius,Theta

def circular_dec(x,c1,c2,c3,c4,d0,d1,d2,d3,d4):
    result = d0+d1*np.cos(x)+d2*np.cos(2*x)+d3*np.cos(3*x)+d4*np.cos(4*x)

    result += c1*np.sin(x)+c2*np.sin(2*x)+c3*np.sin(3*x)+c4*np.sin(4*x)

    return result

def delete_nan(x):
    x[np.isnan(x)]=x[~np.isnan(x)][-1]

    return x

dir=str(np.loadtxt('Input.txt',usecols=(0),dtype=bytes).astype(str))
name_file=str(np.loadtxt('Input.txt',usecols=(1),dtype=bytes).astype(str))
L=np.loadtxt('Input.txt',usecols=(2),dtype=float)

directory=dir+'/'+name_file+'/G-'+name_file[-4:]

ds = yt.load(directory)

grids=ds.refine_by**ds.index.max_level*ds.domain_dimensions[0]

DX=ds.arr(1, 'code_length')
DX.convert_to_units('pc')
dr=DX/grids

NN=int(L/dr)
dA=((L/NN)**2)
print(L/NN)

width = (float(L), 'pc')

dd=ds.all_data()
disk_dd = dd.cut_region(["obj['Disk_H'].in_units('pc') < 1.0e3"])
proj = ds.proj('vc', 2,data_source=disk_dd,weight_field='density')

res = [NN, NN]
frb = proj.to_frb(width, res, center=[0.5,0.5,0.5])

radius  = frb['Disk_Radius'].in_units('pc')
angle   = frb['Disk_Angle']
vel_cir = frb['vc'].in_units('km/s')
vel_rad = frb['radial_velocity'].in_units('km/s')

Redges=np.arange(0*pc,35000*pc,500*pc)
Rcen=0.5*(Redges[1:]+Redges[0:-1])
N=len(Rcen)


B_0 = np.zeros(N)
B_1 = np.zeros(N)
B_2 = np.zeros(N)
B_3 = np.zeros(N)
B_4 = np.zeros(N)

A_1 = np.zeros(N)
A_2 = np.zeros(N)
A_3 = np.zeros(N)
A_4 = np.zeros(N)

D_0 = np.zeros(N)
D_1 = np.zeros(N)
D_2 = np.zeros(N)
D_3 = np.zeros(N)
D_4 = np.zeros(N)

C_1 = np.zeros(N)
C_2 = np.zeros(N)
C_3 = np.zeros(N)
C_4 = np.zeros(N)


for k in range(N):
    ring=(Redges[k]<radius)&(Redges[k+1]>radius)

    if len(ring[ring])<4:

        D_0[k]  = np.nan
        D_1[k]  = np.nan
        D_2[k]  = np.nan
        D_3[k]  = np.nan
        D_4[k]  = np.nan

        C_1[k]  = np.nan
        C_2[k]  = np.nan
        C_3[k]  = np.nan
        C_4[k]  = np.nan

        A_1[k]  = np.nan
        A_2[k]  = np.nan
        A_3[k]  = np.nan
        A_4[k]  = np.nan

        B_0[k]  = np.nan
        B_1[k]  = np.nan
        B_2[k]  = np.nan
        B_3[k]  = np.nan
        B_4[k]  = np.nan

    else:

        vc_ring=vel_cir[ring]
        vr_ring=vel_rad[ring]
        ang_ring=angle[ring]
        Nring=len(vc_ring)
        Nang=int(Nring/200)

        Aedges      = np.linspace(-np.pi,np.pi,Nang)
        Acen        = 0.5*(Aedges[1:]+Aedges[0:-1])

        p0=[np.std(vc_ring),0,0,0,np.mean(vc_ring),np.std(vc_ring),0,0,0]
        v99=np.percentile(vc_ring,99)
        popt, pcov = curve_fit(circular_dec,ang_ring,vc_ring,p0=p0,bounds=(-v99,v99))
        C_1[k],C_2[k],C_3[k],C_4[k],D_0[k],D_1[k],D_2[k],D_3[k],D_4[k] = popt

        p0=[np.std(vr_ring),0,0,0,np.mean(vr_ring),np.std(vr_ring),0,0,0]
        v99=max(np.abs(np.percentile(vr_ring,[1,99])))
        popt, pcov = curve_fit(circular_dec,ang_ring,vr_ring,p0=p0,bounds=(-v99,v99))
        A_1[k],A_2[k],A_3[k],A_4[k],B_0[k],B_1[k],B_2[k],B_3[k],B_4[k] = popt

A_1=delete_nan(A_1)
A_2=delete_nan(A_2)
A_3=delete_nan(A_3)
A_4=delete_nan(A_4)

B_0=delete_nan(B_0)
B_1=delete_nan(B_1)
B_2=delete_nan(B_2)
B_3=delete_nan(B_3)
B_4=delete_nan(B_4)

C_1=delete_nan(C_1)
C_2=delete_nan(C_2)
C_3=delete_nan(C_3)
C_4=delete_nan(C_4)

D_0=delete_nan(D_0)
D_1=delete_nan(D_1)
D_2=delete_nan(D_2)
D_3=delete_nan(D_3)
D_4=delete_nan(D_4)

Rcen=add_extremes(Rcen)
A_1=add_extremes(A_1)
A_2=add_extremes(A_2)
A_3=add_extremes(A_3)
A_4=add_extremes(A_4)

B_0=add_extremes(B_0)
B_1=add_extremes(B_1)
B_2=add_extremes(B_2)
B_3=add_extremes(B_3)
B_4=add_extremes(B_4)

C_1=add_extremes(C_1)
C_2=add_extremes(C_2)
C_3=add_extremes(C_3)
C_4=add_extremes(C_4)

D_0=add_extremes(D_0)
D_1=add_extremes(D_1)
D_2=add_extremes(D_2)
D_3=add_extremes(D_3)
D_4=add_extremes(D_4)

dC_1=derivative(Rcen,C_1)
dC_2=derivative(Rcen,C_2)
dC_3=derivative(Rcen,C_3)
dC_4=derivative(Rcen,C_4)

dD_0=derivative(Rcen,D_0)
dD_1=derivative(Rcen,D_1)
dD_2=derivative(Rcen,D_2)
dD_3=derivative(Rcen,D_3)
dD_4=derivative(Rcen,D_4)

Rcen_a=np.copy(Rcen)
Rcen_a[0]=30

fd_0=interp1d(Rcen, D_0/Rcen_a + dD_0)
fd_1=interp1d(Rcen, D_1/Rcen_a + dD_1)
fd_2=interp1d(Rcen, D_2/Rcen_a + dD_2)
fd_3=interp1d(Rcen, D_3/Rcen_a + dD_3)
fd_4=interp1d(Rcen, D_4/Rcen_a + dD_4)

fc_1=interp1d(Rcen, C_1/Rcen_a + dC_1)
fc_2=interp1d(Rcen, C_2/Rcen_a + dC_2)
fc_3=interp1d(Rcen, C_3/Rcen_a + dC_3)
fc_4=interp1d(Rcen, C_4/Rcen_a + dC_4)

fa_1=interp1d(Rcen,   A_1/Rcen_a)
fa_2=interp1d(Rcen, 2*A_2/Rcen_a)
fa_3=interp1d(Rcen, 3*A_3/Rcen_a)
fa_4=interp1d(Rcen, 4*A_4/Rcen_a)

fb_1=interp1d(Rcen,   B_1/Rcen_a)
fb_2=interp1d(Rcen, 2*B_2/Rcen_a)
fb_3=interp1d(Rcen, 3*B_3/Rcen_a)
fb_4=interp1d(Rcen, 4*B_4/Rcen_a)

Radius,Theta=polar_field(NN,L)

w1 =fd_0(Radius)
w1+=fd_1(Radius)*np.cos(Theta)
w1+=fd_2(Radius)*np.cos(2*Theta)
w1+=fd_3(Radius)*np.cos(3*Theta)
w1+=fd_4(Radius)*np.cos(4*Theta)

w1+=fc_1(Radius)*np.sin(Theta)
w1+=fc_2(Radius)*np.sin(2*Theta)
w1+=fc_3(Radius)*np.sin(3*Theta)
w1+=fc_4(Radius)*np.sin(4*Theta)

w2 =fa_1(Radius)*np.cos(Theta)
w2+=fa_2(Radius)*np.cos(2*Theta)
w2+=fa_3(Radius)*np.cos(3*Theta)
w2+=fa_4(Radius)*np.cos(4*Theta)

w2-=fb_1(Radius)*np.sin(Theta)
w2-=fb_2(Radius)*np.sin(2*Theta)
w2-=fb_3(Radius)*np.sin(3*Theta)
w2-=fb_4(Radius)*np.sin(4*Theta)


np.savetxt('Maps/'+name_file+'_exp.txt',1.0227*(w1+w2)*dA,fmt="%.4e",delimiter='\t')
