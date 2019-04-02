import numpy as np
import math
from common_functions import *
import sys
import yt
from yt.units import kpc,pc,km,second,yr,Myr,Msun
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

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


def V_curve(R,Vo,R1):
    return Vo*np.arctan(R/R1)

def Velocity_curve(R,Vo,R1,R2):
    return Vo*np.arctan(R/R1)*np.exp(-R/R2)

def Vorticity_curve(R,Vo,R1,R2):

    result=Vo*np.exp(-R/R2)/(R+0.1)
    right=np.arctan(R/R1)*(1.0-R/R2)
    left=R*R1/(R1**2+R**2)

    result*= left + right

    return result

cmdargs = sys.argv
L=float(cmdargs[-1])                        # Image size in parsecs
name_file=cmdargs[-2]

dir=str(np.loadtxt('Input.txt',usecols=(0),dtype=bytes).astype(str))

directory=dir+'/'+name_file+'/G-'+name_file[-4:]

ds = yt.load(directory)

grids=ds.refine_by**ds.index.max_level*ds.domain_dimensions[0]

DX=ds.arr(1, 'code_length')
DX.convert_to_units('pc')
dr=DX/grids

NN=int(L/dr)
dA=((L/NN)**2)
print(L/NN)
Disk = ds.disk('c', [0., 0., 1.],(40, 'kpc'), (1, 'kpc'))

radius=Disk['Disk_Radius'].in_units("pc")
vel_cir=Disk['vc'].in_units('km/s')
Redges=np.arange(0*pc,35000*pc,500*pc)
Rcen=0.5*(Redges[1:]+Redges[0:-1])
N=len(Rcen)

V_16 = np.zeros(N)
V_50 = np.zeros(N)
V_84 = np.zeros(N)

for k in range(N):
    ring=(Redges[k]<radius)&(Redges[k+1]>radius)

    if len(ring[ring])<4:
        V_16[k]  = np.nan
        V_50[k]  = np.nan
        V_84[k]  = np.nan

    else:
        vc_ring=vel_cir[ring]

        V_16[k], V_50[k], V_84[k] = np.percentile(vc_ring,[16,50,84])

Rcen=add_extremes(Rcen)
V_16=add_extremes(V_16)
V_50=add_extremes(V_50)
V_84=add_extremes(V_84)

radius  =   Rcen[Rcen<20000]
V_16    =   V_16[Rcen<20000]
V_50    =   V_50[Rcen<20000]
V_84    =   V_84[Rcen<20000]

V_err=0.5*(V_84-V_16)
V_err[0]=V_err[1:].min()

p_0, cov_0 = curve_fit(V_curve,radius,V_50,sigma=V_err,p0=[V_50.max(),np.mean(radius)])
p_0, cov_0 = curve_fit(Velocity_curve,radius,V_50,sigma=V_err,p0=[p_0[0],p_0[1],p_0[1]])

Radius,Theta=polar_field(NN,L)

wort=Vorticity_curve(Radius,*p_0)

np.savetxt('Maps/'+name_file+'_model.txt',1.0227*wort*dA,fmt="%.4e",delimiter='\t')
