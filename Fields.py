import numpy as np
import math
from common_functions import *
import sys
import yt
from yt.units import kpc,pc,km,second,yr,Myr,Msun
from yt.fields.api import ValidateParameter
from yt.utilities.physical_constants import G
from yt.fields.derived_field import \
    ValidateGridType, \
    ValidateParameter, \
    ValidateSpatial, \
    NeedsParameter
from yt.funcs import \
    just_one

from scipy.interpolate import interp1d

sl_left = slice(None, -2, None)
sl_right = slice(2, None, None)
div_fac = 2.0

sl_center = slice(1, -1, None)
ftype='gas'

vort_validators = [ValidateSpatial(1,
                        [(ftype, "velocity_x"),
                         (ftype, "velocity_y"),
                         (ftype, "velocity_z")])]


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

def _Disk_H(field, data):
    center = data.get_field_parameter('center')
    z = data["z"] - center[2]
    return np.abs(z)
yt.add_field("Disk_H",
             function=_Disk_H,
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


def _vorticity_z(field, data):
    f  = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
          data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
          / (div_fac*just_one(data["index", "dx"]))
    f -= (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
          data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
          / (div_fac*just_one(data["index", "dy"]))
    new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                          dtype=np.float64),
                            f.units)
    new_field[sl_center, sl_center, sl_center] = f
    return new_field
yt.add_field((ftype, "vorticity_z"),function=_vorticity_z,
                       units="1/yr",
                       validators=vort_validators,take_log=True)

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

Disk = ds.disk('c', [0., 0., 1.],(L/1.0e3, 'kpc'), (1, 'kpc'))

VC=Disk['vc'].in_units("pc/Myr")
R=Disk['Disk_Radius'].in_units("pc")
Masa=Disk['cell_mass'].in_units("Msun")

R_map=radial_map_N(NN,NN)*L/NN


x,y=Vorticity_profile(VC,R,500,0,0.75*L,np.mean)
fun=interp1d(x,y)
vorticity_map=fun(R_map)
np.savetxt('Maps/'+name_file+'_av_vort.txt',vorticity_map*dA,fmt="%.4e",delimiter='\t')


x,y=Sigma_profile(Masa,R,500,0,0.75*L)
fun=interp1d(x,y)
sigma_map=fun(R_map)
np.savetxt('Maps/'+name_file+'_av_sigma.txt',sigma_map,fmt="%.4e",delimiter='\t')

del Disk,sigma_map,vorticity_map

dd=ds.all_data()
width = (float(L), 'pc')
print('### Projection ###')
disk_dd = dd.cut_region(["obj['Disk_H'].in_units('pc') < 1.0e3"])
proj = ds.proj('vorticity_z', 2,data_source=disk_dd,weight_field='density')

res = [NN, NN]
frb = proj.to_frb(width, res, center=[0.5,0.5,0.5])

print('### Vorticity Map ###')
imvort=frb['vorticity_z'].in_units('1/Myr')*dA

np.savetxt('Maps/'+name_file+'_vort.txt',imvort,fmt="%.4e",delimiter='\t')

proj = ds.proj('density', 2,data_source=disk_dd,method='integrate')
frb = proj.to_frb(width, res, center=[0.5,0.5,0.5])

imsigma=frb['density'].in_units('Msun/pc**2')

np.savetxt('Maps/'+name_file+'_sigma.txt',imsigma,fmt="%.4e",delimiter='\t')
