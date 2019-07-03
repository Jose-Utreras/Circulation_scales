import numpy as np
import glob, os
from astropy.table import Table , Column 

name_file=str(np.loadtxt('Input.txt',usecols=(1,),dtype=bytes).astype(str))

sim_list=sorted(glob.glob('Files/'+name_file+'*.txt'))
sim_list+=sorted(glob.glob('Samplers/'+name_file+'*.txt'))
sim_list+=sorted(glob.glob('Maps/'+name_file+'*.txt'))

"""
for sim_file in sim_list:
    print('rm '+sim_file)
    os.system('rm '+sim_file)
"""
os.system('mkdir Chains/'+name_file)
os.system('mv Chains/'+name_file+'*.txt Chains/'+name_file+'/')

lista1=sorted(glob.glob('Chains/'+name_file+'/*.txt'))

for li in lista1:
    tab=np.loadtxt(li,dtype={'names': ('n1', 'n2', 'pc','vo','alpha'),
                        'formats': ('f4','f4','f4', 'f4', 'f4')})
    new_table=Table()
    new_table['n1']=Column(tab['n1'])
    new_table['n2']=Column(tab['n2'])
    new_table['pc']=Column(tab['pc'])
    new_table['vo']=Column(tab['vo'])
    new_table.write(li.split('.txt')[0],path='data',format='hdf5',overwrite=True)

#os.system('rclone copy Chains/'+name_file+' uchile:Chains/'+name_file)
