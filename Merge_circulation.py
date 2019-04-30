import numpy as np
import glob, os

name_file=str(np.loadtxt('Input.txt',usecols=(1,),dtype=bytes).astype(str))

vort_list=sorted(glob.glob('Files/'+name_file+'_vort_[0-9][0-9][0-9].txt'))
other_list=sorted(glob.glob('Files/'+name_file+'_other_[0-9][0-9][0-9].txt'))

if len(vort_list)>0:
    os.system('cat Files/'+name_file+'_vort_[0-9][0-9][0-9].txt > Files/'+name_file+'_vort.txt')
    for x in vort_list:
        os.system('rm '+x)
if len(other_list)>0:
    os.system('cat Files/'+name_file+'_other_[0-9][0-9][0-9].txt > Files/'+name_file+'_other.txt')
    for x in other_list:
        os.system('rm '+x)
