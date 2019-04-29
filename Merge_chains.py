import numpy as np
import glob, os
import sys
cmdargs = sys.argv

Rmax        =   float(cmdargs[-1])
Rmin        =   float(cmdargs[-2])

name_file=str(np.loadtxt('Input.txt',usecols=(1),dtype=bytes).astype(str))

name_file+='_%05d_%05d' %(int(Rmin),int(Rmax))
chain_list=sorted(glob.glob('Chains/'+name_file+'_chain[0-9][0-9][0-9].txt'))



os.system('cat Chains/'+name_file+'_chain*.txt > Chains/'+name_file+'_chain.txt')

for chain in chain_list:
    os.system('rm '+chain)
