import numpy as np
import os, sys

cmdargs      = sys.argv
new_name     =   cmdargs[-1]

input = np.loadtxt('Input.txt',dtype=bytes).astype(str)

input[1]="%s" %new_name
np.savetxt('Input.txt',np.atleast_2d(input),fmt='%s',delimiter='\t')
