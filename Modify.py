import numpy as np
import os
Nbins=np.loadtxt('Input.txt',usecols=(5),dtype=int)

line="sscanf(line,\""
line2=""
for i in range(Nbins):
    line+=" %lf"
    line2+=", &delta[%d]" %i
line+="\""+line2+");\n"
#print(line)

f_in=open('parameter_space_generator_mpi.c','r')
f_out=open('out.c','w')
switch=False
for x in f_in:
    if switch:
        f_out.write(line)
    else:
        f_out.write(x)
    if x[2:12]=='READSCALES':
        switch=True
    else:
        switch=False

f_in.close()
f_out.close()

os.system('mv out.c parameter_space_generator_mpi.c')
