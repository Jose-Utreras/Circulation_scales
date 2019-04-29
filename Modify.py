import numpy as np
import os
Nbins,Nsteps=np.loadtxt('Input.txt',usecols=(5,12),dtype=int)

line="sscanf(line,\""
line2=""
for i in range(Nbins):
    line+=" %f"
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

lineas=[]
for k in range(8):
    line="\t\t\t\t\t\tsscanf(line,\"%lf %lf %lf"
    line2=""
    for i in range(Nbins):
        if (i%10==0)&(i>0):
            line2+="\n"
            for j in range(10):
                line2+="\t"
        line+=" %lf"
        line2+=", &tab_%d[%d]"%(k+1,i)
    line+="\", \n\t\t\t\t\t\t\t\t\t\t&n1%d,&n2%d,&pc%d \n \t\t\t\t\t\t\t\t\t\t"%(k+1,k+1,k+1)+line2+");}\n"
    lineas.append(line)


f_in=open('mcmc_chains.c','r')
f_out=open('out.c','w')
switch=False
switch2=False
k=0
for x in f_in:

    if switch:
        f_out.write(lineas[k])
        k+=1
    else:
        if switch or switch2:
            pass
        else:
            f_out.write(x)

    if x[2:12]=='READSCALES':
        switch=True
        switch2=True
    else:
        switch=False
    if x[2:8]=='SCALES':
        switch2=False
        f_out.write('//SCALES \n')


f_in.close()
f_out.close()

os.system('mv out.c mcmc_chains.c')



line="sscanf(line,\""
line2=""
line3=""
line4=""
for i in range(Nbins):
    line+=" %f %f %f"
    line2+=", &delta[%d]" %i
    line3+=", &sigma[%d]" %i
    line4+=", &error[%d]" %i
line+="\""+line2+line3+line4+");\n"
#print(line)

f_in=open('mcmc_chains_copy.c','r')
f_out=open('out.c','w')
switch=False
for x in f_in:
    if switch:
        f_out.write(line)
    else:
        f_out.write(x)
    if x[2:12]=='READARRAYS':
        switch=True
    else:
        switch=False

f_in.close()
f_out.close()

os.system('mv out.c mcmc_chains.c')

################################################################
###### Modifiying steps in parameter_space_generator_mpi.c #####
################################################################

line="#define steps %d\n" %Nsteps

f_in=open('parameter_space_generator_mpi.c','r')
f_out=open('out.c','w')

for x in f_in:
    if x[:13]=="#define steps":
        f_out.write(line)
    else:
        f_out.write(x)

f_in.close()
f_out.close()

os.system('mv out.c parameter_space_generator_mpi.c')

f_in=open('mcmc_chains.c','r')
f_out=open('out.c','w')

for x in f_in:
    if x[:13]=="#define steps":
        f_out.write(line)
    else:
        f_out.write(x)

f_in.close()
f_out.close()

os.system('mv out.c mcmc_chains.c')
