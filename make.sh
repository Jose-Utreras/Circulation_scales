mkdir Maps
mkdir Files
mkdir Chains
mkdir Samplers
mpicc -lm Circulation.c -o circ
mpicc -lm Circulation_model_and_exp.c -o circ_me
gcc -lm histogram.c -o histogram
mpicc -lm mcmc_chains.c -o mcmc
mpicc -lm parameter_space_generator_mpi.c -o psgm
python Modify 
