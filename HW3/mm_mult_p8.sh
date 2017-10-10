#!/bin/bash
# Note: for batch run_script queue on jetson
# after correctness has been verified remove
# matrix data output in program but output data size 
#and measured run time from the rank 0 process
# then compile on jetson using
# mpic++ mm_mult_mpi.cpp -o mm_mult_mpi -O3 -lm
# make sure this script has execute rights
#   chmod 744 mm_mult_p8.sh
# Then run this batch script using
#   run_script mm_mult_p8.sh
# where you specify 8 as the number of nodes
#output should be captured in mm_p8.txt file
srun ./mm_mult_mpi 256  >mm_p8.txt
srun ./mm_mult_mpi 512  >>mm_p8.txt
srun ./mm_mult_mpi 768  >>mm_p8.txt 
srun ./mm_mult_mpi 1024 >>mm_p8.txt 
srun ./mm_mult_mpi 1280 >>mm_p8.txt 
srun ./mm_mult_mpi 1536 >>mm_p8.txt 
srun ./mm_mult_mpi 1792 >>mm_p8.txt 
srun ./mm_mult_mpi 2048 >>mm_p8.txt 
srun ./mm_mult_mpi 2304 >>mm_p8.txt 
srun ./mm_mult_mpi 2560 >>mm_p8.txt 
srun ./mm_mult_mpi 2816 >>mm_p8.txt 
srun ./mm_mult_mpi 3072 >>mm_p8.txt 
srun ./mm_mult_mpi 3328 >>mm_p8.txt 
srun ./mm_mult_mpi 3584 >>mm_p8.txt 
