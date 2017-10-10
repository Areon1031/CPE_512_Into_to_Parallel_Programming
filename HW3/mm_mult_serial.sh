#!/bin/bash
# Note: for batch run_script queue on jetson
# after correctness has been verified remove
# matrix data output in program but output data size 
#and measured run time from the rank 0 process
# then compile on jetson using
# mpic++ mm_mult_serial.cpp -o mm_mult_serial -O3 -lm
# make sure this script has execute rights
#   chmod 744 mm_mult_serial.sh
# Then run this batch script using
#   run_script mm_mult_serial.sh
# where you specify 1 as the number of nodes
#output will be captured in mm_mult_serial.txt file
./mm_mult_serial 256  >mm_mult_serial.txt
./mm_mult_serial 512  >>mm_mult_serial.txt
./mm_mult_serial 768  >>mm_mult_serial.txt 
./mm_mult_serial 1024 >>mm_mult_serial.txt 
./mm_mult_serial 1280 >>mm_mult_serial.txt 
./mm_mult_serial 1536 >>mm_mult_serial.txt 
./mm_mult_serial 1792 >>mm_mult_serial.txt 
./mm_mult_serial 2048 >>mm_mult_serial.txt 
./mm_mult_serial 2304 >>mm_mult_serial.txt 
./mm_mult_serial 2560 >>mm_mult_serial.txt 
./mm_mult_serial 2816 >>mm_mult_serial.txt 
./mm_mult_serial 3072 >>mm_mult_serial.txt 
./mm_mult_serial 3328 >>mm_mult_serial.txt 
./mm_mult_serial 3584 >>mm_mult_serial.txt 
