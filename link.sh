#!/bin/bash


COMPILE="-o SuperLU_Dist_bench \
    /home/thomas/Code/SuperLUbench/src/SuperLU_Dist/pddrive_ABglobal.o \
    -std=c11 -O0 -g -fopenmp -lm -Wall -lpthread -latlas -lopenblas -lsuperlu_dist -lmetis \
    -L/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/lib \
    CMakeFiles/SuperLU_Dist_bench.dir/src/Util.c.o  \
    CMakeFiles/SuperLU_Dist_bench.dir/src/dreadMM.c.o"

echo $COMPILE

mpicxx $COMPILE

# mpicxx -c -Wall -std=c11 -O0 -g -I/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/include ../Util.c -o ../Util.o

# es fehlt
#  sluddefs.h
#  
#    -Wl,-rpath,/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/lib \
