CC=gcc
MPICC=mpicxx
CFLAGS=-c -Wall -std=c11 -O0 -g
LDFLAGS=-lm -fopenmp

# SuperLU
#I_SUPERLU=-I/usr/include/SuperLU/
#LD_SUPERLU=-lsuperlu
#I_SUPERLU=-I/home/mf/Downloads/SuperLU_5.2.1/SRC/
#LD_SUPERLU=-L/home/mf/Downloads/SuperLU_5.2.1/build_debug/SRC/ -lsuperlu -Wl,-rpath,/home/mf/Downloads/SuperLU_5.2.1/build_debug/SRC/

# SuperLU_MT
#I_SUPERLUMT=-I/usr/include/SuperLUMT
#LD_SUPERLUMT=-lsuperlumt_d
#I_SUPERLUMT=-I/home/mf/opt/PARADOM/AP3/SuperLU_MT_3.1/SRC/
#LD_SUPERLUMT=-l/home/mf/opt/PARADOM/AP3/SuperLU_MT_3.1/lib/libsuperlu_mt_OPENMP.a

# SuperLU_DIST
#I_SUPERLUDIST=-I/usr/include/openmpi-x86_64/superlu_dist/
#LD_SUPERLUDIST=-L/usr/lib64/openmpi/lib/ -lsuperlu_dist

I_SUPERLUDIST=-I/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/include
LD_SUPERLUDIST=-L/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/lib  -lsuperlu_dist -Wl,-rpath,/home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/lib /home/thomas/Code/libraries/builds/SuperLU_DIST_6.1.0/lib/libblas.a -lmetis
