CC=gcc
MPICC=mpicc
CFLAGS=-c -Wall -std=c99 -O0 -g

# SuperLU
#I_SUPERLU=-I/usr/include/SuperLU/
#LD_SUPERLU=-lsuperlu
I_SUPERLU=-I/home/mf/Downloads/SuperLU_5.2.1/SRC/
LD_SUPERLU=-L/home/mf/Downloads/SuperLU_5.2.1/build_debug/SRC/ -lsuperlu -Wl,-rpath,/home/mf/Downloads/SuperLU_5.2.1/build_debug/SRC/

# SuperLU_MT
I_SUPERLUMT=-I/usr/include/SuperLUMT
LD_SUPERLUMT=-lsuperlumt_d
#I_SUPERLUMT=-I/home/mf/opt/PARADOM/AP3/SuperLU_MT_3.1/SRC/
#LD_SUPERLUMT=-l/home/mf/opt/PARADOM/AP3/SuperLU_MT_3.1/lib/libsuperlu_mt_OPENMP.a

# SuperLU_DIST
#I_SUPERLUDIST=-I/usr/include/openmpi-x86_64/superlu_dist/
#LD_SUPERLUDIST=-L/usr/lib64/openmpi/lib/ -lsuperlu_dist

I_SUPERLUDIST=-I/home/mf/opt/PARADOM/AP3/SuperLU_DIST_5.1.3/BUILD/include
LD_SUPERLUDIST=-L/home/mf/opt/PARADOM/AP3/SuperLU_DIST_5.1.3/BUILD/lib -lsuperlu_dist -Wl,-rpath,/home/mf/opt/PARADOM/AP3/SuperLU_DIST_5.1.3/BUILD/lib
