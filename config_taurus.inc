CC=icc
MPICC=mpicc
CFLAGS=-c -Wall -std=c99 -O2 -march=native

# SuperLU
I_SUPERLU=-I${SUPERLU_INC}
LD_SUPERLU=-L${SUPERLU_LIB} -lsuperlu -Wl,-rpath,${SUPERLU_LIB}

# SuperLU_MT
I_SUPERLUMT=-I${SUPERLU_MT_INC}
LD_SUPERLUMT=${SUPERLU_MT_LIB}/superlu_OPENMP.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#I_SUPERLUMT=-I/home/mf/opt/PARADOM/AP3/SuperLU_MT_3.1/SRC/
#LD_SUPERLUMT=-l/home/mf/opt/PARADOM/AP3/SuperLU_MT_3.1/lib/libsuperlu_mt_OPENMP.a

# SuperLU_DIST
I_SUPERLUDIST=-I${SUPERLU_DIST_INC}
LD_SUPERLUDIST=-L${SUPERLU_DIST_LIB} -lsuperlu_dist

#I_SUPERLUDIST=-I/home/mf/opt/PARADOM/AP3/SuperLU_DIST_5.1.3/BUILD/include
#LD_SUPERLUDIST=-L/home/mf/opt/PARADOM/AP3/SuperLU_DIST_5.1.3/BUILD/lib -lsuperlu_dist -Wl,-rpath,/home/mf/opt/PARADOM/AP3/SuperLU_DIST_5.1.3/BUILD/lib
