#!/bin/bash
exe="pddrive_ABglobal"
if [ -f "$exe" ]; then
  OMP_NUM_THREADS=1 mpirun -np 2 ./$exe -r 1 -c 2 -A ../../example/A.mtx -b ../../example/b.mtx  -x ../../example/x.mtx -v
else
  echo "Error: Executable $exe not found. Call make first."
fi
