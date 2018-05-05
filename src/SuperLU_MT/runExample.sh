#!/bin/bash
exe="pdlinsolx"
if [ -f "$exe" ]; then
  OMP_NUM_THREADS=2 ./$exe -p 2 -A ../../example/A.mtx -b ../../example/b.mtx  -x ../../example/x.mtx -v
else
  echo "Error: Executable $exe not found. Call make first."
fi
