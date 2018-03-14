#!/bin/bash
exe="dlinsolx"
if [ -f "$exe" ]; then
  ./$exe -A ../../example/A.mtx -b ../../example/b.mtx  -x ../../example/x.mtx -v
else
  echo "Error: Executable $exe not found. Call make first."
fi
