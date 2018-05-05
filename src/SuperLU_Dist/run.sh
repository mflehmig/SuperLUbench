#!/bin/bash

exe="pddrive_ABglobal"
declare -a ORDERINGS=("MMD_AT_PLUS_A")
#declare -a ORDERINGS=("NATURAL" "COLAMD" "MMD_ATA" "MMD_AT_PLUS_A")

if [ "$#" -ne 4 ]; then
  echo "Usage: sh bench.sh FILE_PREFIX NUM_ITER MAX_THREADS MAX_MPI"
else

prefix=$1
A=${prefix}A.mtx
b=${prefix}b.mtx
x=${prefix}x.mtx
R=$2
threads=$3
mpi=$4

echo "         date: ($date)"
echo "            A: ${A}"
echo "            b: ${b}"
echo "            x: ${x}"
echo "            R: ${R}"
echo "      Threads: ${threads}"
echo "    MPI Procs: ${mpi}"

# Executable is compiled?
if [ ! -f "$exe" ]; then
  echo "Error: Executable $exe not found. Call make first."
  exit 1
fi


# Is srun command available? If so, use it!
if $(command -v srun >/dev/null 2>&1 ); then
  echo "srun: yes"
  use_srun=true
else
  echo "srun: no"
  use_srun=false
fi


for i in "${ORDERINGS[@]}"; do
   echo -e "\n\nORDERING: $i"
   export ORDERING=${i}
   if [ "$use_srun" = true ]; then
      srun --ntasks ${mpi} ./${exe} -A ${A} -b ${b} -x ${x} -R ${R}
   else
      mpirun -np ${np} ./${exe} -A ${A} -b ${b} -x ${x} -R ${R}
      # -c -r
   fi
done
export ORDERING=""

fi
