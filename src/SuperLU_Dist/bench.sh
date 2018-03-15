#!/bin/bash

exe="pddrive_ABglobal"
declare -a ORDERINGS=("COLAMD" "MMD_ATA" "MMD_AT_PLUS_A")
#declare -a ORDERINGS=("NATURAL" "COLAMD" "MMD_ATA" "MMD_AT_PLUS_A")


prefix=../linearSystems/hqp3_60/FullUMF_AMD/FullUMF_3_
A=${prefix}A.mtx
b=${prefix}b.mtx
x=${prefix}x.mtx
R=3
max_threads=2
max_mpi=2

echo "         date: ($date)"
echo "            A: ${A}"
echo "            b: ${b}"
echo "            x: ${x}"
echo "            R: ${R}"
echo "      max Th.: ${max_threads}"
echo "max MPI Procs: ${max_mpi}"

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
   for np in `seq 1 ${max_mpi}`; do
      echo "j : $np"
      if [ "$use_srun" = true ]; then
        srun -c ${np} ./${exe} -A ${A} -b ${b} -x ${x} -R ${R}
      else
        mpirun -np ${np} ./${exe} -A ${A} -b ${b} -x ${x} -R ${R}
        # -c -r
      fi
   done
done
export ORDERING=""
