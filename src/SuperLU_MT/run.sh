#!/bin/bash

exe="pdlinsolx"
declare -a ORDERINGS=("COLAMD")
#declare -a ORDERINGS=("NATURAL" "COLAMD" "MMD_ATA" "MMD_AT_PLUS_A")

if [ "$#" -ne 3 ]; then
  echo "Usage: sh bench.sh FILE_PREFIX NUM_ITER THREADS"
else

prefix=$1
A=${prefix}A.mtx
b=${prefix}b.mtx
x=${prefix}x.mtx
R=$2
threads=$3

echo "      A: ${A}"
echo "      b: ${b}"
echo "      x: ${x}"
echo "      R: ${R}"
echo "Threads: ${threads}"

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
   echo "j : $th"
   if [ "$use_srun" = true ]; then
     srun -c ${threads} ./${exe} -A ${A} -b ${b} -x ${x} -R ${R} -p ${threads} | grep "==="
   else
     ./${exe} -A ${A} -b ${b} -x ${x} -R ${R} -p ${th} | grep "==="
   fi
done
export ORDERING=""

fi
