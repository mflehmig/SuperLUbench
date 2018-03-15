#!/bin/bash

exe="dlinsolx"
declare -a ORDERINGS=("NATURAL" "COLAMD" "MMD_ATA" "MMD_AT_PLUS_A")


prefix=../linearSystems/hqp3_60/FullUMF_AMD/FullUMF_3_
A=${prefix}A.mtx
b=${prefix}b.mtx
x=${prefix}x.mtx
R=5

echo "   A: ${A}"
echo "   b: ${b}"
echo "   x: ${x}"
echo "   R: ${R}"

# Executable is compiled?
if [ ! -f "$exe" ]; then
  echo "Error: Executable $exe not found. Call make first."
  exit 1
fi


# Is srun command available? If so, use it!
if $(command -v srun >/dev/null 2>&1 ); then
  echo "srun: yes"
  exec_cmd="srun -c 1"
else
  echo "srun: no"
  exec_cmd=""
fi


for i in "${ORDERINGS[@]}"; do
   echo -e "\n\nORDERING: $i"
   export ORDERING=${i}
   ${exec_cmd} ./${exe} -A ${A} -b ${b} -x ${x} -R ${R}
done
export ORDERING=""
