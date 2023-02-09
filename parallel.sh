#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Illegal number of parameters, please give num of nodes, num of cpus/node, dataset, k"
    exit 2
fi

echo "Number of nodes: $1"
echo "Number of cpus/node: $2"
echo "Dataset path: $3"
echo "K: $4"

file="script.sh"
executable="executable"


echo "#!/bin/bash" > $file
echo "" >> $file
echo "#PBS -l select=$1:ncpus=$2:mem=2gb" >> $file
echo "#PBS -l walltime=0:03:00" >> $file
echo "#PBS -q short_cpuQ" >> $file
echo "#PBS -l place=scatter:excl" >> $file
echo "" >> $file
echo "module load mpich-3.2" >> $file
echo "mpirun.actual -n $1 ./$executable $3 $2 $4" >> $file
chmod +x ./$file
echo ""
echo "Queuing execution, ID:"
qsub $file
