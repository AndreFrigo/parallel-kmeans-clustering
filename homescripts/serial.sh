#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, please give dataset, k"
    exit 2
fi

echo "Dataset path: $1"
echo "K: $1"

file="script.sh"
executable="executable"

echo "#!/bin/bash" > $file
echo "" >> $file
echo "#PBS -l select=1:ncpus=1:mem=2gb" >> $file
echo "#PBS -l walltime=0:05:00" >> $file
echo "#PBS -q short_cpuQ" >> $file
echo "" >> $file
echo "./$executable $1 $2" >> $file
chmod +x ./$file
echo ""
echo "Queuing execution, ID:"
qsub $file
