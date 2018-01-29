echo "Launching " $1 " jobs of configuration " $2

for (( c=1; c<=$1; c++ ))
do  
   echo "cd ${PWD}
./nnfit ${c} ${2}" > ${c}_${2}.run
   qsub -q generic -l walltime=24:00:00 ${c}_${2}.run
done
