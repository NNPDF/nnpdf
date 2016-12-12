echo "Launching " $1 " jobs of configuration " $2

for (( c=1; c<=$1; c++ ))
do  
   echo "cd ${PWD}
./nnfit ${c} ${2}" > $c.run
   qsub -q long -l walltime=72:00:00 $c.run
done
