#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --export=none  # Ensures job gets a fresh login environment
#SBATCH --mail-user=robohyun66@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=subsample
#SBATCH --output=./logs2/subsample.%A_%a.log
#SBATCH --error=./logs2/subsample.%A_%a.err

# Set the session up to use R
module load gcc/8.3.0
module load openblas/0.3.8
module load r/3.6.3

echo This is SLURM task $SLURM_ARRAY_TASK_ID

Rscript subsample.R\
	arraynum=$SLURM_ARRAY_TASK_ID\
	arraynum_max=$SLURM_ARRAY_TASK_MAX\
  	isim=$isim\
  	nsim=$nsim\
	subsample_size=$subsample_size\
	summ=$summ\

exit

## Use this to run the main experiments
# The arguments, in order, are (nsim, subsample_size, arraynum_max)
run (){
    nsim=$1
    subsample_size=$2
    arraynum_max=$3
    # for isim in 96; do
    for isim in $(seq 1 1 $nsim); do
	# file_count=$(ls -f ~/scratchdir/output/subsample-random-time/subsample-b-"$subsample_size"/sim-"$isim" | wc -l)
	# if(($file_count == "2945")); then
	# 	continue
	# fi
  	printf "$isim "
  	jid1=$(sbatch --export=isim=$isim,nsim=$nsim,subsample_size=$subsample_size,summ=0 --array=1 run-subsample.slurm)
  	jid2=$(sbatch --export=isim=$isim,nsim=$nsim,subsample_size=$subsample_size,summ=0 --array=1-$arraynum_max --dependency=afterany:${jid1##* } run-subsample.slurm)
  	# jid2=$(sbatch --export=isim=$isim,nsim=$nsim,subsample_size=$subsample_size,summ=0 --array=1-$arraynum_max run-subsample.slurm)
  	jid3=$(sbatch --export=isim=$isim,nsim=$nsim,subsample_size=$subsample_size,summ=1 --dependency=afterany:${jid2##* } --array=1 run-subsample.slurm)
    done
    return 0
}

run 200 30 20

run 200 60 20

run 200 100 20

run 200 150 20

run 200 200 10

run 200 250 10


for sz in 30 60 100 200 250; do
    run 100 $sz 10
done


# Example settings
isim=1
nsim=1
arraynum_max=1
arraynum=1
summ=TRUE
subsample_size=150
