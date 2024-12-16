#!/bin/bash

#SBATCH --job-name=ArrayJob
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-10
#SBATCH --qos=short


# load module
module load GCCcore/11.2.0 # bwa, GATK & R module dependency
module load BWA/0.7.17
module load GATK/4.2.5.0-Java-11
module load picard/2.25.1-Java-11
module load GCC/11.2.0 # samtools dependency
module load SAMtools/1.14
module load OpenMPI/4.1.1 # R module dependency
module load R/4.1.2


# reads list: c1 head "arrayTaskID", body 1 2 3...   c2 R1 c3 R2, out (name of prefix. (inclune some header in all columns)
rl=veluto.csv
# AJ4 + BMV58.gff
refGen= comb.fasta 



# extract name of r1
r1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rl)
echo this is r1 $r1


# extract name of r2
r2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $rl)
echo this is r2 $r2

# name of the outputs
pref=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4'} $rl)


# extract name of output
out=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
echo this is out $out

python variant_calling_1.py -p $pref -numth $numth -refgen $refgen -r1 $r1 -r2 $r2

