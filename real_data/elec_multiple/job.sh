#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=24:00:00
#SBATCH --array=0-147
#SBATCH --output=out/slurm_%j.out
#SBATCH --error=err/slurm_%j.err
#SBATCH --partition=epyc-64

sh parameter.sh > parameter.txt
module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.2.3

for i in $(seq 1 $SLURM_NTASKS | head -n  $SLURM_NTASKS)
do
    line_to_read=$((SLURM_ARRAY_TASK_ID*$SLURM_NTASKS+$i))
    var=$(awk "NR==${line_to_read}" parameter.txt)
    read -r a <<<$(echo $var)
    #Rscript -e 'install.packages(c("usedist"), repos="https://cloud.r-project.org")'
    srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK Rscript elec_multiple.R $a &
done
wait


