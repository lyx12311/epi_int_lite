#!/bin/csh
# SLURM docs: http://www.glue.umd.edu/hpcc/help/jobs.html#queues

#SBATCH --ntasks=1
#SBATCH --switches=1 # if possible!
#SBATCH -t 20:00:00 # wallclock limit SS or MM:SS or HH:MM:SS, DD-HH, etc.

hostname
date
pwd

echo JOBID "$SLURM_JOBID"
echo JOB_NAME "$SLURM_JOB_NAME"
echo SUBMIT_DIR "$SLURM_SUBMIT_DIR"
echo JOB_NODELIST "$SLURM_JOB_NODELIST"
echo SUBMIT_HOST "$SLURM_SUBMIT_HOST"
echo JOB_NUM_NODES "$SLURM_JOB_NUM_NODES"
echo CPUS_ON_NODE "$SLURM_CPUS_ON_NODE"
echo NTASKS "$SLURM_NTASKS"

./epi_int_lite.py > outputFile
