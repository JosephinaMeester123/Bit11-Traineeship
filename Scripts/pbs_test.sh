#!/usr/bin/env bash
#PBS -d /home/jmeester/WGS_Jeannette
#PBS -l nodes=1:ppn=4
#PBS -l mem=12gb
#PBS -l walltime=48:00:00
#PBS -N myscript
#PBS -o /home/jmeester/WGS_Jeannette/myscript.stdout.$PBS_JOBID
#PBS -e /home/jmeester/WGS_Jeannette/myscript.stderror.$PBS_JOBID
#PBS -m abe
#PBS -q batch
#PBS -M josephina.meester@uantwerpen.be
#PBS -V
#PBS -A default

echo 'Running on : ' `hostname`
echo 'Start Time : ' `date`
echo 'Command:'
echo '========'
echo ‘Test job, please ignore’

sleep 10

echo 'End Time : ' `date`
printf 'Execution Time = %dh:%dm:%ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60))


