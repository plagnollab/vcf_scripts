#!/bin/bash
#$ -S /bin/bash
#$ -o /dev/null
#$ -e /dev/null
#$ -cwd
#$ -V
#$ -R y
#$ -l tmem=2G,h_vmem=1G
#$ -l h_rt=1:0:0
#$ -t 1-25
#$ -tc 25
set -u
set -x
let "SGE_TASK_ID=$SGE_TASK_ID-1"
scriptname=split
mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err
args=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
chr=${args[$SGE_TASK_ID]}
exec >${scriptname}.qsub.out/${scriptname}_${chr}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${chr}_${JOB_ID}.err
rm -f chr$chr.vcf.gz*
file=$1
cat vcf_header.txt <(tabix $file $chr) | bgzip -c -f > chr$chr.vcf.gz
tabix -f -p vcf chr$chr.vcf.gz


