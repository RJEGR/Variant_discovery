#!/bin/bash
#SBATCH --job-name=qtrim
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

NPROCS=$SLURM_NPROCS

TRIMM=/LUSTRE/apps/bioinformatica/Trimmomatic-0.40/
trimmomatic="java -jar $TRIMM/trimmomatic-0.40-rc1.jar"

for f in $(ls *_L001_R1_001.fastq.gz)
do
basename=${f##*/}
bs="${basename%_L001_R1_001.fastq.gz}"
infile="${f%_L001_R1_001.fastq.gz}"

left_file=${infile}_L001_R1_001.fastq.gz
right_file=${infile}_L001_R2_001.fastq.gz

call="$trimmomatic PE -threads $NPROCS -phred33 \
    -summary ${bs}-summary.txt \
    $left_file $right_file \
    ${bs}_R1.P.qtrim.fq.gz ${bs}_R1.UP.qtrim.fq.gz \
    ${bs}_R2.P.qtrim.fq.gz ${bs}_R2.UP.qtrim.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36 LEADING:5 TRAILING:5;"

# echo $call

eval $call

unlink $left_file
unlink $right_file

done

FASTQC=/LUSTRE/apps/bioinformatica/FastQC_v0.12.1/
export PATH=$PATH:$FASTQC

mkdir -p fastqc

fastqc *P.qtrim.fq.gz -t $NPROCS --nogroup -o fastqc

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

mkdir -p MULTIQC_VIZ_DIR

$WDM/multiqc fastqc -o MULTIQC_VIZ_DIR --config multiqc_info.conf

exit
