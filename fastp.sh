#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00

FASTQC=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software
export PATH=$PATH:$FASTQC

NPROCS=$SLURM_NPROCS


mkdir -p MULTIQC_VIZ_DIR

mkdir -p FASTP_OUT_DIR

mkdir -p CHKPNT_DIR

for f in $(ls *_L001_R1_001.fastq.gz)
do
basename=${f##*/}
bs="${basename%_L001_R1_001.fastq.gz}"
infile="${f%_L001_R1_001.fastq.gz}"

left_file=${infile}_L001_R1_001.fastq.gz
right_file=${infile}_L001_R2_001.fastq.gz

if [ ! -f "CHKPNT_DIR/${bs}.chkpt" ]; then

    call="fastp --thread $NPROCS --detect_adapter_for_pe \
    --json MULTIQC_VIZ_DIR/${bs}_fastp.json \
    --html MULTIQC_VIZ_DIR/${bs}_fastp.html \
    -i $left_file -I $right_file \
    -o FASTP_OUT_DIR/${bs}_R1.fq.gz  -O FASTP_OUT_DIR/${bs}_R2.fq.gz"

    echo "Running fastp in '${bs}' sample."

    echo $call

    eval $call

    unlink $left_file
    unlink $right_file

    touch  CHKPNT_DIR/${bs}.chkpt

else
    echo "'CHKPNT_DIR/${bs}.chkpt' already exists."
    echo "Continue with next sample."
fi


done

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

$WDM/multiqc --module fastp MULTIQC_VIZ_DIR  -o MULTIQC_VIZ_DIR --config multiqc_info.conf


exit
