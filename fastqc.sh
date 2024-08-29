#!/bin/sh
## Directivas
#SBATCH --job-name=fastqc
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

NPROCS=$SLURM_NPROCS

FASTQC=/LUSTRE/apps/bioinformatica/FastQC_v0.12.1/
export PATH=$PATH:$FASTQC

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/


mkdir -p FASTQC_DIR
mkdir -p FASTQC_DIR/MULTIQC_DIR
for i in $(ls *.fastq.gz);
do
bs="${i%.fastq.gz}"

if [ ! -f "FASTQC_DIR/${bs}_fastqc.zip" ]; then

    call="fastqc $i -t $NPROCS --nogroup -o FASTQC_DIR"

    echo "Running fastqc in '${bs}' sample."

    echo $call

    eval $call

else
    echo "'FASTQC_DIR/${bs}_fastqc.zip' already exists."
    echo "Continue with next sample."
fi

done

$WDM/multiqc $PWD/FASTQC_DIR -o $PWD/FASTQC_DIR/MULTIQC_DIR --config multiqc_info.conf

exit