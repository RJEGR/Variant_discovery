#!/bin/sh
## Directivas
#SBATCH --job-name=Mutect2
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH -t 6-00:00:00

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/gatk-4.6.0.0
export PATH=$PATH:$EXPORT

thread_count=$SLURM_NPROCS

reference=$1

mkdir -p S4_GATK_Mutect2_DIR


for i in $(ls S2_GATK_DIR/*.sort.dup.bqsr.bam);
do
basename=${i##*/}
bs="${basename%.sort.dup.bqsr.bam}"

if [ ! -f "CHKPNT_DIR/${bs}_Mutect2.chkpt" ]; then

    call="gatk --java-options "-Xmx50g" Mutect2 \
           --native-pair-hmm-threads $thread_count \
           -R $reference \
           -I S2_GATK_DIR/${bs}.sort.dup.bqsr.bam \
           -O S4_GATK_Mutect2_DIR/${bs}.vcf.gz"

    echo $call
    eval $call

    touch  CHKPNT_DIR/${bs}_Mutect2.chkpt

else
    echo "file already exists. Omitting Mutect2"
fi 

echo "Mutect2 was done for $bs group"

done

exit