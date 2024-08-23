#!/bin/sh
## Directivas
#SBATCH --job-name=gatk
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --propagate=STACK

# In addition to --propagate, 
# increase the max number of open files with Linux command ulimit before running STAR

# ulimit -s unlimited

# EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/STAR-2.7.11b/bin/Linux_x86_64_static/
# export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/bwa/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/gatk-4.6.0.0
export PATH=$PATH:$EXPORT


star_mem=100 # gb
thread_count=$SLURM_NPROCS
reference=$1

mkdir -p BWA_index

if [ ! -f "BWA_index/${reference%.f}.bwt" ]; then
        bwa index -p ${reference%.f} $reference
else
    echo "BWA index for '$reference' already exists."
fi

if [ ! -f "${reference%.*}.dict" ]; then
    gatk CreateSequenceDictionary -R ${reference}
else
    echo "gatk dir for '$reference' already exists."
fi

# 1)


for i in $(ls *1P.fq);
do
base_name="${i%_1P.fq}"

left_file=${base_name}_1P.fq
right_file=${base_name}_2P.fq

#read_group="@RG\tID:SRR622461.7\tSM:NA12878\tLB:ERR194147\tPL:ILLUMINA"
read_group="@RG\tID:${base_name}\tSM:${base_name}\tLB:${base_name}\tPL:ILLUMINA"


bwa mem -M -t $thread_count –R $read_group $reference $left_file $right_file > ${base_name}.sam
# if alt file
# bwa mem chr19chr19KI270866v1alt.fasta 8017read1.fq 8017read2.fq > 8017bwamem.sam
# continue w/ https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38#3
#SORTEDBAM=${base_name}.Aligned.sortedByCoord.out.bam

# Step 2, pre-process bam for gatk


gatk MergeBamAlignment \
            --REFERENCE_SEQUENCE ${reference} \
            --ALIGNED_BAM ${SORTEDBAM%.bam}.al.bam \
            --UNMAPPED_BAM ${SORTEDBAM%.bam}.uT.bam \
            --SORT_ORDER coordinate \
            --OUTPUT $PWD/${base_name}.merged.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT
 
            
done

exit

gatk MarkDuplicates \
 	        --INPUT ${base_name}.merged.bam \
 	        --OUTPUT ${base_name}.merged.dup.bam  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE ${base_name}.metrics