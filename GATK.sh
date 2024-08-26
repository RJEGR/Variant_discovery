#!/bin/sh
## Directivas
#SBATCH --job-name=GATK
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --propagate=STACK

ulimit -s unlimited

EXPORT=/LUSTRE/apps/bioinformatica/bwa/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/gatk-4.6.0.0
export PATH=$PATH:$EXPORT


mem=100 # gb

thread_count=$SLURM_NPROCS

#reference=$1

reference=`basename $1`

REF_PREFIX=${reference%.f*}

mkdir -p BWA_index

if [ ! -f "BWA_index/${REF_PREFIX}.bwt" ]; then
        bwa index -p BWA_index/$REF_PREFIX $reference
        
        ln -s $PWD/$reference BWA_index
else
    echo "BWA index for '${REF_PREFIX}' already exists."
fi

if [ ! -f "${REF_PREFIX}.dict" ]; then
    gatk CreateSequenceDictionary -R ${reference}
else
    echo "gatk dictionary for '${REF_PREFIX}' already exists."
fi

# Stage 1 (S1)

mkdir -p S1_BWA_SAM_BAM_FILES

for i in $(ls *_R1_001.fastq.gz);
do
bs="${i%_R1_001.fastq.gz}"

if [ ! -f "S1_BWA_SAM_BAM_FILES/${bs}.bam" ]; then

left_file=${bs}_R1_001.fastq.gz
right_file=${bs}_R2_001.fastq.gz

# because errors in  –R $read_group, omit and use AddOrReplaceReadGroups

read_group="'@RG\tID:${bs}\tSM:${bs}\tLB:${bs}\tPL:ILLUMINA'"

call="bwa mem -M -t $thread_count \
    BWA_index/${REF_PREFIX} \
    $left_file $right_file | \
    samtools view -bh -@ $thread_count -m 12G -o S1_BWA_SAM_BAM_FILES/${bs}.bam"
    
    # samtools view -bh -@ $thread_count -m 12G - | \
	# samtools sort -@ $thread_count -m 12G -O bam -o S1_BWA_SAM_BAM_FILES/${bs}.sorted.bam

    # using picard version for sort

echo $call;eval $call

else
    echo "Alignment for 'S1_BWA_SAM_BAM_FILES/${bs}.sorted.bam' already exists."
fi

# continue w/ https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38#3

# Step 2, pre-process bam for gatk


gatk --java-options "-Xmx7g" AddOrReplaceReadGroups \
    --I S1_BWA_SAM_BAM_FILES/${bs}.bam \
    --O S1_BWA_SAM_BAM_FILES/${bs}.RG.bam \
    --RGID ${bs} --RGSM ${bs} --RGLB ${bs} --RGPU NOVOGEN --RGPL ILLUMINA

mkdir -p S2_GATK_DIR

gatk SortSam \
    --I S1_BWA_SAM_BAM_FILES/${bs}.RG.bam \
    --O S2_GATK_DIR/${bs}.sorted.bam \
    --VALIDATION_STRINGENCY LENIENT \
    --SO coordinate \
    --MAX_RECORDS_IN_RAM 3000000 \
    --CREATE_INDEX true    

call="samtools flagstat S2_GATK_DIR/${bs}.sorted.bam > S2_GATK_DIR/${bs}_flagstat.txt"

echo $call;eval $call

gatk MarkDuplicates \
 	        --INPUT S2_GATK_DIR/${bs}.sorted.bam \
 	        --OUTPUT S2_GATK_DIR/${bs}.sorted.dup.bam  \
 	        --METRICS_FILE S2_GATK_DIR/${bs}.MarkDuplicates.metrics

# step 1  - Build the model
WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Human/hg38/ftp.broadinstitute.org/bundle/hg38/

gatk BaseRecalibrator \
    -I S2_GATK_DIR/${bs}.sorted.dup.bam \
    -R $reference \
    --known-sites ${WD}/dbsnp_146.hg38.vcf.gz \
    -O S2_GATK_DIR/${bs}_recal_data.table

# Step 2: Apply the model to adjust the base quality scores

gatk ApplyBQSR \
    -I S2_GATK_DIR/${bs}.sorted.dup.bam \
    -R $reference \
    --bqsr-recal-file S2_GATK_DIR/${bs}_recal_data.table \
    -O S2_GATK_DIR/${bs}.sort.dup.bqsr.bam

# And collect metrics

gatk CollectMultipleMetrics \
    -R $reference \
    -I S2_GATK_DIR/${bs}.sort.dup.bqsr.bam \
    -O S2_GATK_DIR/${bs}.sort.dup.bqsr.CollectMultipleMetrics

echo "The sample group ${bs} has been pre-processed\n" 
echo "BAM file S2_GATK_DIR/${bs}.sort.dup.bqsr.bam is ready for variant calling"

# Here you can run multiqc in S2_GATK_DIR

# Step 3: Apply HaplotypeCaller
# HaplotypeCaller is the focal tool within GATK4 to simultaneously call germline SNVs and small Indels using local de novo assembly of haplotype regions.
mkdir -p S3_GATK_Haplotype_DIR

gatk --java-options "-Xmx7g" HaplotypeCaller \
    -I S2_GATK_DIR/${bs}.sort.dup.bqsr.bam \
    -R $reference \
    -ERC GVCF \
    -O S3_GATK_Haplotype_DIR/${bs}.g.vcf.gz


done

# One or more VCF files containing variants  This argument must be specified at least once. Require

vcf_files=`ls S3_GATK_Haplotype_DIR/*g.vcf.gz`
vcf_files=`echo -V $vcf_files`

gatk --java-options "-Xmx7g" CombineGVCFs \
    -R $reference \
    $vcf_files \
    -O S3_GATK_Haplotype_DIR/CombineGVCF.g.vcf.gz

exit

