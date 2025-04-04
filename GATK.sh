#!/bin/sh
## Directivas
#SBATCH --job-name=GATK
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH -t 6-00:00:00

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
mkdir -p CHKPNT_DIR

mkdir -p S1_BWA_SAM_BAM_FILES

PATTERN="_R1.fq.gz"

for i in $(ls *${PATTERN});
do
bs="${i%$PATTERN}"

if [ ! -f "CHKPNT_DIR/${bs}_bwaMem.chkpt" ]; then

left_file=${bs}_R1.fq.gz
right_file=${bs}_R2.fq.gz

# because errors in  –R $read_group, omit and use AddOrReplaceReadGroups

read_group="'@RG\tID:${bs}\tSM:${bs}\tLB:${bs}\tPL:ILLUMINA'"

call="bwa mem -M -t $thread_count \
    BWA_index/${REF_PREFIX} \
    $left_file $right_file | \
    samtools view -bh -@ $thread_count -m 12G -o S1_BWA_SAM_BAM_FILES/${bs}.bam"
    
    # samtools view -bh -@ $thread_count -m 12G - | \
	# samtools sort -@ $thread_count -m 12G -O bam -o S1_BWA_SAM_BAM_FILES/${bs}.sorted.bam

    # using picard version for sort

echo $call

eval $call

unlink $left_file

unlink $right_file

touch  CHKPNT_DIR/${bs}_bwaMem.chkpt

else
    echo "Alignment for 'S1_BWA_SAM_BAM_FILES/${bs}.sorted.bam' already exists."
fi

# continue w/ https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38#3

# Step 2, pre-process bam for gatk

# Is working now

if [ ! -f "S1_BWA_SAM_BAM_FILES/${bs}.RG.bam" ]; then

gatk --java-options "-Xmx20g" AddOrReplaceReadGroups \
    --I S1_BWA_SAM_BAM_FILES/${bs}.bam \
    --O S1_BWA_SAM_BAM_FILES/${bs}.RG.bam \
    --RGID ${bs} --RGSM ${bs} --RGLB ${bs} --RGPU NOVOGEN --RGPL ILLUMINA

else
    echo "file already exists. Omitting AddOrReplaceReadGroups"
fi

mkdir -p S2_GATK_DIR


if [ ! -f "S2_GATK_DIR/${bs}.sorted.bam" ]; then

    gatk SortSam \
        --I S1_BWA_SAM_BAM_FILES/${bs}.RG.bam \
        --O S2_GATK_DIR/${bs}.sorted.bam \
        --VALIDATION_STRINGENCY LENIENT \
        --SO coordinate \
        --MAX_RECORDS_IN_RAM 3000000 \
        --CREATE_INDEX true    


else
    echo "file already exists. Omitting SortSam"
fi


call="samtools flagstat S2_GATK_DIR/${bs}.sorted.bam > S2_GATK_DIR/${bs}_flagstat.txt"


if [ ! -f "S2_GATK_DIR/${bs}_flagstat.txt" ]; then

    echo $call
    eval $call

else
    echo "file already exists. Omitting samtools flagstat"
fi


if [ ! -f "S2_GATK_DIR/${bs}.sorted.dup.bam" ]; then

    gatk MarkDuplicates \
 	    --INPUT S2_GATK_DIR/${bs}.sorted.bam \
 	    --OUTPUT S2_GATK_DIR/${bs}.sorted.dup.bam  \
 	    --METRICS_FILE S2_GATK_DIR/${bs}.MarkDuplicates.metrics

else
    echo "file already exists. Omitting MarkDuplicates"
fi


# step 1  - Build the model

WD=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Human/hg38/ftp.broadinstitute.org/bundle/hg38/

if [ ! -f "S2_GATK_DIR/${bs}_recal_data.table" ]; then

    gatk BaseRecalibrator \
        -I S2_GATK_DIR/${bs}.sorted.dup.bam \
        -R $reference \
        --known-sites ${WD}/dbsnp_146.hg38.vcf.gz \
        -O S2_GATK_DIR/${bs}_recal_data.table

else
    echo "file already exists. Omitting BaseRecalibrator"
fi 

# Step 2: Apply the model to adjust the base quality scores

if [ ! -f "S2_GATK_DIR/${bs}.sort.dup.bqsr.bam" ]; then

    gatk ApplyBQSR \
        -I S2_GATK_DIR/${bs}.sorted.dup.bam \
        -R $reference \
        --bqsr-recal-file S2_GATK_DIR/${bs}_recal_data.table \
        -O S2_GATK_DIR/${bs}.sort.dup.bqsr.bam

else
    echo "file already exists. Omitting ApplyBQSR"
fi 


# And collect metrics

if [ ! -f "S2_GATK_DIR/${bs}.sort.dup.bqsr.CollectMultipleMetrics" ]; then

    gatk CollectMultipleMetrics \
        -R $reference \
        -I S2_GATK_DIR/${bs}.sort.dup.bqsr.bam \
        -O S2_GATK_DIR/${bs}.sort.dup.bqsr.CollectMultipleMetrics

else
    echo "file already exists. Omitting CollectMultipleMetrics"
fi 

touch  CHKPNT_DIR/${bs}_ApplyBQSR.chkpt

echo "The sample group ${bs} has been pre-processed\n" 
echo "BAM file S2_GATK_DIR/${bs}.sort.dup.bqsr.bam is ready for variant calling"

done

# Alternative step:

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/
mkdir -p MULTIQC_VIZ_DIR
$WDM/multiqc $PWD/S1_BWA_SAM_BAM_FILES/ $PWD/S2_GATK_DIR -o MULTIQC_VIZ_DIR --config multiqc_info.conf


# Step 3: Apply HaplotypeCaller
# HaplotypeCaller is the focal tool within GATK4 to simultaneously call germline SNVs and small Indels using local de novo assembly of haplotype regions.
mkdir -p S3_GATK_Haplotype_DIR

echo "######################################################\n"
echo "Running HaplotypeCaller"
echo "######################################################\n"

for i in $(ls S2_GATK_DIR/*.sort.dup.bqsr.bam);
do
basename=${i##*/}
bs="${basename%.sort.dup.bqsr.bam}"

if [ ! -f "CHKPNT_DIR/${bs}_HaplotypeCaller.chkpt" ]; then

    call="gatk --java-options "-Xmx50g" HaplotypeCaller \
        --native-pair-hmm-threads $thread_count \
        -I S2_GATK_DIR/${bs}.sort.dup.bqsr.bam \
        -R $reference \
        -ERC GVCF \
        -O S3_GATK_Haplotype_DIR/${bs}.g.vcf.gz"


    echo $call
    eval $call

    touch  CHKPNT_DIR/${bs}_HaplotypeCaller.chkpt

else
    echo "file already exists. Omitting HaplotypeCaller"
fi 

echo "HaplotypeCaller was done for $bs group"

done

echo "gatk analyses for this directory were done"
echo "continue with filter and annotate variant step"


exit

# One or more VCF files containing variants  This argument must be specified at least once. Require

vcf_files=`ls S3_GATK_Haplotype_DIR/*g.vcf.gz`
vcf_files=`echo -V $vcf_files`

gatk --java-options "-Xmx7g" CombineGVCFs \
    -R $reference \
    $vcf_files \
    -O S3_GATK_Haplotype_DIR/CombineGVCF.g.vcf.gz

exit


vcf_files=`ls S4_GATK_Mutect2_DIR/*.f.vcf.gz`
vcf_files=`echo -V $vcf_files`

call="gatk --java-options "-Xmx7g" CombineGVCFs \
    -R $reference \
    $vcf_files \
    -O S4_GATK_Mutect2_DIR/CombineGVCF.f.vcf.gz"
