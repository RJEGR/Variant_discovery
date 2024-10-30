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

if [ ! -f "S4_GATK_Mutect2_DIR/${bs}.f.vcf.gz" ]; then

    call="gatk FilterMutectCalls -R $reference -V S4_GATK_Mutect2_DIR/${bs}.vcf.gz -O S4_GATK_Mutect2_DIR/${bs}.f.vcf.gz"

    echo $call
    eval $call

else
    echo "Filtered file already exists. Omitting FilterMutectCalls"
fi

echo "Mutect2 and FilterMutectCalls was done for $bs group"

# done

# exit

# extract SNPs & INDELS

# SNPs

mkdir -p S5_GATK_SelectVariants_DIR

if [ ! -f "S5_GATK_SelectVariants_DIR/${bs}.f.snp.vcf.gz" ]; then

    call="gatk SelectVariants -R $reference -V S4_GATK_Mutect2_DIR/${bs}.f.vcf.gz --select-type SNP -O S5_GATK_SelectVariants_DIR/${bs}.f.snp.vcf.gz"

    echo $call
    eval $call

else
    echo "Filtered file already exists. Omitting SelectVariants (SNPs) and Funcotator"
fi

# INDELS

if [ ! -f "S5_GATK_SelectVariants_DIR/${bs}.f.indels.vcf.gz" ]; then

    call="gatk SelectVariants -R $reference -V S4_GATK_Mutect2_DIR/${bs}.f.vcf.gz --select-type INDEL -O S5_GATK_SelectVariants_DIR/${bs}.f.indels.vcf.gz"

    echo $call
    eval $call

else
    echo "Filtered file already exists. Omitting SelectVariants (Indels) and Funcotator"
fi

# Annotation

dataSourcesFolder=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Human/hg38/funcotator/funcotator_dataSources.v1.8.hg38.20230908s

if [ ! -f "$dataSourcesFolder/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz" ]; then

    call="gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download --hg38"
    eval $call
else
    echo "FuncotatorDataSource already was downloaded"
fi

# Choose Variant type
variant_type=indels

if [ ! -f "S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.vcf.gz" ]; then

    call="gatk Funcotator \
	--variant S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.vcf.gz \
	--reference $reference \
	--ref-version hg38 \
	--data-sources-path $dataSourcesFolder \
	--output S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.vcf.gz \
	--output-file-format VCF"
    
    eval $call

    call="gatk VariantsToTable -V S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.vcf.gz -F AC -F AN -F DP -F AF -F FUNCOTATION -O S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.table"

    eval $call
else
    echo "Funcotator already was run for ${bs}.f.${variant_type}.functotated.vcf.gz"
fi


variant_type=snp

if [ ! -f "S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.vcf.gz" ]; then

    call="gatk Funcotator \
	--variant S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.vcf.gz \
	--reference $reference \
	--ref-version hg38 \
	--data-sources-path $dataSourcesFolder \
	--output S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.vcf.gz \
	--output-file-format VCF"
    
    eval $call

    call="gatk VariantsToTable -V S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.vcf.gz -F AC -F AN -F DP -F AF -F FUNCOTATION -O S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.functotated.table"

    eval $call


    
else
    echo "Funcotator already was run for ${bs}.f.${variant_type}.functotated.vcf.gz"
fi

done

exit
