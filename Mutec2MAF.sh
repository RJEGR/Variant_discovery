#!/bin/sh
## Directivas
#SBATCH --job-name=Funcotator
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH -t 6-00:00:00


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/gatk-4.6.0.0
export PATH=$PATH:$EXPORT

thread_count=$SLURM_NPROCS

reference=$1

dataSourcesFolder=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Human/hg38/funcotator/


if [ ! -f "$dataSourcesFolder/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz" ]; then

    call="gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download --hg38"
    eval $call
else
    echo "FuncotatorDataSource already was downloaded"
fi

mkdir -p S5_GATK_SelectVariants_MAF_DIR

for i in $(ls S2_GATK_DIR/*.sort.dup.bqsr.bam);
do
basename=${i##*/}
bs="${basename%.sort.dup.bqsr.bam}"

# Choose Variant type
variant_type=indels

if [ ! -f "S5_GATK_SelectVariants_MAF_DIR/${bs}.f.${variant_type}.functotated.vcf.gz" ]; then

    call="gatk Funcotator \
	--variant S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.vcf.gz \
	--reference $reference \
	--ref-version hg38 \
	--data-sources-path $dataSourcesFolder \
	--output S5_GATK_SelectVariants_MAF_DIR/${bs}.f.${variant_type}.functotated.maf.gz \
	--output-file-format MAF"
    
    eval $call

else
    echo "Funcotator already was run for ${bs}.f.${variant_type}.functotated.vcf.gz"
fi


variant_type=snp

if [ ! -f "S5_GATK_SelectVariants_MAF_DIR/${bs}.f.${variant_type}.functotated.vcf.gz" ]; then

    call="gatk Funcotator \
	--variant S5_GATK_SelectVariants_DIR/${bs}.f.${variant_type}.vcf.gz \
	--reference $reference \
	--ref-version hg38 \
	--data-sources-path $dataSourcesFolder \
	--output S5_GATK_SelectVariants_MAF_DIR/${bs}.f.${variant_type}.functotated.maf.gz \
	--output-file-format MAF"
    
    eval $call

else
    echo "Funcotator already was run for ${bs}.f.${variant_type}.functotated.vcf.gz"
fi

done

exit
