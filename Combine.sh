#!/bin/sh
## Directivas
#SBATCH --job-name=CombineGVCFs
#SBATCH -N 1
#SBATCH --mem=50GB
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

#vcf_files=`ls S4_GATK_Mutect2_DIR/*.f.vcf.gz`

ls -1 S4_GATK_Mutect2_DIR/*.f.vcf.gz > MergeList

vcf_files=`awk '{print "-V",$1}' MergeList`

call="gatk --java-options "-Xmx40g" CombineGVCFs -R $reference $vcf_files -O S4_GATK_Mutect2_DIR/CombineGVCF.f.vcf.gz"

echo $call

eval $call

call="gatk --java-options "-Xmx40g" Funcotator \
	--variant S4_GATK_Mutect2_DIR/CombineGVCF.f.vcf.gz \
	--reference $reference \
	--ref-version hg38 \
	--data-sources-path $dataSourcesFolder/funcotator_dataSources.v1.8.hg38.20230908s \
	--output S4_GATK_Mutect2_DIR/CombineGVCF.f.functotated.maf.gz \
	--output-file-format MAF"
    
echo $call
eval $call

exit