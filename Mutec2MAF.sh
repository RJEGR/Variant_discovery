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

done

dataSourcesFolder=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Human/hg38/funcotator/

if [ ! -f "$dataSourcesFolder/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz" ]; then

    call="gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download --hg38"
    eval $call
else
    echo "FuncotatorDataSource already was downloaded"
fi

mkdir -p S5_GATK_SelectVariants_MAF_DIR

for i in $(ls S4_GATK_Mutect2_DIR/*.f.vcf.gz);
do
basename=${i##*/}
bs="${basename%.f.vcf.gz}"

if [ ! -f "S5_GATK_SelectVariants_MAF_DIR/${bs}.f.functotated.maf.gz" ]; then

    call="gatk --java-options "-Xmx20g" Funcotator \
	--variant S4_GATK_Mutect2_DIR/${bs}.f.vcf.gz \
	--reference $reference \
	--ref-version hg38 \
	--data-sources-path $dataSourcesFolder/funcotator_dataSources.v1.8.hg38.20230908s \
	--output S5_GATK_SelectVariants_MAF_DIR/${bs}.f.functotated.maf.gz \
	--output-file-format MAF"
    
    eval $call

else
    echo "Funcotator already was run for S5_GATK_SelectVariants_MAF_DIR/${bs}.f.functotated.maf.gz"
fi

done

exit