# Prognostic insights for breast cancer. 
# Author
Ricardo Gomez-Reyes
# Background
<div align="justify">
Breast cancer is characterized by various genomic disruptions, including altered gene expression and splicing, gene mutation. 
These genomic and transcriptomic disruptions often interact and contribute to the hallmarks of breast cancer, including sustained proliferation, evasion of growth suppressors, resistance to cell death, and activation of invasion and metastasis.
</div>


| Genomic disruption    | Common |
| -------- | ------- |
| Altered Gene Expression  | Overexpression of oncogenes (e.g., MYC, CCND1), and Downregulation of tumor suppressors (e.g., TP53, PTEN)     |
| Altered splicing of genes  | Genes like CD44, BRCA1, and TP53.     |
| Common Gene Mutations (SNPs)    | - TP53 (tumor suppressor gene), PIK3CA (oncogene in PI3K pathway), BRCA1 and BRCA2 (DNA repair genes), PTEN (tumor suppressor gene) and AKT1, GATA3, CDH1 (various roles in cell processes)  |
| Copy Number Variations (Gene Insertion/Deletions) | Amplification of HER2/ERBB2 (defines HER2+ subtype), MYC, CCND1, and FGFR1, Deletion of RB1, PTEN, and CDKN2A |
|  |  |



# Bioinformatic insights

<div align="justify">

The intention for this insilico lab is to perform variant discovery in breast-cancer samples analysis using tools such as GATK. The GATK Best Practices provide a detailed, well-tested workflow for germline short variant discovery in DNA sequencing data and RNAseq short variant per-sample calling. By following this pipeline, you can identify high-quality SNPs and indels from sequencing data aligned to a reference genome.

According to recent benchmarks, new methods show best performance in discover variant (Barbitoff,et al., 2022). Including [deepvariant](https://github.com/google/deepvariant/blob/r1.6.1/docs/deepvariant-rnaseq-case-study.md) software, but complicate to install in a server without root access

</div>

Barbitoff,et al., 2022 https://doi.org/10.1186/s12864-022-08365-3
# Reference details
## Features of GRCh38/hg38
For human data, the GATK protocol author recommend to: 
- Use one of the reference genome builds that gatk provide in their Resource Bundle or in Terra, their cloud-based analysis portal. The authors currently support GRCh38/hg38 and b37 (and to a lesser extent, hg19). For more [information](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19) on the human genome reference builds. 

- Switching to GRCh38/hg38 if you are working with human sequence data (Strongly recommend). In addition to adding many alternate contigs, GRCh38 corrects thousands of small sequencing artifacts that cause false SNPs and indels to be called when using the GRCh37 assembly (b37/Hg19).

- Considere use alternate or ALT contigs to represent common complex variation, including HLA loci (Human leukocyte antigen, some HLA-mediated diseases are directly involved in the promotion of cancer). Alternate haplotypes are highly important to SNPs variants between [populations](https://www.perplexity.ai/search/what-is-a-human-haplotype-nuZaaXQFS9mrTWHhqb.WZg), for examples Human haplotypes blocks and HapMap.

## GATK bundle of reference
This is a truth sets of variants for hg38 reference needed for short variant discovery in WGS. It can be found in [resource boundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) section, and can be download easly by:

```bash
wget -m ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/
```

## RNAseq short variant discovery (SNPs + Indels)
In contrast to other alignment methods, this protocol recommend to use STAR aligner because it increased sensitivity (especially for INDELS). This [protocol](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels) use STAR’s two-pass mode to get better alignments around novel splice junctions.



Finally, IA perplexity will help us to summarises how to perform variant discovery:

**1. Prepare Input Files**
- Quality Control: Begin with quality control of raw sequencing data using tools like FastQC to assess the quality metrics.
- Align sequencing reads to the human reference genome using a tool like BWA.
- Sort and index the resulting BAM file.
- Add read group information to the BAM file using Picard tools.
- Index the human reference genome using samtools and create a sequence dictionary using Picard.
- Post-Alignment Processing: Use tools like Picard to sort and mark duplicates in the aligned BAM files. This helps to reduce biases in variant calling.

**2. Perform Base Quality Score Recalibration (BQSR)**
- Run GATK BaseRecalibrator to analyze patterns of covariation in the sequence dataset.
- Apply the recalibration using GATK ApplyBQSR to generate a recalibrated BAM file.

**3. Call Variants**
- Run GATK HaplotypeCaller on the recalibrated BAM file to call genetic variants.
- This performs local de novo assembly of haplotypes in regions with variation.
- Emits both SNPs and indels simultaneously via local re-assembly of haplotypes.
- Joint Genotyping: For multi-sample analyses, combine gVCF files from individual samples using CombineGVCFs and perform joint genotyping with GenotypeGVCFs.

**4. Filter Variants**
- Apply variant quality score recalibration (VQSR) using GATK VariantRecalibrator and ApplyVQSR.
- Builds a Gaussian mixture model from the annotations of a set of trusted variants.
- Assigns well-calibrated probability estimates to variant calls.
- Alternatively, apply hard filters based on specific criteria like quality, depth, strand bias, etc.

**5. Annotate Variants**
- Annotate the filtered VCF file with additional information from databases like dbSNP, ClinVar, etc.
- Variant Annotation: Annotate the filtered variants using tools like ANNOVAR or SnpEff to provide biological context, such as functional impacts and associations with known variants.

# 1) Quality control
Prior to perform the variant discovery analysis, lets to screen data by fastqc, and preprocessin since it is hoghly recommended ([Pfeifer, 2017](https://doi.org/10.1038/hdy.2016.102))
```bash
FASTQC=/LUSTRE/apps/bioinformatica/FastQC_v0.12.1/
export PATH=$PATH:$FASTQC

mkdir -p fastqc_raw

srun fastqc *.fastq.gz -t 20 --nogroup -o fastqc_raw &> fastqc.log &

WDM=/LUSTRE/apps/Anaconda/2023/miniconda3/bin/

$WDM/multiqc $PWD/fastqc_raw -o $PWD/fastqc_raw
```

# 2) Base Quality Score Recalibration (BQSR)
I wrote a code to easly perform well-tested workflow for variant calling discovery using GATK Best Practices, and matching the experience of the Melbourne Bioinformatics, University of Melbourne ([Mahmood K, 2024](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/#variant-calling-using-gatk4)).  

```bash
sbatch GATK.sh Homo_sapiens_assembly38.fasta
```

# 3) Call Variants
## HaplotypeCaller vs Mutect2

HaplotypeCaller and Mutect2 are both tools from the Genome Analysis Toolkit (GATK) used for variant calling, but they are designed for different types of variants and have distinct methodologies. In summary, the choice between HaplotypeCaller and Mutect2 should be guided by the nature of the samples and the type of variants being studied: use HaplotypeCaller for germline variants and Mutect2 for somatic variants in cancer research. 

For calling variants in breast cancer samples, the recommended strategy is to use both HaplotypeCaller and Mutect2 from the GATK toolkit, as they serve different purposes. By combining the outputs of HaplotypeCaller and Mutect2, you can obtain a comprehensive list of germline predisposition variants and somatic mutations specific to the breast tumor, enabling a more complete understanding of the genetic landscape of the cancer. For example, Shin et al., 2024 describes Variant calling for germline small variants utilized HaplotypeCaller and Strelka2, while somatic small variant detection employed Strelka2 and Mutect2 ([Shin et al., 2024]([text](https://doi.org/10.1159/000536087)))

|               | HaplotypeCaller | Mutect2 |
| :---------------: | ------: | ----: |
| Purpose        |   HaplotypeCaller is primarily used for calling germline variants. It is optimized for identifying single nucleotide variants (SNVs) and insertions/deletions (indels) in diploid organisms, assuming a diploid model for variant calling.   | Mutect2 is specifically designed for calling somatic variants, primarily in cancer genomics. It focuses on identifying mutations that occur in tumor cells but are not present in normal cells.|
| Methodology           |   It employs a haplotype-based approach, utilizing local assembly of haplotypes and pair-HMM (Hidden Markov Model) alignment to improve the accuracy of variant calls. HaplotypeCaller can generate GVCF (Genomic VCF) files, which are useful for joint genotyping across multiple samples.ue   | Uses a similar local assembly and alignment approach as HaplotypeCaller but incorporates somatic-specific genotyping and filtering techniques. It can operate in a tumor-only mode, which allows for variant calling without a matched normal sample, although this is less supported than tumor-normal comparisons. |
| Use Cases    |  It is suitable for whole-genome or exome sequencing data where the focus is on identifying common and rare variants across populations or within individuals.   | It is optimized for detecting low-frequency somatic mutations, which are common in cancer due to the heterogeneous nature of tumor cells. Mutect2's algorithms account for factors like subclonality and copy number variations, making it more suitable for cancer studies. |

```bash
https://multiqc.info/modules/gatk/
```