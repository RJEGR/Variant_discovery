# Author
Ricardo Gomez-Reyes
# Description
The intention for this lab is to perform variant discovery using GATK for model and non model species. The GATK Best Practices provide a detailed, well-tested workflow for germline short variant discovery in DNA sequencing data and RNAseq short variant per-sample calling. By following this pipeline, you can identify high-quality SNPs and indels from sequencing data aligned to a reference genome.

According to recent benchmarks, new methods show best performance in detect disc variant (Barbitoff,et al., 2022). Including deepvariant software (https://github.com/google/deepvariant/blob/r1.6.1/docs/deepvariant-rnaseq-case-study.md), but complicate to install in a server without root access


Barbitoff,et al., 2022 https://doi.org/10.1186/s12864-022-08365-3
# Reference details
## Features of GRCh38/hg38
If you are working with human data, we recommend you use one of the reference genome builds that we provide in our Resource Bundle or in Terra, our cloud-based analysis portal. We currently support GRCh38/hg38 and b37 (and to a lesser extent, hg19). For more [information](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19) on the human genome reference builds. We strongly recommend switching to GRCh38/hg38 if you are working with human sequence data. In addition to adding many alternate contigs, GRCh38 corrects thousands of small sequencing artifacts that cause false SNPs and indels to be called when using the GRCh37 assembly (b37/Hg19).

Also considere use alternate or ALT contigs to represent common complex variation, including HLA loci (Human leukocyte antigen, some HLA-mediated diseases are directly involved in the promotion of cancer). Alternate haplotypes are highly important to SNPs variants between [populations](https://www.perplexity.ai/search/what-is-a-human-haplotype-nuZaaXQFS9mrTWHhqb.WZg), for examples Human haplotypes blocks and HapMap

Perplexity summarises how to perform variant discovery:

**1. Prepare Input Files**
- Quality Control: Begin with quality control of raw sequencing data using tools like FastQC to assess the quality metrics.
- Align sequencing reads to the human reference genome using a tool like BWA.
- Sort and index the resulting BAM file.
- Add read group information to the BAM file using Picard tools.
- Index the human reference genome using samtools and create a sequence dictionary using Picard.
- Post-Alignment Processing: Use tools like Picard to sort and mark duplicates in the aligned BAM files. This helps to reduce biases in variant calling.
2. Perform Base Quality Score Recalibration (BQSR)
- Run GATK BaseRecalibrator to analyze patterns of covariation in the sequence dataset.
- Apply the recalibration using GATK ApplyBQSR to generate a recalibrated BAM file.
3. Call Variants
- Run GATK HaplotypeCaller on the recalibrated BAM file to call genetic variants.
- This performs local de novo assembly of haplotypes in regions with variation.
- Emits both SNPs and indels simultaneously via local re-assembly of haplotypes.
- Joint Genotyping: For multi-sample analyses, combine gVCF files from individual samples using CombineGVCFs and perform joint genotyping with GenotypeGVCFs.
4. Filter Variants
- Apply variant quality score recalibration (VQSR) using GATK VariantRecalibrator and ApplyVQSR.
- Builds a Gaussian mixture model from the annotations of a set of trusted variants.
- Assigns well-calibrated probability estimates to variant calls.
- Alternatively, apply hard filters based on specific criteria like quality, depth, strand bias, etc.
5. Annotate Variants
- Annotate the filtered VCF file with additional information from databases like dbSNP, ClinVar, etc.
- Variant Annotation: Annotate the filtered variants using tools like ANNOVAR or SnpEff to provide biological context, such as functional impacts and associations with known variants.