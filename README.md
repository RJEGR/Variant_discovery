# Author
Ricardo Gomez-Reyes
# Description
The intention for this lab is to perform variant discovery using GATK for model and non model species. The GATK Best Practices provide a detailed, well-tested workflow for germline short variant discovery in DNA sequencing data and RNAseq short variant per-sample calling. By following this pipeline, you can identify high-quality SNPs and indels from sequencing data aligned to a reference genome.


# 
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