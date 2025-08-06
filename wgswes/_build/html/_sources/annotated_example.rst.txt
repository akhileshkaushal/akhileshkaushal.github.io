Annotated Examples & Interpretation
====================================

This section provides step-by-step annotated examples of WGS/WES workflows for **tumor–normal**, **tumor-only**, and **germline** analyses, along with notes on **result interpretation** and **clinical caveats**.

.. contents::
   :local:
   :depth: 2

Example: Tumor–Normal Variant Calling
-------------------------------------

**Scenario:**  
Paired tumor–normal workflow using WES data (extendable to WGS by adjusting depth, PoN selection, and compute resources).

**Inputs:**
- **Tumor**: `tumor_R1.fq.gz`, `tumor_R2.fq.gz`
- **Normal**: `normal_R1.fq.gz`, `normal_R2.fq.gz`
- **Reference**: GRCh38/hg38, associated indices, and known sites
- **Resources**: PoN VCF, germline resource (gnomAD AF), exome BED

**Steps:**

1. **QC and Trimming**

   .. code-block:: bash

      fastp \
        -i tumor_R1.fq.gz -I tumor_R2.fq.gz \
        -o tumor_trim_R1.fq.gz -O tumor_trim_R2.fq.gz \
        --html tumor_fastp.html --thread 16

      fastp \
        -i normal_R1.fq.gz -I normal_R2.fq.gz \
        -o normal_trim_R1.fq.gz -O normal_trim_R2.fq.gz \
        --html normal_fastp.html --thread 16

2. **Alignment**

   .. code-block:: bash

      bwa-mem2 mem -t 16 ref.fa tumor_trim_R1.fq.gz tumor_trim_R2.fq.gz \
        | samtools sort -o tumor.bam

      bwa-mem2 mem -t 16 ref.fa normal_trim_R1.fq.gz normal_trim_R2.fq.gz \
        | samtools sort -o normal.bam

3. **Duplicate Marking & BQSR**

   .. code-block:: bash

      gatk MarkDuplicates -I tumor.bam -O tumor.markdup.bam -M tumor.metrics.txt
      gatk MarkDuplicates -I normal.bam -O normal.markdup.bam -M normal.metrics.txt

      gatk BaseRecalibrator -I tumor.markdup.bam -R ref.fa \
        --known-sites dbsnp.vcf --known-sites mills.vcf -O tumor.recal.table
      gatk BaseRecalibrator -I normal.markdup.bam -R ref.fa \
        --known-sites dbsnp.vcf --known-sites mills.vcf -O normal.recal.table

      gatk ApplyBQSR -R ref.fa -I tumor.markdup.bam \
        --bqsr-recal-file tumor.recal.table -O tumor.bqsr.bam
      gatk ApplyBQSR -R ref.fa -I normal.markdup.bam \
        --bqsr-recal-file normal.recal.table -O normal.bqsr.bam

4. **Somatic Calling**

   .. code-block:: bash

      gatk Mutect2 \
        -R ref.fa \
        -I tumor.bqsr.bam -tumor TUMOR \
        -I normal.bqsr.bam -normal NORMAL \
        --germline-resource af-only-gnomad.vcf.gz \
        --panel-of-normals pon.vcf.gz \
        -O somatic.vcf.gz

5. **Filtering**

   .. code-block:: bash

      gatk FilterMutectCalls \
        -V somatic.vcf.gz -R ref.fa \
        -O somatic.filtered.vcf.gz

6. **Annotation**

   .. code-block:: bash

      vep \
        -i somatic.filtered.vcf.gz \
        -o somatic.annotated.vcf \
        --cache --everything --vcf --assembly GRCh38 --offline

**Interpretation Notes:**
- Variants labeled **PASS** are high-confidence somatic calls.
- Check the ``FILTER`` field for failure reasons (e.g., ``contamination``, ``panel_of_normals``).
- **VAF:** Low VAF (<5%) may be subclonal; interpret with tumor purity.
- Consequence terms (missense, nonsense, frameshift) inform functional impact.
- Cross-reference COSMIC/OncoKB for known drivers; **use orthogonal validation** before any clinical decision-making.

Example: Tumor-Only Variant Calling
-----------------------------------

**Scenario:**  
No matched normal available—requires aggressive artifact and germline suppression.

**Inputs:**
- **Tumor**: `tumor_R1.fq.gz`, `tumor_R2.fq.gz`
- PoN VCF, gnomAD germline resource

**Steps:**

1. **QC and Trimming**

   .. code-block:: bash

      fastp \
        -i tumor_R1.fq.gz -I tumor_R2.fq.gz \
        -o tumor_trim_R1.fq.gz -O tumor_trim_R2.fq.gz \
        --html tumor_fastp.html --thread 16

2. **Alignment**

   .. code-block:: bash

      bwa-mem2 mem -t 16 ref.fa tumor_trim_R1.fq.gz tumor_trim_R2.fq.gz \
        | samtools sort -o tumor.bam

3. **Duplicate Marking & BQSR**

   .. code-block:: bash

      gatk MarkDuplicates -I tumor.bam -O tumor.markdup.bam -M tumor.metrics.txt

      gatk BaseRecalibrator -I tumor.markdup.bam -R ref.fa \
        --known-sites dbsnp.vcf --known-sites mills.vcf \
        -O tumor.recal.table

      gatk ApplyBQSR -R ref.fa -I tumor.markdup.bam \
        --bqsr-recal-file tumor.recal.table -O tumor.bqsr.bam

4. **Somatic Calling (Tumor-Only Mode)**

   .. code-block:: bash

      gatk Mutect2 \
        -R ref.fa -I tumor.bqsr.bam -tumor TUMOR \
        --germline-resource af-only-gnomad.vcf.gz \
        --panel-of-normals pon.vcf.gz \
        -O tumoronly.vcf.gz

5. **Filtering**

   .. code-block:: bash

      gatk FilterMutectCalls \
        -V tumoronly.vcf.gz -R ref.fa \
        -O tumoronly.filtered.vcf.gz

6. **Annotation**

   .. code-block:: bash

      vep \
        -i tumoronly.filtered.vcf.gz \
        -o tumoronly.annotated.vcf \
        --cache --everything --vcf --assembly GRCh38 --offline

**Interpretation Notes:**
- Expect higher **false positives**; rare germline often remains even after filtering.
- Use gnomAD AF thresholds (e.g., exclude AF > 0.001 for rare cancers).
- Orthogonal validation (Sanger, ddPCR) is **strongly recommended**.
- Review strand bias and read-level evidence in IGV for low-VAF calls.

Example: Germline Variant Calling
---------------------------------

**Scenario:**  
Single-sample germline variant discovery using GATK Best Practices (extendable to joint calling for cohorts).

**Inputs:**
- **Sample**: `sample_R1.fq.gz`, `sample_R2.fq.gz`
- Known sites (dbSNP, Mills, 1000G indels)

**Steps:**

1. **QC and Trimming**

   .. code-block:: bash

      fastp \
        -i sample_R1.fq.gz -I sample_R2.fq.gz \
        -o trim_R1.fq.gz -O trim_R2.fq.gz \
        --html sample_fastp.html --thread 16

2. **Alignment**

   .. code-block:: bash

      bwa-mem2 mem -t 16 ref.fa trim_R1.fq.gz trim_R2.fq.gz \
        | samtools sort -o sample.bam

3. **Duplicate Marking**

   .. code-block:: bash

      gatk MarkDuplicates -I sample.bam -O sample.markdup.bam -M sample.metrics.txt

4. **BQSR**

   .. code-block:: bash

      gatk BaseRecalibrator -I sample.markdup.bam -R ref.fa \
        --known-sites dbsnp.vcf --known-sites mills.vcf \
        -O sample.recal.table

      gatk ApplyBQSR -R ref.fa -I sample.markdup.bam \
        --bqsr-recal-file sample.recal.table \
        -O sample.bqsr.bam

5. **Variant Calling**

   .. code-block:: bash

      gatk HaplotypeCaller -R ref.fa -I sample.bqsr.bam \
        -O sample.g.vcf.gz -ERC GVCF

6. **Joint Genotyping (multi-sample)**

   .. code-block:: bash

      gatk CombineGVCFs -R ref.fa \
        -V sample1.g.vcf.gz \
        -V sample2.g.vcf.gz \
        -O cohort.g.vcf.gz

      gatk GenotypeGVCFs -R ref.fa -V cohort.g.vcf.gz \
        -O cohort.vcf.gz

7. **Filtering**

   Large cohorts: VQSR  
   Small cohorts: hard filters

   Example hard filters:

   .. code-block:: bash

      gatk VariantFiltration -V cohort.vcf.gz \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        -O cohort.filtered.vcf.gz

8. **Annotation**

   .. code-block:: bash

      vep -i cohort.filtered.vcf.gz \
          -o cohort.annotated.vcf \
          --cache --everything --vcf --assembly GRCh38 --offline

**Interpretation Notes:**
- Filter on quality metrics (``QD``, ``FS``, ``MQ``) to reduce false positives.
- Use gnomAD to flag common benign variants.
- Classify by ACMG/AMP where appropriate (research context): **pathogenic**, **likely pathogenic**, **VUS**, **likely benign**, **benign**.
- For research publications, emphasize **novel**, **rare**, and **functionally significant** variants.

.. note::
   Regardless of workflow type, validate high-impact variants with orthogonal methods before clinical use, and interpret in the context of phenotype, tumor histology (if applicable), and supporting molecular data.
