Advanced Analysis
=================

This section covers advanced workflows and analysis modules extending the standard WGS/WES pipeline. These include **family-based (trio/quad) analyses**, **joint-calling strategies for large cohorts**, **copy-number and mitochondrial genome analysis**, and **somatic tumor purity assessment and adjustment**.

.. contents::
   :local:
   :depth: 2

Trio and Family-Based Analysis
------------------------------

**Overview**

Family-based designs (e.g., trios: proband + parents; quads: proband + parents + sibling) enhance the detection of de novo variants, phase genotypes, and improve filtering accuracy.

**Workflow**

1. **Joint Genotyping**  
   Run `HaplotypeCaller` in **GVCF** mode for all family members → `GenotypeGVCFs` → single multi-sample VCF.

2. **Relationship Verification**  
   Use **KING**, **PLINK**, or **Somalier** to verify reported relationships and detect swaps.

3. **De Novo Detection**  
   Tools:
   - **GATK VariantAnnotator** with `-A PossibleDeNovo`
   - **DeNovoGear**
   - **TrioDeNovo**

   Example (GATK de novo annotation):

   .. code-block:: bash

      gatk VariantAnnotator \
        -R refs/GRCh38.fa \
        -V family.joint.vcf.gz \
        -A PossibleDeNovo \
        --pedigree family.ped \
        -O family.denovo_annotated.vcf.gz

4. **Mendelian Inconsistency Checks**  
   Use **bcftools +mendelian** plugin to identify violations.

5. **Filtering and Prioritization**  
   - Restrict to high-confidence de novo calls (PASS, depth >10 in all samples).
   - Annotate with **VEP** or **ANNOVAR** to assess functional impact and allele frequency.

**Tips**
- Include only variants in high-confidence regions (e.g., GIAB).
- For exomes, restrict analysis to targets ± padding.

Joint Calling for Large Cohorts
-------------------------------

**Advantages**
- Consistent genotyping across samples.
- Improved variant quality metrics (VQSR modeling).
- Better sensitivity for rare variant discovery.

**Best Practices**
- Always call variants in **GVCF** mode per sample.
- Use `GenomicsDBImport` or `CombineGVCFs` for aggregation.
- For >500 samples, use **scatter/gather** and **sharded GenomicsDB workspaces**.

Example:

.. code-block:: bash

   gatk GenomicsDBImport \
     --genomicsdb-workspace-path cohort_gdb \
     --sample-name-map gvcf_samples.map \
     --reader-threads 8 \
     --batch-size 50

   gatk GenotypeGVCFs \
     -R refs/GRCh38.fa \
     -V gendb://cohort_gdb \
     -O cohort.unfiltered.vcf.gz

.. note::
   VQSR requires several thousand high-quality SNPs/indels for stable model training. For small cohorts, use hard filters.

Copy Number Variation (CNV) Analysis
------------------------------------

**Germline CNVs**
- Tools: **GATK gCNV**, **XHMM**, **ExomeDepth** (WES-specific).
- Requires matched control cohort or PoN.

Example (GATK gCNV WES mode):

.. code-block:: bash

   gatk PreprocessIntervals \
     -R refs/GRCh38.fa \
     --interval-merging-rule OVERLAPPING_ONLY \
     -L refs/targets_exome.bed \
     -O targets.interval_list

   gatk AnnotateIntervals \
     -R refs/GRCh38.fa \
     -L targets.interval_list \
     -O targets.annotated.tsv

   gatk CollectReadCounts \
     -I SAMPLE.recal.bam \
     -L targets.interval_list \
     --interval-merging-rule OVERLAPPING_ONLY \
     -O SAMPLE.counts.hdf5

**Somatic CNVs**
- Tool: **GATK ModelSegments** with matched tumor-normal or tumor-only modes.
- Input: read counts + allele counts from common SNPs.

Example (somatic CNV):

.. code-block:: bash

   gatk CollectReadCounts \
     -I tumor.bam -L genome_intervals.interval_list \
     --interval-merging-rule OVERLAPPING_ONLY \
     -O tumor.counts.hdf5

   gatk CollectAllelicCounts \
     -I tumor.bam -R refs/GRCh38.fa \
     -L snps.interval_list \
     -O tumor.allelicCounts.tsv

   gatk ModelSegments \
     --denoised-copy-ratios tumor.denoisedCR.tsv \
     --allelic-counts tumor.allelicCounts.tsv \
     -O cnv_segments/

   gatk CallCopyRatioSegments \
     -I cnv_segments/called_copy_ratios.seg \
     -O cnv_segments/cnv_calls.seg

Mitochondrial Variant Analysis
------------------------------

Mitochondrial variants can be called from WGS (and sometimes WES) using GATK Mutect2 in mitochondrial mode.

Example:

.. code-block:: bash

   gatk Mutect2 \
     -R refs/GRCh38.fa \
     -I SAMPLE.recal.bam \
     --mitochondria-mode \
     -O SAMPLE.mitochondria.unfiltered.vcf.gz

   gatk FilterMutectCalls \
     -V SAMPLE.mitochondria.unfiltered.vcf.gz \
     -O SAMPLE.mitochondria.filtered.vcf.gz

.. note::
   For mitochondrial analysis, use a dedicated circularized mtDNA reference to minimize alignment artifacts.

Tumor Purity Estimation & Adjustment
------------------------------------

**Why it matters**
- Low purity can reduce somatic variant allele fractions (VAFs), impacting sensitivity.

**Tools**
- **PureCN** (integrates CNV + VAF).
- **FACETS**, **ABSOLUTE**.
- **Sequenza** (tumor/normal only).

Example (FACETS tumor-normal):

.. code-block:: bash

   Rscript run_facets.R \
     --tumor tumor.bam \
     --normal normal.bam \
     --genome hg38 \
     --output facets_results/

**Adjustment Strategies**
- Lower `--min-af` in Mutect2 for low-purity tumors.
- Use CNV-aware filtering to retain real variants with low VAF in amplified regions.
- Flag likely subclonal variants for separate interpretation.

Integration and Reporting
-------------------------

- Combine SNVs/indels, CNVs, and mtDNA calls into a unified **multi-omic report**.
- Summarize per-sample key metrics: purity, ploidy, mutation burden, % genome altered.
- Annotate variants with clinical relevance (COSMIC, ClinVar, OncoKB) — **interpretation only**.

Next Steps
----------

- See :doc:`structure_and_containerisation` for reproducibility and execution setup.
- Refer to :doc:`references` for literature on family-based designs, CNV/mtDNA detection, and tumor purity modeling.
