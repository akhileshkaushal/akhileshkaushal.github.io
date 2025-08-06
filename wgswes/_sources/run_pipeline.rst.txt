Run the WGS/WES Pipeline
========================

This section describes an end-to-end, **config-driven** workflow for WGS/WES data. The same pipeline supports **germline** and **somatic** analyses and can be switched via a single YAML configuration file.

.. contents::
   :local:
   :depth: 2

Pipeline Inputs & Layout
-------------------------

Recommended project layout:

.. code-block:: text

    project/
    ├── config/
    │   └── settings.yaml
    ├── refs/
    │   ├── GRCh38.fa
    │   ├── GRCh38.fa.fai
    │   ├── GRCh38.dict
    │   ├── dbsnp.vcf.gz
    │   ├── dbsnp.vcf.gz.tbi
    │   ├── Mills_and_1000G_gold_standard.indels.vcf.gz
    │   ├── Mills_and_1000G_gold_standard.indels.vcf.gz.tbi
    │   ├── 1000G_phase1.indels.vcf.gz
    │   ├── 1000G_phase1.indels.vcf.gz.tbi
    │   ├── af-only-gnomad.vcf.gz
    │   ├── af-only-gnomad.vcf.gz.tbi
    │   └── targets_exome.bed  # WES only (+ optional padded BED)
    ├── fastq/
    │   ├── SAMPLE1_R1.fastq.gz
    │   ├── SAMPLE1_R2.fastq.gz
    │   └── ...
    ├── work/      # intermediate outputs
    └── results/   # final QC, BAM/CRAM, VCF, annotation

Example **settings.yaml** (edit to your environment):

.. code-block:: yaml

    mode: germline  # one of: germline | somatic_tn | somatic_tumor_only
    reference: refs/GRCh38.fa
    known_sites:
      dbsnp: refs/dbsnp.vcf.gz
      mills: refs/Mills_and_1000G_gold_standard.indels.vcf.gz
      indels: refs/1000G_phase1.indels.vcf.gz
    population_af: refs/af-only-gnomad.vcf.gz
    wes:
      targets_bed: refs/targets_exome.bed
      padded_bed: refs/targets_exome.padded100.bed  # optional
    resources:
      threads: 16
      memory_gb: 64
    containers:
      use_apptainer: true
      image: docker://broadinstitute/gatk:latest
    bams:
      compress_to_cram: true
      cram_ref: refs/GRCh38.fa
    scatter:
      intervals: 22  # chromosomes or interval shards for parallelization

.. note::

    For **WES**, use the vendor-provided BED; consider generating a padded BED (+/- 50–100 bp) to capture near-target indels.

Step-by-Step Workflow
---------------------

1) Raw Read QC
~~~~~~~~~~~~~~

- Run **FastQC** on all FASTQ pairs; aggregate with **MultiQC**.

.. code-block:: bash

    mkdir -p results/qc/fastqc
    fastqc fastq/*_R1.fastq.gz fastq/*_R2.fastq.gz -o results/qc/fastqc
    multiqc results/qc -o results/qc

**Gate**: No catastrophic adapter content; base qualities reasonable; overrepresented sequences understood (e.g., adapters).

2) Adapter/Quality Trimming (optional but recommended for WES)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Use **fastp** (recommended) or **Trimmomatic**. Enable poly-G trimming for NovaSeq data if needed.

.. code-block:: bash

    fastp \
      -i fastq/SAMPLE_R1.fastq.gz -I fastq/SAMPLE_R2.fastq.gz \
      -o work/SAMPLE_R1.trim.fq.gz -O work/SAMPLE_R2.trim.fq.gz \
      --detect_adapter_for_pe --thread 16 \
      --html results/qc/SAMPLE.fastp.html --json results/qc/SAMPLE.fastp.json

**Outputs**: Trimmed FASTQs + fastp QC reports (included in MultiQC).

3) Alignment to Reference
~~~~~~~~~~~~~~~~~~~~~~~~~

- Align with **BWA-MEM2** (or BWA-MEM) against GRCh38/hg38 (ALT-aware). Always set **Read Group (RG)** tags.

.. code-block:: bash

    RG='@RG\tID:SAMPLE\tSM:SAMPLE\tPL:ILLUMINA\tLB:LIB1\tPU:FLOWCELL.LANE'
    bwa-mem2 mem -t 16 refs/GRCh38.fa \
      work/SAMPLE_R1.trim.fq.gz work/SAMPLE_R2.trim.fq.gz \
      | samtools view -b - \
      | samtools sort -@ 8 -o work/SAMPLE.sorted.bam
    samtools index work/SAMPLE.sorted.bam

**QC**:

- Alignment stats: ``samtools flagstat`` and (WES) **Picard CollectHsMetrics** / (WGS) **Picard CollectWgsMetrics**.
- Coverage profiling: **mosdepth** or **bedtools genomecov**.

.. code-block:: bash

    samtools flagstat work/SAMPLE.sorted.bam > results/qc/SAMPLE.flagstat.txt
    # WES capture metrics
    gatk CollectHsMetrics \
      -I work/SAMPLE.sorted.bam -O results/qc/SAMPLE.hs_metrics.txt \
      -R refs/GRCh38.fa --BAIT_INTERVALS refs/targets_exome.bed \
      --TARGET_INTERVALS refs/targets_exome.bed
    # WGS metrics
    gatk CollectWgsMetrics \
      -I work/SAMPLE.sorted.bam -O results/qc/SAMPLE.wgs_metrics.txt -R refs/GRCh38.fa

4) Duplicate Marking & BQSR
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Mark duplicates (UMI-aware if applicable). Perform **Base Quality Score Recalibration** (BQSR) using known sites.

.. code-block:: bash

    gatk MarkDuplicatesSpark \
      -I work/SAMPLE.sorted.bam \
      -O work/SAMPLE.dedup.bam \
      -M results/qc/SAMPLE.dup_metrics.txt
    samtools index work/SAMPLE.dedup.bam
    gatk BaseRecalibrator \
      -I work/SAMPLE.dedup.bam -R refs/GRCh38.fa \
      --known-sites refs/dbsnp.vcf.gz \
      --known-sites refs/Mills_and_1000G_gold_standard.indels.vcf.gz \
      --known-sites refs/1000G_phase1.indels.vcf.gz \
      -O work/SAMPLE.recal_data.table
    gatk ApplyBQSR \
      -I work/SAMPLE.dedup.bam -R refs/GRCh38.fa \
      --bqsr-recal-file work/SAMPLE.recal_data.table \
      -O work/SAMPLE.recal.bam
    samtools index work/SAMPLE.recal.bam

.. tip::

    To reduce storage, write **CRAM** instead of BAM:
    ``gatk PrintReads -I work/SAMPLE.recal.bam -O results/bam/SAMPLE.cram -R refs/GRCh38.fa``

5) Variant Calling
~~~~~~~~~~~~~~~~~~

**Germline (per-sample GVCF)**

- Call each sample in **GVCF** mode for joint genotyping later.

.. code-block:: bash

    gatk HaplotypeCaller \
      -R refs/GRCh38.fa -I work/SAMPLE.recal.bam \
      -O results/gvcf/SAMPLE.g.vcf.gz -ERC GVCF

**Joint genotyping (cohort)**

.. code-block:: bash

    gatk GenomicsDBImport \
      --genomicsdb-workspace-path work/cohort_gdb \
      --sample-name-map config/gvcf_samples.map \
      --reader-threads 8
    gatk GenotypeGVCFs \
      -R refs/GRCh38.fa \
      -V gendb://work/cohort_gdb \
      -O results/vcf/cohort.unfiltered.vcf.gz

.. note::

    For small cohorts, you may replace GenomicsDB with ``CombineGVCFs`` but GenomicsDB scales better.

**Somatic (Tumor–Normal)**

- Use **Mutect2** with matched normal, a **germline resource** (gnomAD), and a **Panel of Normals (PoN)** if available.

.. code-block:: bash

    # Calling
    gatk Mutect2 \
      -R refs/GRCh38.fa \
      -I work/TUMOR.recal.bam -tumor TUMOR \
      -I work/NORMAL.recal.bam -normal NORMAL \
      --germline-resource refs/af-only-gnomad.vcf.gz \
      --panel-of-normals refs/pon.vcf.gz \
      -O work/TUMOR.mutect2.unfiltered.vcf.gz
    # Orientation bias model (important for FFPE)
    gatk LearnReadOrientationModel \
      -I work/TUMOR.f1r2.tar.gz \
      -O work/TUMOR.read-orientation-model.tar.gz
    # Contamination estimation
    gatk GetPileupSummaries \
      -I work/TUMOR.recal.bam -V refs/af-only-gnomad.vcf.gz \
      -L refs/targets_exome.bed -O work/TUMOR.pileups.table
    gatk GetPileupSummaries \
      -I work/NORMAL.recal.bam -V refs/af-only-gnomad.vcf.gz \
      -L refs/targets_exome.bed -O work/NORMAL.pileups.table
    gatk CalculateContamination \
      -I work/TUMOR.pileups.table \
      -matched work/NORMAL.pileups.table \
      -O work/TUMOR.contamination.table \
      --tumor-segmentation work/TUMOR.segments.table
    # Filtering
    gatk FilterMutectCalls \
      -R refs/GRCh38.fa \
      -V work/TUMOR.mutect2.unfiltered.vcf.gz \
      --contamination-table work/TUMOR.contamination.table \
      --ob-priors work/TUMOR.read-orientation-model.tar.gz \
      -O results/vcf/TUMOR.somatic.filtered.vcf.gz

**Somatic (Tumor-only)**

- Call without a matched normal. Use **PoN** and gnomAD to suppress artifacts and likely germline.

.. code-block:: bash

    gatk Mutect2 \
      -R refs/GRCh38.fa \
      -I work/TUMOR.recal.bam -tumor TUMOR \
      --germline-resource refs/af-only-gnomad.vcf.gz \
      --panel-of-normals refs/pon.vcf.gz \
      -O work/TUMOR.to.unfiltered.vcf.gz
    gatk LearnReadOrientationModel \
      -I work/TUMOR.f1r2.tar.gz \
      -O work/TUMOR.read-orientation-model.tar.gz
    gatk GetPileupSummaries \
      -I work/TUMOR.recal.bam -V refs/af-only-gnomad.vcf.gz \
      -O work/TUMOR.pileups.table
    gatk CalculateContamination \
      -I work/TUMOR.pileups.table -O work/TUMOR.contamination.table
    gatk FilterMutectCalls \
      -R refs/GRCh38.fa \
      -V work/TUMOR.to.unfiltered.vcf.gz \
      --contamination-table work/TUMOR.contamination.table \
      --ob-priors work/TUMOR.read-orientation-model.tar.gz \
      -O results/vcf/TUMOR.somatic.filtered.vcf.gz

.. tip::

    **Building a PoN**: Run Mutect2 on ≥30 normal samples (as “tumors” without normal), then
    ``gatk CreateSomaticPanelOfNormals -V normals.list -O refs/pon.vcf.gz``

6) Variant Filtering
~~~~~~~~~~~~~~~~~~~~

**Germline**

- Prefer **VQSR** for sufficiently large cohorts. Otherwise use **hard filters**.

**VQSR (cohort)**

.. code-block:: bash

    # SNP model
    gatk VariantRecalibrator \
      -R refs/GRCh38.fa -V results/vcf/cohort.unfiltered.vcf.gz \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 refs/hapmap.vcf.gz \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 refs/1000G_omni.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 refs/1000G_high_conf.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 refs/dbsnp.vcf.gz \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      -mode SNP -O work/cohort.snp.recal --tranches-file work/cohort.snp.tranches
    gatk ApplyVQSR \
      -V results/vcf/cohort.unfiltered.vcf.gz \
      --recal-file work/cohort.snp.recal \
      --tranches-file work/cohort.snp.tranches \
      -mode SNP -O work/cohort.snp.recalibrated.vcf.gz
    # INDEL model
    gatk VariantRecalibrator \
      -R refs/GRCh38.fa -V work/cohort.snp.recalibrated.vcf.gz \
      --resource:mills,known=false,training=true,truth=true,prior=12.0 refs/Mills_and_1000G_gold_standard.indels.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 refs/dbsnp.vcf.gz \
      -an QD -an ReadPosRankSum -an FS -an SOR \
      -mode INDEL -O work/cohort.indel.recal --tranches-file work/cohort.indel.tranches
    gatk ApplyVQSR \
      -V work/cohort.snp.recalibrated.vcf.gz \
      --recal-file work/cohort.indel.recal \
      --tranches-file work/cohort.indel.tranches \
      -mode INDEL -O results/vcf/cohort.filtered.vcf.gz

**Hard filters (small cohorts)**

.. code-block:: bash

    # SNPs (example thresholds)
    gatk VariantFiltration \
      -V results/vcf/cohort.unfiltered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "FS > 60.0" --filter-name "FS60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ40" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRS-12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS-8" \
      --filter-expression "SOR > 3.0" --filter-name "SOR3" \
      -O work/cohort.snps.hardfiltered.vcf.gz
    # INDELs (example thresholds)
    gatk VariantFiltration \
      -V work/cohort.snps.hardfiltered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "FS > 200.0" --filter-name "FS200" \
      --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS-20" \
      --filter-expression "SOR > 10.0" --filter-name "SOR10" \
      -O results/vcf/cohort.filtered.vcf.gz

**Somatic**

- Use **FilterMutectCalls** output; optionally apply **FilterByOrientationBias** for specific FFPE artifacts:

.. code-block:: bash

    gatk FilterByOrientationBias \
      -V results/vcf/TUMOR.somatic.filtered.vcf.gz \
      -P G/T -P C/T \
      -O results/vcf/TUMOR.somatic.ffpe_filtered.vcf.gz

7) Annotation
~~~~~~~~~~~~~

Annotate for consequence, population frequency, and clinical context.

**VEP (recommended)**

.. code-block:: bash

    vep \
      -i results/vcf/*.vcf.gz \
      -o results/annotation/annotated.vep.vcf.gz \
      --vcf --compress_output bgzip \
      --assembly GRCh38 --offline --cache \
      --fork 8 --everything

**ANNOVAR (alternative)**

.. code-block:: bash

    table_annovar.pl results/vcf/*.vcf.gz humandb/ \
      -buildver hg38 -out results/annotation/annovar \
      -remove -protocol refGene,gnomad211_exome,clinvar_202403 \
      -operation g,f,f -nastring . -polish -vcfinput

.. warning::

    **ClinVar** and cancer knowledge bases are for **interpretation**. Do not use them as hard filters without a defined policy.

8) Cohort-Level QC & Reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Aggregate QC with **MultiQC** (FastQC, alignment metrics, coverage, duplication).
- Confirm sample identity/sex and detect swaps with **Somalier** or **VerifyBamID2**.
- Summarize callable bases and % bases ≥20×/30× (WES) or ≥30× (WGS).

Parallelization & Performance
-----------------------------

- **Scatter/gather** by chromosomes or interval shards; keep read group boundaries intact.
- Use **MarkDuplicatesSpark** for speed (ensure adequate I/O).
- Tune BWA-MEM2 threads to CPU/NUMA; avoid oversubscription.

Example Slurm snippet:

.. code-block:: bash

    #!/bin/bash
    #SBATCH -J wgswes
    #SBATCH -N 1
    #SBATCH -c 16
    #SBATCH --mem=64G
    #SBATCH -t 24:00:00
    module load apptainer
    apptainer exec docker://broadinstitute/gatk:latest \
      gatk --version
    # run HaplotypeCaller / Mutect2 shards here ...

WES-Specific Notes
------------------

- Report **on-target %**, **fold-80 base penalty**, **mean target coverage**, and **% targets ≥ 20×/30×**.
- Use padded BED for calling; revert to strict BED for coverage/QC.
- Capture-kit updates may change target definitions; pin versions in the repository.

Expected Outputs
----------------

- **BAM/CRAM** + index (post-BQSR), duplication metrics, alignment/coverage metrics.
- **VCF/BCF**: Germline joint genotyped and filtered; or somatic filtered calls.
- **Annotation** files (VEP/ANNOVAR outputs).
- **QC**: MultiQC report, per-sample metrics, contamination estimates.

Troubleshooting
---------------

- **Low on-target (WES)**: Verify BED compatibility, library quality, and adapter trimming.
- **High contamination**: Re-evaluate sample swaps; remove contaminated samples from PoN.
- **Over-filtering (germline)**: Revisit hard-filter thresholds or use VQSR with adequate cohort size.
- **Somatic false positives (tumor-only)**: Ensure PoN is technology-matched; tighten population AF thresholds; apply orientation bias filtering.

Next Steps
----------

- See :doc:`advanced_analysis` for trios, joint calling strategies, CNV/mtDNA modules, and tumor purity considerations.
- See :doc:`structure_and_containerisation` for exact containers, pinned versions, and YAML-driven execution wrappers.
- Consult :doc:`references` for methods papers and best-practice guidance.