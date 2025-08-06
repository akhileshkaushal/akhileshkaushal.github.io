WGS/WES Variant Calling Pipeline
================================

:Author: Akhilesh Kaushal
:Version: 1.0.0
:Keywords: WGS, WES, Germline, Somatic, GATK, Apptainer, Docker, Singularity

Overview
--------

This documentation introduces a **comprehensive, reproducible** pipeline for variant discovery from **Whole Genome Sequencing (WGS)** and **Whole Exome Sequencing (WES)** data. It supports both **germline** and **somatic** analyses and can operate **with or without matched normal** samples. The workflow aligns with **GATK Best Practices**, emphasizing high-quality read preprocessing, robust variant discovery, calibrated filtering, and rigorous quality control (QC), followed by standards-compliant annotation for downstream interpretation.

Use Cases
---------

- **Germline**: population-scale studies, rare disease genetics, trio/quad analyses, and cohort joint genotyping.
- **Somatic (Tumor/Normal)**: cancer genomics using matched tumor–normal pairs; optimized for sensitivity with artifact suppression.
- **Somatic (Tumor-only)**: analyses without normal tissue using **Panel of Normals (PoN)** and population allele frequency filtering.

Supported Data Types
--------------------

- **WGS**: 30× coverage typical (higher for somatic); uniform genome-wide capture.
- **WES**: capture-based enrichment (e.g., Agilent SureSelect, IDT xGen). Requires target BED (and optional “padding” for near-target indels).

When to Choose WGS vs WES
-------------------------

- **WGS**: comprehensive coverage (coding + noncoding), superior indel/SV detection, better uniformity; higher cost and storage.
- **WES**: cost-effective for coding regions, higher effective coverage over exons, limited noncoding insight and capture biases.

Design Principles
-----------------

- **Standards-aligned**: adheres to GATK4 Best Practices; no indel realignment (deprecated); **BQSR** performed using known sites.
- **Reproducible**: containerized with **Apptainer/Singularity** or **Docker**; deterministic configuration captured in environment files.
- **Modular**: pluggable steps for alignment, calling, filtering, and annotation; optional mitochondrial and CNV modules.
- **Scalable**: HPC/Cloud-friendly; parallelization by sample/chromosome/scatter-gather.
- **Auditable**: rich QC (FastQC/MultiQC, alignment metrics, coverage profiling) and provenance tracking.

High-Level Workflow
-------------------

**Common Preprocessing**
 - Input FASTQ (gz)
 - Adapter/quality assessment (**FastQC**, **MultiQC**)
 - Alignment (**BWA-MEM2** or **BWA-MEM**) to **GRCh38/hg38** (ALT-aware, decoys if applicable)
 - Sorting, duplicate marking (**Picard/GATK MarkDuplicates(Spark)**; UMI-aware dedup possible)
 - Base Quality Score Recalibration (**BQSR**) using known variant sites (dbSNP, Mills indels, 1000G indels)

**Germline Path**
 - Per-sample calling: **GATK HaplotypeCaller** in **GVCF** mode
 - Cohort integration: **GenomicsDBImport** (or CombineGVCFs) → **GenotypeGVCFs**
 - Variant filtering: **VQSR** (recommended for sufficiently large cohorts) or **hard filters** (small cohorts)
 - Annotation: **VEP/Funcotator/ANNOVAR** with gnomAD, dbSNP, dbNSFP, ClinVar (interpretation only), ±HGVS

**Somatic (Tumor/Normal) Path**
 - Calling: **GATK Mutect2** with matched normal, **germline resource** (e.g., gnomAD), and optional **PoN**
 - Post-processing: **LearnReadOrientationModel**, **CalculateContamination**, **FilterMutectCalls**, optional **FilterByOrientationBias** (FFPE)
 - Optional modules: mitochondrial calling (Mutect2 mitochondria mode), CNV (GATK ModelSegments/CallCopyRatioSegments)

**Somatic (Tumor-only) Path**
 - Calling: **Mutect2** without matched normal
 - Artifact suppression: **PoN** (ideally ≥30 normals, technology-matched), gnomAD-based germline filtering, orientation bias model, contamination estimation
 - **Caveat**: higher false-positive risk and limited certainty of somatic status; careful post-filters and orthogonal validation recommended.

Inputs and Outputs
------------------

**Inputs**
 - Paired-end **FASTQ** files (R1/R2), consistent sample naming
 - **Reference** genome (GRCh38/hg38 recommended), indices and dictionary
 - Known sites for BQSR (dbSNP; Mills/1000G indels)
 - **WES**: target BED (+ optional padded BED)
 - **Somatic**: matched normal BAM/FASTQ when available; **PoN** VCF (or cohort of normals to build one); gnomAD resource

**Outputs**
 - Recalibrated **BAM/CRAM** + index and metrics
 - **VCF/BCF** (germline joint-called, or somatic filtered calls), with indices
 - QC reports (MultiQC), coverage summaries (e.g., **mosdepth**, **Qualimap**)
 - Annotated variant files (e.g., VEP/Funcotator outputs) for interpretation

Quality, Coverage, and QC
-------------------------

- **WGS**: ≥30× mean depth for germline; somatic studies often require higher tumor depth (e.g., 60–100×) and ≥30× normal.
- **WES**: ≥80–120× mean on-target; report on-target %, uniformity, and % bases ≥20×/30×.
- **Contamination**: estimate via **VerifyBamID2** (germline) or **CalculateContamination** (somatic); flag swaps with **Somalier**.
- **Sex chromosomes**: handle PAR/non-PAR; report sex inference for consistency checks.
- **UMIs**: if present, use UMI-aware collapsing/dedup (e.g., **fgbio**) to reduce PCR artifacts.

Filtering Strategy Notes
------------------------

- **VQSR** is preferred for large cohorts (stable tranche modeling); for small studies, use **hard filters** with empirically chosen cutoffs.
- **Somatic** filtering combines: artifact modeling (orientation/FFPE), contamination estimates, PoN, and population AF thresholds.
- **ClinVar** is for **interpretation**, not for hard filtering criteria.

Reference and Resource Recommendations
--------------------------------------

- **Genome**: GRCh38/hg38 with ALT contigs and decoys where applicable.
- **Known sites**: dbSNP (latest), Mills and 1000G indels.
- **Population AF**: **gnomAD** (exomes + genomes) as a germline resource in somatic workflows.
- **Cancer knowledge bases** (annotation only): **COSMIC**, **CIViC**, **OncoKB** (license terms may apply).

Reproducibility & Execution
---------------------------

- Provided **Apptainer/Singularity** and **Docker** images for stable toolchains.
- YAML/JSON configuration captures references, intervals, parameters, and resource sizing.
- Scales on Slurm/SGE/PBS or cloud batch services; supports scatter/gather by interval list.

Data Stewardship & Compliance
-----------------------------

- Maintain sample sheets with immutable IDs, checksums for all inputs/outputs, and complete metadata.
- Remove or avoid embedding PHI in filenames/headers; adhere to **HIPAA/GDPR** and institutional IRB policies.
- Store large artifacts in object storage with lifecycle policies; keep a minimal, query-friendly variant store.

Limitations & Non-Goals
-----------------------

- Tumor-only analyses cannot definitively distinguish all germline from somatic variation; interpret with caution.
- Structural variant and repeat expansion calling are **optional** modules and may require specialized callers and validation.
- Clinical reporting requires orthogonal confirmation and domain-specific review beyond this pipeline.

What’s Next
-----------

- See :doc:`run_pipeline` for exact commands, parameters, and file layouts.
- Explore :doc:`advanced_analysis` for joint calling, trios, tumor purity handling, CNV/mtDNA options, and best-practice filters.
- Review :doc:`structure_and_containerisation` for environment reproducibility.
- Consult :doc:`references` for standards and key literature.

.. note::
   This pipeline assumes high-quality libraries and appropriate experimental design. For FFPE, ultra-low input, or single-cell protocols, additional artifact controls and validation are recommended.

