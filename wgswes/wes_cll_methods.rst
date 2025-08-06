Stage-to-Tool Mapping With Example Commands
-------------------------------------------

This section provides a comprehensive mapping of bioinformatics pipeline stages to commonly used tools, drawing from the workflows and methodologies highlighted in the two case studies: Yang et al. (2024) on pharmacogenomic profiling of intra-tumor heterogeneity in liver cancer, and Ljungström et al. (2016) on whole-exome sequencing in relapsing chronic lymphocytic leukemia (CLL). The tools are selected based on their relevance to whole-exome sequencing (WES), whole-genome sequencing (WGS), variant calling, copy number analysis, clonal evolution modeling, transcriptomic integration, and mutation signature analysis.

The table below outlines key stages, recommended tools (with latest versions as of August 2025, verified via web searches), and detailed notes on their application. Tools from the case studies are explicitly referenced, with modern equivalents or updates provided for reproducibility and performance. For instance, Yang et al. emphasized bwa-mem2 (v2.0 at the time, now v2.2.1) for alignment, GATK (v4.1.2.0, now v4.5.0.0) for variant calling, CNVkit (v0.9.7.b1, now v0.9.10) for CNAs, and Salmon/edgeR for RNA-seq. Ljungström et al. used BWA (v0.7.5, now v0.7.18), GATK (v2.8, now v4.5.0.0), ANNOVAR (2014-07-14, now 2024May25), FACETS (v0.3.9, now v0.6.2), SciClone (v1.0, now v1.1.0), DeconstructSigs (v1.6.0, now v1.9.0), MutSigCV (v1.41), and dNdScv (v0.1.0).

Where possible, tools are updated to their latest versions (e.g., GATK 4.5.0.0, ANNOVAR 2024May25) to leverage improvements in speed, accuracy, and compatibility. For ITH-aware and clonal evolution analyses, tools like PyClone-VI (v0.1.6), SciClone (v1.1.0), and FACETS (v0.6.2) are highlighted, with adaptations for multi-region or paired relapse samples. RNA-seq tools (e.g., STAR v2.7.11b, Salmon v1.10.3, edgeR v3.46.1, DESeq2 v1.44.0) are included for transcriptomic integration, as in Yang et al. Mutation signature tools (e.g., SigProfiler v1.1.25, DeconstructSigs v1.9.0, MutSigCV v1.41, dNdScv v0.1.0) are added for downstream analyses.

The example commands use Bash scripts with modern versions, assuming a Linux/Unix environment (e.g., Ubuntu 24.04). Commands are verbose, with explanations inline. Use containers (Docker/Singularity) for reproducibility, as in the case studies (e.g., GATK Docker images). All paths are placeholders; replace with actual values.

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Stage
     - Common Tools (Latest Versions)
     - Notes
   * - Inputs & Metadata
     - Sample sheet (CSV/TSV), PoN VCF (e.g., from GATK v4.5.0.0), gnomAD resource (v4.1), capture BED (e.g., Agilent SureSelect v6)
     - Prepare metadata with sample IDs, tumor-normal pairs, timepoints (pre/post-relapse in CLL), regions (multi-region in liver cancer). Include genome build (GRCh38/hg38 preferred). Verify FASTQ naming for paired-end reads. From studies: Yang used multi-regional samples; Ljungström used paired pre/post-relapse. Use CSV for compatibility with tools like GATK and DESeq2.
   * - QC of Raw Reads
     - FastQC (v0.12.1), MultiQC (v1.30)
     - Assess per-base quality, adapters, overrepresented sequences, duplication. MultiQC aggregates across samples. In studies: Yang and Ljungström emphasized QC before alignment. For RNA-seq, check GC bias. Use --threads for speed.
   * - Trimming (Optional)
     - fastp (v1.0.1), Cutadapt (v5.1)
     - Trim adapters/low-quality tails. fastp is all-in-one (QC+trim+filter). From Yang: Used for RNA-seq preprocessing. Skip if data is clean; record params for ITH consistency.
   * - Alignment
     - bwa-mem2 (v2.2.1), STAR (v2.7.11b), Salmon (v1.10.3) for quasi-mapping
     - Align to GRCh38. bwa-mem2 for WES/WGS (Yang/Ljungström). STAR for RNA-seq (Yang). Salmon for transcript quasi-mapping. Handle multi-region/relapse by aligning per sample. Mark duplicates post-alignment.
   * - Post-Alignment Processing
     - GATK (v4.5.0.0), Picard (v3.2.0), samtools (v1.22.1)
     - Mark duplicates, BQSR, index BAMs. GATK for recalibration (both studies). samtools for sorting/indexing. Essential for variant calling accuracy in heterogeneous samples.
   * - Somatic SNV/Indel Calling
     - GATK Mutect2 (v4.5.0.0), Strelka2 (v2.9.10)
     - Tumor-normal paired mode (prioritized in both studies). Filter with PoN/gnomAD. Yang used MuTect v2.0; Ljungström v1.1.4. Adjust VAF thresholds for subclones.
   * - CNAs
     - CNVkit (v0.9.10), FACETS (v0.6.2), GISTIC (v2.0.23)
     - CNVkit for segmentation (Yang). FACETS for purity/ploidy (Ljungström). GISTIC for focal events. ITH-aware: Quantify CNA-ITH (Yang). Use matched normals.
   * - Clonality & Evolution
     - PyClone-VI (v0.1.6), SciClone (v1.1.0), FACETS (v0.6.2)
     - SciClone for subclonal clustering (Ljungström). PyClone-VI for ITH (modern alternative). Integrate VAF/CN for CCF. Phylogenetic trees with MEGA (v12.0.13, both studies).
   * - Annotation
     - ANNOVAR (2024May25), Ensembl VEP (v110)
     - Functional impact, COSMIC/ClinVar. ANNOVAR in both studies (2019-10/2014-07). VEP for modern updates. Handle multi-allelic variants.
   * - Mutation Signatures
     - SigProfiler (v1.1.25), DeconstructSigs (v1.9.0), MutSigCV (v1.41), dNdScv (v0.1.0)
     - DeconstructSigs/MutSigCV/dNdScv for signatures (Ljungström). SigProfiler for de novo extraction. Analyze per-region/timepoint for evolution.
   * - Transcriptomics (RNA-seq)
     - STAR (v2.7.11b), Salmon (v1.10.3), edgeR (v3.46.1), DESeq2 (v1.44.0)
     - STAR for alignment, Salmon for quasi-mapping (Yang). edgeR/DESeq2 for DE/signatures. Integrate with WES for multi-omics ITH.
   * - Reporting & Visualization
     - MultiQC (v1.30), R/ggplot2 (v3.5.1)
     - Aggregate QC/reports. Custom plots for phylogenies (MEGA), signatures. Studies used R for visualizations.

A key aspect of the pipeline, particularly relevant to the case studies, is handling intra-tumor heterogeneity (ITH) and clonal evolution. This builds directly on the Clonality & Evolution stage, providing specialized adaptations for multi-region (liver cancer) and paired relapse (CLL) analyses to classify trunk/branch variants, compute cancer cell fractions (CCF), and reconstruct phylogenies.

ITH/Evolution Adaptations: Intra-tumor heterogeneity (ITH) and clonal evolution are critical in cancers like primary liver cancer (PLC) and chronic lymphocytic leukemia (CLL), as highlighted in the case studies. In Yang et al. (2024), multi-region sampling revealed branch-dominant ITH (region-specific mutations/CNAs) linked to drug resistance, with trunk mutations shared across regions (early events) and branch mutations private (late events). Trunk/branch classification uses phylogenetic reconstruction to infer evolutionary trees. In Ljungström et al. (2016), paired pre/post-relapse samples showed subclonal expansion (e.g., TP53 clones), quantified via cancer cell fraction (CCF) to track dynamics.

For multi-region analysis (e.g., liver cancer from Yang et al.), first consolidate VCFs from individual regions using bcftools to create a multi-sample VCF, then classify variants as trunk (shared across all regions, CCF >0.9 in all), branch (shared in subsets), or private (region-specific). Use phylogenetic tools like MEGA (v12.0.13, command-line MEGA-CC for automation) for tree building based on binary presence/absence matrices, or PyClone-VI (v0.1.6) for probabilistic clustering integrating VAF and CN. Compute trunk ratio (TR = trunk mutations / branch mutations) in R/Python for ITH quantification.

For paired relapse analysis (e.g., CLL from Ljungström et al.), compute CCF with FACETS (v0.6.2) to adjust VAF for purity/ploidy, enabling subclonal inference (e.g., expansion if CCF increases post-relapse). FACETS outputs segments and fits for downstream tools like SciClone.

**Bash Scripts for Multi-Region (Liver Cancer - Consolidate and Classify)**

Prepare a presence/absence matrix from consolidated VCF (e.g., using bcftools query to extract variants per region), then build trees.

.. code-block:: bash

   # Assume VCFs: region1.vcf.gz, region2.vcf.gz, region3.vcf.gz
   # Step 1: Consolidate VCFs (merge multi-sample)
   bcftools merge -Oz -o multi_region.vcf.gz region1.filtered.vcf.gz region2.filtered.vcf.gz region3.filtered.vcf.gz
   bcftools index multi_region.vcf.gz

   # Step 2: Extract presence/absence (custom Python/R; example with bcftools + awk for binary matrix)
   # Create binary matrix: rows=variants, columns=regions (1 if ALT>0, 0 otherwise)
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' multi_region.vcf.gz | awk 'BEGIN {print "Variant\tRegion1\tRegion2\tRegion3"} {var=$1":"$2":"$3":"$4; r1=($5!="0/0" && $5!="./.")?1:0; r2=($6!="0/0" && $6!="./.")?1:0; r3=($7!="0/0" && $7!="./.")?1:0; print var"\t"r1"\t"r2"\t"r3}' > presence_matrix.tsv

   # Step 3: Phylogenetic tree with MEGA-CC (command-line; assumes MEGA-CC installed)
   # Convert to MEGA format (.meg) via custom script, then:
   megacc -a neighbor_joining.mao -d presence_matrix.meg -o tree.nwk  # Neighbor-joining tree

   # Alternative: PyClone-VI (v0.1.6) for clustering (prepare input with VAF/CN per region)
   # Input: tsv with mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn per region
   pyclone-vi fit --in-file multi_region_input.tsv --out-file pyclone_results.h5 --num-clusters 5 --num-restarts 100 --seed 42
   pyclone-vi write-results-file --in-file pyclone_results.h5 --out-file pyclone_posterior.tsv

   # Step 4: Classify trunk/branch (custom R/Python; example R snippet)
   Rscript -e 'matrix <- read.table("presence_matrix.tsv", header=TRUE); trunk <- rowSums(matrix[,2:4]) == 3; branch <- rowSums(matrix[,2:4]) > 1 & rowSums(matrix[,2:4]) < 3; private <- rowSums(matrix[,2:4]) == 1; cat("TR =", sum(trunk)/sum(branch), "\n")'

   # Explanation: Merge VCFs, create binary matrix, build tree (MEGA for parsimony), or cluster (PyClone-VI Bayesian). TR quantifies ITH (Yang: TR >1 trunk-dominant).

**Bash Scripts for Paired Relapse (CLL - Compute CCF with FACETS)**

FACETS is R-based; wrap in Rscript or use facets-suite wrapper.

.. code-block:: bash

   # Install FACETS if needed (in R: install.packages("facets"))
   # Run FACETS (v0.6.2) via Rscript
   Rscript -e 'library(facets); snp_pileup("tumor.bam", "normal.bam", "snp_positions.gz", output_file="pileup.dat.gz"); preProcSample("pileup.dat.gz", snp.nbhd=250, ndepth=35, cval=25); procSample(rcmat, cval=150); fit <- emcncf(xx); write.table(fit$cncf, "facets_cncf.tsv", sep="\t", quote=FALSE, row.names=FALSE)'

   # Alternative using facets-suite wrapper (if installed: https://github.com/mskcc/facets-suite)
   run-facets-wrapper.R --tumor-bam tumor.bam --normal-bam normal.bam --output-prefix output --snp-vcf snp_positions.vcf.gz --cval-preproc 25 --cval-proc 150 --ndepth 35

   # Compute CCF (post-FACETS; custom)
   # Use fit$cncf for segments; for variants, adjust VAF: CCF = (VAF * purity) / ( (2 * (1 - purity)) + (purity * total_cn) )
   Rscript -e 'variants <- read.table("variants.tsv", header=TRUE); purity <- 0.8; # from FACETS; ccf <- (variants$VAF * purity) / ((2 * (1 - purity)) + (purity * variants$total_cn)); write.table(ccf, "ccf.tsv")'

   # Explanation: FACETS estimates purity/ploidy/CN segments (Ljungström). Compute CCF for each variant to detect expansion (e.g., if CCF_post > CCF_pre). Integrate with SciClone for clustering. For relapse, compare pre/post CCF to identify selected subclones (e.g., RPS15/TP53).

Example Commands (Bash)
~~~~~~~~~~~~~~~~~~~~~~~

**0) Environment and References**

Set up a Conda environment for reproducibility, incorporating tools from studies.

.. code-block:: bash

   # Create Conda env with key tools
   conda create -n cancer_pipeline -y python=3.11 r-base=4.4.1
   conda activate cancer_pipeline
   conda install -c bioconda -c conda-forge bwa-mem2=2.2.1 gatk4=4.5.0.0 annovar=2024may25 cnvkit=0.9.10 pyclone-vi=0.1.6 facets=0.6.2 sciclone=1.1.0 sigprofiler=1.1.25 deconstructsigs=1.9.0 mutsigcv=1.41 dndscv=0.1.0 star=2.7.11b salmon=1.10.3 bioconductor-edger=3.46.1 bioconductor-deseq2=1.44.0 fastqc=0.12.1 multiqc=1.30 fastp=1.0.1 bcftools=1.22.1 samtools=1.22.1
   conda install -c r r-ggplot2=3.5.1

   # References (GRCh38)
   REF=GRCh38.fa  # From Ensembl/UCSC
   DBSNP=dbsnp_156.vcf.gz  # Latest dbSNP v156
   MILLS=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
   GNOMAD_AF=gnomad.v4.1.vcf.gz  # Latest gnomAD v4.1
   PON=pon.hg38.mutect2.vcf.gz  # Panel of normals
   BAIT_BED=exome_targets_hg38.bed  # e.g., Agilent SureSelect v6

**1) QC and Optional Trimming**

Perform QC with FastQC/MultiQC; trim with fastp if needed.

.. code-block:: bash

   # FastQC and MultiQC
   mkdir qc
   fastqc tumor_R1.fq.gz tumor_R2.fq.gz normal_R1.fq.gz normal_R2.fq.gz -o qc/ --threads 16
   multiqc qc/ -o qc/ --verbose

   # Trimming with fastp (all-in-one QC+trim)
   fastp -i tumor_R1.fq.gz -I tumor_R2.fq.gz -o tumor.trim.R1.fq.gz -O tumor.trim.R2.fq.gz --detect_adapter_for_pe --html qc/tumor_fastp.html --thread 16 --report_title "Tumor QC Report"

   # Explanation: fastp (from Yang) is faster than Cutadapt; use for adapters/low-quality. MultiQC aggregates for multi-sample ITH (Yang/Ljungström).

**2) Alignment and Sorting**

Align with bwa-mem2 (DNA) or STAR (RNA); sort/index with samtools.

.. code-block:: bash

   # DNA alignment (bwa-mem2 from both studies)
   R1=tumor.trim.R1.fq.gz
   R2=tumor.trim.R2.fq.gz
   bwa-mem2 mem -t 16 -R "@RG\tID:TUMOR\tSM:TUMOR\tPL:ILLUMINA" $REF $R1 $R2 | samtools sort -@ 16 -o tumor.sorted.bam -
   samtools index -@ 16 tumor.sorted.bam

   # RNA-seq alignment (STAR from Yang)
   mkdir star_out
   STAR --runThreadN 16 --genomeDir star_index --readFilesIn rna_R1.fq.gz rna_R2.fq.gz --readFilesCommand zcat --outFileNamePrefix star_out/ --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic

   # Explanation: bwa-mem2 for WES/WGS (faster than BWA v0.7.18). STAR for spliced RNA alignment. Handle relapse pairs separately (Ljungström).

**3) Post-Alignment Processing (GATK4)**

Mark duplicates, recalibrate.

.. code-block:: bash

   gatk MarkDuplicates -I tumor.sorted.bam -O tumor.dedup.bam -M tumor.mkdup.txt --CREATE_INDEX true
   gatk BaseRecalibrator -R $REF -I tumor.dedup.bam --known-sites $DBSNP --known-sites $MILLS -O tumor.recal.table
   gatk ApplyBQSR -R $REF -I tumor.dedup.bam --bqsr-recal-file tumor.recal.table -O tumor.bqsr.bam

   # Similar for normal BAM
   # Explanation: GATK (v4.5.0.0 update from studies) for BQSR. Essential for accurate variant calling in heterogeneous tumors.

**4) Somatic SNV/Indel Calling (Mutect2)**

Call variants.

.. code-block:: bash

   gatk Mutect2 -R $REF -I tumor.bqsr.bam -tumor TUMOR -I normal.bqsr.bam -normal NORMAL --germline-resource $GNOMAD_AF --panel-of-normals $PON --intervals $BAIT_BED -O tumor_vs_normal.unfiltered.vcf.gz

   gatk FilterMutectCalls -R $REF -V tumor_vs_normal.unfiltered.vcf.gz --contamination-table contamination.table -O tumor_vs_normal.filtered.vcf.gz

   # Explanation: Mutect2 (GATK) for somatic calls (Yang/Ljungström). Use PoN for artifacts. Adjust for low-VAF subclones.

**5) CNAs (CNVkit & FACETS)**

Detect CNAs.

.. code-block:: bash

   # CNVkit (Yang)
   cnvkit.py batch tumor.bqsr.bam --normal normal.bqsr.bam --targets $BAIT_BED --fasta $REF --output-reference cnvkit_ref.cnn --output-dir cnvkit_out
   cnvkit.py call cnvkit_out/tumor.bqsr.cns -o cnvkit_out/tumor.call.cns

   # FACETS (Ljungström)
   Rscript facets.R --tumor tumor.bqsr.bam --normal normal.bqsr.bam --output facets_out

   # Explanation: CNVkit for segmentation; FACETS for purity/ploidy in subclonal analysis. Quantify CNA-ITH (Yang).

**6) Clonality (PyClone-VI & SciClone)**

Infer subclones.

.. code-block:: bash

   # SciClone (Ljungström)
   Rscript -e "library(sciClone); sciClone(vafs=..., copyNumberCalls=..., sampleNames=..., outputFile='clones')"

   # PyClone-VI (modern)
   pyclone-vi fit --in-file pyclone_input.tsv --out-file pyclone_results.tsv --seed 13

   # Explanation: SciClone for VAF clustering; PyClone-VI for Bayesian ITH (update). Use multi-timepoint (CLL) or regions (liver).

**7) Annotation (ANNOVAR & VEP)**

Annotate variants.

.. code-block:: bash

   # ANNOVAR (both studies)
   perl table_annovar.pl filtered.vcf humandb/ -buildver hg38 -out annotated -protocol refGene,cosmic70,clinvar -operation g,f,f

   # VEP (modern)
   vep -i filtered.vcf.gz -o annotated.vep.vcf --assembly GRCh38 --cache --offline --vcf --everything --fork 16

   # Explanation: ANNOVAR for functional impact; VEP for detailed consequences. Handle branch-private variants.

**8) Mutation Signatures (SigProfiler & DeconstructSigs)**

Extract signatures.

.. code-block:: bash

   # DeconstructSigs (Ljungström)
   Rscript -e "library(deconstructSigs); whichSignatures(sigs.input=..., signatures.ref=signatures.cosmic)"

   # SigProfiler (modern)
   python -m sigProfilerExtractor sigprofiler_extract input output project

   # Explanation: DeconstructSigs for known signatures; SigProfiler for de novo. Analyze per-subclone (ITH-aware).

**9) Transcriptomics (STAR, Salmon, edgeR, DESeq2)**

Process RNA-seq.

.. code-block:: bash

   # STAR alignment
   STAR --runThreadN 16 --genomeDir star_index --readFilesIn rna_R1.fq.gz rna_R2.fq.gz --outSAMtype BAM SortedByCoordinate

   # Salmon quantification
   salmon quant -i transcripts_index -l A -1 rna_R1.fq.gz -2 rna_R2.fq.gz -o quant_out --gcBias

   # DE with DESeq2
   Rscript -e "library(DESeq2); dds <- DESeqDataSetFromMatrix(countData=..., colData=..., design=~condition); dds <- DESeq(dds); res <- results(dds)"

   # Explanation: STAR/Salmon for alignment/quant (Yang). edgeR/DESeq2 for signatures/DE. Integrate for multi-omic ITH.

**10) Reporting (MultiQC & Custom)**

Aggregate results.

.. code-block:: bash

   multiqc . -o multiqc_report --verbose

   # Custom R for phylogenies/signatures
   Rscript visualize.R

   # Explanation: MultiQC for QC aggregation. Custom scripts for ITH plots (MEGA, ggplot2 from studies).

Reproducibility Checklist
~~~~~~~~~~~~~~~~~~~~~~~~~

- Use Conda/Docker for envs (e.g., GATK Docker).
- Record versions (conda env export > env.yml).
- Pin references (MD5 checksums).
- Emit metrics (dup rate, coverage, TMB).
- Store logs/MultiQC in run folder. For ITH, track per-region metrics.