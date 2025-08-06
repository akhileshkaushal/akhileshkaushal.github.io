Case Study: ITH-aware Somatic Interpretation for WGS/WES
========================================================

**Focus**: Sampling strategy, when to prioritize tumor-normal over tumor-only, ITH-aware filtering and reporting, RNA-based response signatures, and mechanism-informed combinations (c-Jun/JNK/β-catenin axis).

**Primary reference**:  
Yang *et al.*, "Pharmacogenomic profiling of intra-tumor heterogeneity using a large organoid biobank of liver cancer," *Cancer Cell* (2024), DOI: 10.1016/j.ccell.2024.03.004.

This seminal study, published in April 2024, has since been cited in emerging research on organoid-based models for liver cancer, including works exploring pharmacogenomic vulnerabilities and advances in tumor heterogeneity modeling.



 The authors established a biobank of 399 organoids from 144 primary liver cancer (PLC) patients, quantified intra-tumor heterogeneity (ITH), developed predictive multi-gene expression signatures for drug responses, and elucidated c-Jun-mediated resistance to lenvatinib, proposing a rational combination therapy strategy.

What Was Built
--------------

**Living PLC organoid biobank**

Intra-tumor heterogeneity (ITH) refers to the genetic, epigenetic, and phenotypic variations within a single tumor, arising from subclonal evolution and environmental pressures. To model this in primary liver cancer (PLC, encompassing hepatocellular carcinoma [HCC], intrahepatic cholangiocarcinoma [ICC], and combined HCC-ICC [cHCC-ICC]), the researchers created a "living biobank" of patient-derived organoids (PDOs).

- **Scale**: 399 tumor organoids derived from 144 patients, sourced from 528 multi-regional samples of surgical resections (e.g., center, edge, and intermediate regions of tumors).
- **Sampling**: Multi-regional approach to capture spatial ITH, including inter-tumor (between patients) and intra-tumor (within the same tumor) variations. This enabled simultaneous genomic and functional analyses.
- **Recapitulation**: Organoids faithfully mirrored parental tumors in histology (e.g., H&E staining showing trabecular patterns in HCC), genomics (mutations and copy number alterations [CNAs]), and transcriptomics (high correlations in RNA-seq profiles for markers like AFP, GPC3 for HCC).
- **Utility**: Ideal for high-throughput drug sensitivity screening, mechanistic studies (e.g., knockdown/overexpression), and xenograft validations, bridging ex vivo modeling with clinical phenotypes.

**Establishment details**

Organoid culture success is critical for biobanking, with rates varying by tumor type and processing.

- Overall establishment rate: **75.6%** (399/528 regions), higher for HCC (80.5%) than ICC (64.5%) or cHCC-ICC (68.2%).
- Success factors:
  - Proportion of viable cells at tissue acquisition (assessed via trypan blue exclusion).
  - Time from resection to processing (ideally <2 hours to minimize ischemia).
  - Enzymatic digestion method (e.g., collagenase-based for gentle dissociation).

- Culture protocol (based on established methods like Broutier et al., 2016):
  - Basal medium: Advanced DMEM/F12 supplemented with penicillin/streptomycin, Glutamax, HEPES, B27 (without vitamin A), N2, N-acetyl-L-cysteine, nicotinamide.
  - Growth factors: Recombinant human EGF, FGF10, HGF, gastrin, forskolin.
  - Wnt pathway agonists: Noggin, R-spondin1 (Rspo1), Wnt3a.
  - Inhibitors: A83-01 (TGF-β), Y-27632 (ROCK).
  - Embedded in 5-10% Matrigel; passaged every 1-3 weeks via mechanical/trypsin dissociation.
  - Cryopreservation: Viable recovery post-thawing, enabling long-term storage.

This biobank (detailed in Figure 1) allowed comprehensive profiling, revealing ITH's impact on drug resistance.

Cohort Profiling and Heterogeneity
----------------------------------

**Concordance and ITH**

PDOs showed strong fidelity to tumors:

- Histological: IHC for markers like HepPar1, AFP (HCC), KRT19 (ICC).
- Genomic: High concordance in variant allele frequencies (VAFs) for COSMIC cancer genes (e.g., TP53, CTNNB1).
- Transcriptomic: Pearson correlations >0.8 between paired tumor-organoid RNA-seq (n=99 pairs; Figure 1H).

Multi-regional sampling uncovered extensive ITH, with higher mutation/CNA ITH linked to poorer overall survival (OS; hazard ratios via Cox proportional hazards, p<0.05; Figure 2E).

**Trunk/branch structure**

Tumors evolve via branching phylogenies, where "trunk" mutations are shared (early events) and "branch" mutations are region-specific (late events).

- **Trunk-dominant tumors**: 12/32 patients (trunk ratio [TR] >1), with most drivers (e.g., TP53) in the trunk.
- **Branch-dominant tumors**: 20/32 patients (TR ≤1), featuring region-private alterations (e.g., divergent drivers in different organoids; Figure 2B).
- **Clinical implication**: Branch-dominant ITH correlated with lenvatinib resistance (e.g., variable target expression like FGFRs; Figure 2F) and worse OS. Single-biopsy analyses may miss branch events, underestimating resistance.

**Genomic Profiling Methods – Whole-Exome Sequencing (WES)**

Whole-exome sequencing (WES) targets the protein-coding regions (~1–2% of the genome) and is well-suited for detecting somatic mutations and copy-number alterations (CNAs) in cancer driver genes.  
In this study, tumor–normal paired WES was emphasized to accurately distinguish somatic from germline variants.

**WES Sequencing**

- **Samples**: Genomic DNA was extracted from:
  - Organoids (**n** = 399)
  - Matched tumors (subset)
  - Matched normal samples from adjacent liver or blood (**n** = 144)
- **Library preparation**: Exons captured using *Agilent SureSelectXT Human All Exon V6* (~60 Mb target size)
- **Sequencing**: Illumina NovaSeq 6000 platform  
  - 150 bp paired-end reads  
  - Target coverage: ~100–200× to enable detection of low-VAF subclones

**Bioinformatics Pipeline** (see STAR Methods for full reproducibility)

- **Alignment**:  
  - `bwa-mem2` (v2.0) to GRCh38/hg38 reference  
  - Duplicate reads marked with *GATK MarkDuplicates*

- **Variant Calling**:  
  - *GATK* (v4.1.2.0) best practices  
  - Somatic SNVs/indels called with *MuTect2* in tumor–normal mode  
  - Germline filtering via Panel of Normals (PoN)  
  - Post-calling filters: VAF ≥ 0.05, depth ≥ 20, COSMIC annotation presence

- **CNA Detection**:  
  - *CNVkit* (v0.9.7.b1) – read-depth segmentation with CBS algorithm  
  - *GISTIC2* (v2.0) – focal and arm-level CNA significance (q < 0.25)  
  - CNA intratumoral heterogeneity (CNA-ITH) = variance in copy ratios across regions

- **Phylogenetic Reconstruction**:  
  - *MEGA5* – neighbor-joining trees based on shared/private mutation presence/absence

- **ITH Quantification**:  
  - Mutation-ITH = fraction of private mutations  
  - CNA-ITH = Jaccard distance on CNA segments  
  - TR (trunk ratio) = trunk / branch mutations  
  - Correlations: Spearman’s ρ showed positive association between mutation- and CNA-ITH (p < 0.05)  
  - Survival analysis: Kaplan–Meier with log-rank test (two-sided)

- **Advantages of Tumor–Normal Design**:  
  - Reduces false positives (e.g., via *GATK CalculateContamination*)  
  - Improves VAF accuracy for subclonal inference  
  - Tumor-only fallback possible using PoN/germline resources, but may miss low-VAF events

**Key Findings from WES**

- Recovered known PLC drivers:  
  - *TP53* (~50%)  
  - *CTNNB1* (~30%)
- ITH levels alone were not predictive of specific drug responses
- Branch-dominant phylogenies were enriched for resistance phenotypes
- Findings support multi-region sampling and complementary transcriptomic profiling


Drug Screening and Response Definition
---------------------------------------

Drug screening in PDOs assesses functional ITH, where phenotypic responses (e.g., viability) vary across regions.

**Agents tested**

- Clinically relevant: Four tyrosine kinase inhibitors (TKIs; lenvatinib, sorafenib, regorafenib, apatinib), anti-VEGF (bevacizumab), FGFR inhibitor (pemigatinib), IDH1 inhibitor (ivosidenib).
- Mechanism-focused: c-Jun inhibitors (veratramine, SR11302, NY2267); novel conjugate PKUF-01 (lenvatinib-veratramine; synthesized in-house).
- Additional: Chemotherapy (gemcitabine, cisplatin) for context.

**Response metrics**

- Assays: Organoids seeded in 96/384-well plates (ultra-low attachment; ~100 organoids/well in 5% Matrigel). 7-point dose-response (5-10-fold dilutions; 72h exposure). Viability via CellTiter-Glo (ATP-based luminescence).
- Metrics: IC50 (half-maximal inhibitory concentration); normalized area under curve (AUC; 0-1 scale, lower = sensitive). High correlation (Spearman ρ >0.79; Figure S3A).
- Patient-level: "Worst-region" rule—most resistant organoid dictates sensitivity (reflecting ITH-driven escape; Figure 3D).
- Thresholds: Percentile-based on clinical objective response rates (ORRs; e.g., lenvatinib 24.1% → top 24.1% sensitive; Figure 3E).

**Passage stability**

- Tested in 16 organoids: 12/16 stable lenvatinib sensitivity (early vs. late passages; Figure S3B).
- 4/16 showed drift (resistance gain), attributable to subclone selection—underscoring need for early-passage testing.

**Clinical alignment**

- Validated in 14 relapsed patients: Organoid sensitivities matched outcomes (e.g., lenvatinib responders sensitive in all regions; progressors had resistant regions; Figure S3H).
- Xenografts: Multi-region PDOs from resistant patients confirmed heterogeneous tumor growth under lenvatinib (Figure 3H).

Expression Signatures Predicting Response
------------------------------------------

Transcriptomic modeling integrates RNA-seq with drug data for predictive biomarkers, as genomic ITH alone lacked specificity.

- **RNA-seq**: Illumina paired-end; Salmon/edgeR for quantification/differential expression. Variance analysis (ANOVA) for ITH.

**Lenvatinib signature**

- Cohorts: Training (n=106 organoids), validation (n=106).
- Workflow: 254 genes correlated with AUC (Spearman, p<0.05; e.g., positive: JUN, WNT pathway; negative: cell cycle genes; Figure 4B,C).
- Refined: 13-gene panel via machine learning (elastic net regression; STAR Methods). AUROC 0.86 (train), 0.81 (valid; Figure 4F).
- External: TCGA-LIHC (n=371) clustered into ~25% "sensitive" (matching ORR; unsupervised k-means; Figure 4D).

**Other TKIs**

- Sorafenib: 199-gene signature, AUROC >0.9 (Figure 4L).
- Regorafenib/Apatinib: Similar (AUROCs 0.8-0.7).
- Alignment: Predicted clinical responses in relapsed cases (e.g., signature scores vs. progression; Figure 4L).

Signatures (Table S5) enable patient stratification, with ITH-aware aggregation (worst-region score).

Mechanism of Lenvatinib Resistance & Combination Strategy
----------------------------------------------------------
Lenvatinib targets VEGFR/FGFR; resistance often via bypass signaling.

**Upstream targets**

- FGFR1-4 highly expressed (RNA-seq); siRNA knockdown reduced sensitivity (Figure S5A).

**c-Jun axis**

- c-Jun (transcription factor) overexpressed in resistant PDOs (Figure 5A).
- Manipulation: shRNA knockdown sensitized (IC50 drop in 21/6 resistant organoids; Figure 5G); overexpression conferred resistance (Figure 5S).
- Upstream: JNK (activator) and Wnt/β-catenin (stabilizer). Knockdown/activation modulated sensitivity reciprocally (Figures 5K, S5).
- Mechanism: c-Jun activates JNK/β-catenin loop, bypassing FGFR inhibition (Western blots for phospho-c-Jun; Figure 5M).

**PKUF-01**

- Design: Lenvatinib-veratramine conjugate (dual FGFR/c-Jun targeting; synthesis in STAR Methods).
- Efficacy: Sensitized 27/144 patients (59% of resistant; Figure S5E); reduced IC50s ex vivo; superior xenograft control (Figure 6H).
- Limitation: Poor oral bioavailability in vivo, suggesting need for formulation optimization.

ITH-aware Somatic Interpretation & Recommendations
---------------------------------------------------

**Why ITH matters**

ITH drives therapy failure; multi-region analysis dropped predicted lenvatinib benefit from 72.9% (per-organoid) to 37.1% (per-patient; Figure 3I).

**Sampling strategy**

- Multi-region (≥3 sites/tumor) plus metastases if available.
- Prioritize tumor-normal WES for somatic confidence (vs. tumor-only, which relies on databases like gnomAD for germline filtering).

**ITH-aware filtering**

- Consolidate multi-region VCFs; classify variants (trunk/branch/private) via phylogenetics (e.g., MEGA5).
- "Worst-region" for sensitivity; integrate WES/RNA-seq (e.g., VAF-CN-expression triangles for subclone validation).

**RNA-based response signatures**

- Pre-treatment RNA-seq from WES regions; compute scores per-region, aggregate conservatively.
- Research-use only; enrich for trials (e.g., NCT numbers for TKI combos).

**Mechanism-aware reporting**

- Flag c-Jun/JNK/β-catenin activation (e.g., via IHC/RNA).
- Recommend combinations (lenvatinib + JNK inhibitors like SP600125; PKUF-01 analogs).
- Include trials (e.g., lenvatinib combos in HCC).

**Implementation tips**

- Matched normals for GATK; maintain PoN/germline VCFs.
- Log versions (e.g., GATK v4.1.2.0) for reproducibility; containerize (Docker/Singularity).

Limitations
-----------

- Small clinical validation (n=14 relapsed); signatures research-grade, require prospective trials.
- Passage drift in ~25% lines (subclonal evolution).
- PKUF-01 bioavailability issues.
- WES limits: Misses non-coding variants; WGS could enhance for structural variants/intratumoral fusions. Recent studies advocate integrated WGS/RNA for fuller ITH mapping.

