.. WES_CLL documentation master file

Bioinformatics Pipeline for Whole-Exome Sequencing (WES) and Whole-Genome Sequencing (WGS) in Cancer Genomics
===========================================================================

This documentation summarizes bioinformatics workflows, updated WGS/WES pipeline recommendations, and key findings from two pivotal studies: Yang et al. (2024) on pharmacogenomic profiling of intra-tumor heterogeneity in liver cancer using organoid biobanks, and Ljungström et al. (2016) on whole-exome sequencing in relapsing chronic lymphocytic leukemia revealing recurrent RPS15 mutations. It includes study contexts, major results, and interpretation guidelines for ITH-aware somatic variant detection, clonal evolution analysis, and clinical implications in heterogeneous cancers.

.. note::
	The original studies employed sequencing pipelines tailored to their experimental contexts. Yang et al. (2024) utilized modern tools like bwa-mem2 and GATK v4 for multi-region organoid profiling in primary liver cancer (PLC), emphasizing transcriptomic integration for drug response signatures. Ljungström et al. (2016) used earlier versions (e.g., BWA v0.7, GATK v2) for paired pre/post-relapse WES in CLL, focusing on subclonal shifts and novel drivers like RPS15.
	
	In this documentation, workflows have been adapted and modernized for general-purpose WGS/WES bioinformatics, incorporating current best practices (as of 2025) such as BWA-MEM2 alignment to GRCh38, GATK4 Mutect2 for somatic calling, FACETS for purity/ploidy estimation, SciClone for subclonal inference, and lower VAF thresholds (e.g., ≥0.05) to enhance sensitivity, reproducibility, and applicability to cancers with high ITH or evolutionary dynamics.


Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Sections:
   
   
   intro
   run_pipeline
   advanced_analysis
   annotated_example
   structure_and_containerisation
   wes_cll_methods
   casestudy1
   casestudy2
   references


