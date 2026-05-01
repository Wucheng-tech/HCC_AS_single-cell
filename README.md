# Title: A depth-enhanced single-cell framework reveals cell-type-specific alternative splicing programs in hepatocellular carcinoma associated with tumorigenesis and sorafenib resistance
Alternative splicing (AS) is a major source of transcriptomic diversity and has been increasingly implicated in cancer progression, including hepatocellular carcinoma (HCC). However, accurately comparing AS variation within defined cell types at single-cell resolution remains challenging due to limited transcript coverage and data sparsity. To address this gap, we developed scSplicing, a depth-enhanced analytical framework that integrates cell-type-specific pseudo-bulk aggregation with gene expression similarity-guided dynamic KNN imputation for robust splicing inference from sparse single-cell data. Using scFAST-seq, a full-length single-cell transcriptome platform, we applied scSplicing to systematically characterize AS landscapes in HCC.

Importantly, scSplicing is broadly applicable beyond scFAST-seq and can be readily adapted to droplet-based single-cell RNA-seq platforms (e.g., 10x Genomics) as well as other sequencing technologies, providing a generalizable strategy for cell-type-resolved alternative splicing analysis across diverse physiological and disease contexts.
# Usage
We provide an overview of the associated analysis workflow, which includes the construction of AS profiles from bulk RNA-seq data (PSI matrices,HCC_bulk_rmats_PSI.sh, Huh7_SR_rmats_PSI.sh), as well as AS profiling from scFAST-seq data. 

The single-cell AS analysis workflow involves the following steps:

1，Cell type annotation based on gene expression profiles – Assigning each single cell to a defined cell type using transcriptomic signatures. (step0_scFAST_generate_metadata.R)

2，Splitting BAM files by cell identity – Dividing raw single-cell BAM files according to annotated cell types. (step1_scFAST_prepare_cell_bam.sh)

3，Pseudo-bulk binning by sample and cell type – Aggregating BAM files of the same sample and cell type into bins to increase coverage for AS analysis. (step2_scFAST_merge_bam_by_celltype_bin.sh） 

4，AS quantification using rMATS – Calculating PSI values for splicing events in each pseudo-bulk bin. (step3_rmats_compare.sh，step4_rmats_filter.sh，step5_rmats_PSI.sh) 

5，Imputation of missing values – Applying similarity-based strategies to fill in sparse AS measurements and reduce data dropout. (AS_imputation_dynamic.py)

6，Single-cell differential AS analysis – Performing differential alternative splicing analysis between tumor and normal cells within defined cell types.  (AS_FindASMarkers.r)

This workflow enables systematic, cell type-specific investigation of AS patterns at single-cell resolution and facilitates the identification of splicing events associated with tumorigenesis and chemoresistance in HCC.

# Author
Cheng Wu, Wenliang Zhang.

# Contact
Cheng Wu wucheng820@163.com

