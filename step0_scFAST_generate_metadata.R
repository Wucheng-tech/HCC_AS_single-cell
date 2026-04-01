# Step 0 — Generate cell metadata from Seurat
############################################################
library(Seurat)
library(plyr)
HCC <- readRDS("/public/home/wucheng/Analysis/scFAST/RDS/HCC_scFAST.rds")
# Map cluster ID to major cell type
cluster_from <- c("2","1","0","7","10","3","6","8","9","4","5")
cluster_to   <- c("Hepatocyte","Hepatocyte","T","T","NK","Mye","B","pB","Fib","Endo","Pro")
HCC@meta.data$celltype <- plyr::mapvalues(
    x = HCC@meta.data$Cluster,
    from = cluster_from,
    to = cluster_to
)
table(HCC$celltype)
#Hepatocyte          T        Pro         NK        Mye          B         pB 
#     62210      47669       6100       4095      22171       5908       5137 
#       Fib       Edno 
#      4344      10024 

meta <- HCC@meta.data[,c("Sample","Lesion_size","celltype","sub_celltype")] # Extract required metadata
meta$barcode_full <- rownames(meta) # Add original cell barcode
meta$cell <- substr(meta$barcode_full, 1, 17) # barcode (first 17 bp)
# Map lesion ID to unified patient ID
lesion_from <- c("S243091_N","2_S243091_M","6_S243091_T1","S262103_N","7_S262103_T1","5_S262103_M","S264663_N","10_S264663_T1","1_S264663_M","S267661_N","3_S267661_T1","4_S267661_T2","S267835_N","9_S267835_T1","8_S267835_T2")
lesion_to <- c("P243091_PN","P243091_M","P243091_T","P262103_N","P262103_T1","P262103_M1","P264663_N","P264663_T","P264663_M","P267661_N","P267661_T1","P267661_T2","P267835_N","P267835_T1","P267835_T2")
meta$ID <- plyr::mapvalues(
    x = meta$Lesion_size,
    from = lesion_from,
    to = lesion_to
)
meta_out <- meta[,c("Sample","celltype","sub_celltype","cell","ID")]
write.table( meta_out, "/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/sample_cell.txt",sep="\t",quote=FALSE,row.names=FALSE)
head(meta_out)
#                    Sample   celltype sub_celltype              cell        ID
#CGACAGCCGAAAGCAAT_1 P01_T2          T        Other CGACAGCCGAAAGCAAT P264663_M
#GGCAGTCATGACTAAGT_1 P01_T2 Hepatocyte        Hep_N GGCAGTCATGACTAAGT P264663_M
#AAAGGTACGATTGTGCA_1 P01_T2          T        Other AAAGGTACGATTGTGCA P264663_M
#AACCACACGATACTGAC_1 P01_T2          T        Other AACCACACGATACTGAC P264663_M
#ACCCTTGTACTGCAGCG_1 P01_T2        Mye        Other ACCCTTGTACTGCAGCG P264663_M
#ACTTATCATGGAAGTCC_1 P01_T2 Hepatocyte        Hep_M ACTTATCATGGAAGTCC P264663_M
