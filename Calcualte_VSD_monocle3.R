#library(monocle3)
library (DESeq2)

uniref_counts <- read.csv(file="../Filtering/Bowtie2_map_uniRefGeneCounts_Filter7.csv", sep = ",", row.names = 1) 
#uniref_counts <- read.csv(file="CPM_Bowtie2_map_uniRefGeneCounts_Filter7.csv", sep = ",", row.names = 1) #CPM counts of UniRefs which are present in all Patients # Last column in a Sum

#--Sample Size Less then 100 MB
FileSize <- read.csv(file="FileSize.txt", sep = "\t", header = T) # For all R1 and R2
FileSize_filtered <- subset(FileSize, SizeMB > 100)
SamplesSize <- unique(FileSize_filtered$FileName) # 1276 Samples have higher then 100 MB Size and 62 Samples lesser then 100 MB are removed

colnames(uniref_counts)
uniref_counts_filtered <- uniref_counts[SamplesSize]

#-- To remove samples which have zero UniRefs (11 Samples to be removed)
Sample_sum <- colSums(uniref_counts_filtered)
Sample_sum <- data.frame(Sample_sum)
row.names(Sample_sum) <- colnames(uniref_counts_filtered)
Samples_to_be_retained <- row.names (subset(Sample_sum, Sample_sum > 10000)) #These are the samples which have sum of UniRef Count more than 10k

uniref_counts_filtered <- uniref_counts_filtered[Samples_to_be_retained]

#######################################################
#Check how much read left
Total_unirefCounts_Filter7 <- data.frame(colSums(uniref_counts))
#--Load Stats of complete data
total_Stats <- read.csv(file="GeneCounts_Stats.csv", sep = ",", row.names = 1)# Total Count stats on All UniRef (complete data file) without Filtering
Total_unirefCounts <- data.frame(total_Stats$TotalUnirefCount)
rownames(Total_unirefCounts) <- rownames(total_Stats)

Total_Counts <- merge(Total_unirefCounts, Total_unirefCounts_Filter7, by=0, all=F)
rownames(Total_Counts) <- Total_Counts$Row.names; Total_Counts$Row.names <- NULL
colnames(Total_Counts) <- c("totalMappedReads","AfterFiltering")
Total_Counts_filtered <- data.frame(Total_Counts[rownames(Total_Counts) %in% SamplesSize, ])
Proportion_Left <- data.frame(Total_Counts_filtered$AfterFiltering/Total_Counts_filtered$totalMappedReads)*100
rownames(Proportion_Left) <- rownames(Total_Counts_filtered)
#Proportion_Left_Filtered <- data.frame(Proportion_Left[rownames(Proportion_Left) %in% SamplesSize, ])  # Extract rows from data

########################################################
# Use SD/AVE for all uniRefs across Patients --> Use least variable UniRefs for first UMAP and then Keep trying
#-- Load Metadata
meta <- read.csv(file="hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)
meta <- meta[order(row.names(meta)), ]
row.names(meta)
meta_filtered <- meta [rownames(meta) %in% SamplesSize, ]
meta_filtered <- meta_filtered [rownames(meta_filtered) %in% Samples_to_be_retained, ]

  
#pd <- new("AnnotatedDataFrame", data = meta_filtered)

# To reorder: https://hbctraining.github.io/Intro-to-R-flipped/lessons/09_reordering-to-match-datasets.html
genomic_idx <- match(rownames(meta_filtered), colnames(uniref_counts_filtered))
uniref_counts_ordered  <- uniref_counts_filtered[ , genomic_idx]

#--- Data Sparsity
sum(uniref_counts_ordered == 0)/(dim(uniref_counts_ordered)[1]*dim(uniref_counts_ordered)[2])
#[1] 0.06168271

UniRef_names <- row.names(uniref_counts_ordered)
UniRef_names <- data.frame(UniRef_names)
row.names(UniRef_names) <- UniRef_names$UniRef_names
colnames(UniRef_names) <- "gene_short_name"

#---- Try the same Monocle3 on VSD data
library(DESeq2)
Deseq_uniref_counts_ordered <- DESeqDataSetFromMatrix(uniref_counts_ordered, colData=meta_filtered,design= ~Participant.ID)
Deseq_uniref_counts_ordered <- DESeq(Deseq_uniref_counts_ordered)

Deseq_uniref_counts_ordered_vsd <- varianceStabilizingTransformation(Deseq_uniref_counts_ordered,blind=FALSE)
Deseq_uniref_counts_ordered_vsd_df <- data.frame(assay(Deseq_uniref_counts_ordered_vsd))
write.csv(Deseq_uniref_counts_ordered_vsd_df, file = "Deseq_uniref_counts_ordered_vsd_df.csv")
save.image("Filter7_VSD.RData")

#UniRef_counts_Dataset <- new_cell_data_set(as.matrix(uniref_counts_ordered),
#                                           cell_metadata = meta_filtered,
#                                           gene_metadata = UniRef_names)

#UniRef_counts_Dataset <- new_cell_data_set(as.matrix(Deseq_uniref_counts_ordered_vsd_df),
 #                                          cell_metadata = meta_filtered,
  #                                         gene_metadata = UniRef_names)


#UniRef_counts_Dataset <- preprocess_cds(UniRef_counts_Dataset, num_dim = 50)
#UniRef_counts_Dataset <- preprocess_cds(UniRef_counts_Dataset, num_dim = 50, norm_method = "none")

#UniRef_counts_Dataset <- align_cds(UniRef_counts_Dataset, alignment_group = "Participant.ID") # To remove batch effect

## Step 3: Reduce the dimensions using UMAP
#UniRef_counts_Dataset <- reduce_dimension(UniRef_counts_Dataset)

#plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE,  color_cells_by = "Participant.ID")

#plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE,  color_cells_by = "diagnosis")

#Cluster the cells
#UniRef_counts_Dataset <- cluster_cells(UniRef_counts_Dataset)
#plot_cells(UniRef_counts_Dataset, color_cells_by = "partition")

## Step 5: Learn a graph
#UniRef_counts_Dataset <- learn_graph(UniRef_counts_Dataset)

## Step 6: Order cells
# UniRef_counts_Dataset <- order_cells(UniRef_counts_Dataset) Don't run this command if you are Running this on Cluster

# a helper function to identify the root principal points:
#get_earliest_principal_node <- function(UniRef_counts_Dataset, time_bin="CD"){
#  cell_ids <- which(colData(UniRef_counts_Dataset)[, "diagnosis"] == time_bin)
  
#  closest_vertex <-
#    UniRef_counts_Dataset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#  closest_vertex <- as.matrix(closest_vertex[colnames(UniRef_counts_Dataset), ])
#  root_pr_nodes <-
#    igraph::V(principal_graph(UniRef_counts_Dataset)[["UMAP"]])$name[as.numeric(names
#                                                                                (which.max(table(closest_vertex[cell_ids,]))))]
  
#  root_pr_nodes
#}
#UniRef_counts_Dataset <- order_cells(UniRef_counts_Dataset, root_pr_nodes=get_earliest_principal_node(UniRef_counts_Dataset))

#jpeg("VSD_Filter7_Trajectory_diagnosis.jpg", height = 5, width = 5, units = 'in', res = 600)
#plot_cells(UniRef_counts_Dataset, color_cells_by = "diagnosis")
#dev.off ()
#jpeg("VSD_Filter7_Trajectory_pseudotime.jpg", height = 5, width = 5, units = 'in', res = 600)
#plot_cells(UniRef_counts_Dataset, color_cells_by = "pseudotime")
#dev.off ()
#jpeg("VSD_Filter7_Trajectory_PID.jpg", height = 7, width = 7, units = 'in', res = 600)
#plot_cells(UniRef_counts_Dataset, color_cells_by = "Participant.ID")
#dev.off ()
