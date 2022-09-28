library(monocle3)
uniref_counts <- read.csv(file="Final_filtered_UniRef20k_and_SamplesLowdepth_proportionLeft.csv", sep = ",", row.names = 1) 
#uniref_counts <- read.csv(file="CPM_Bowtie2_map_uniRefGeneCounts_Filter7.csv", sep = ",", row.names = 1) #CPM counts of UniRefs which are present in all Patients # Last column in a Sum
#-- Samples less then 50% count proportion left removed
#-- Unirefs removed whose total count is less than 20k

########################################################
# Use SD/AVE for all uniRefs across Patients --> Use least variable UniRefs for first UMAP and then Keep trying
#-- Load Metadata
meta <- read.csv(file="metaTrascriptomic.csv", header = T, row.names = 1)
Samples_to_be_retained <- colnames(uniref_counts)
meta_filtered <- meta[rownames(meta) %in% Samples_to_be_retained, ]
samples_in_meta <- rownames(meta_filtered)

#--to get same samples are were in the metadata file
uniref_counts <- uniref_counts[samples_in_meta]

# To reorder: https://hbctraining.github.io/Intro-to-R-flipped/lessons/09_reordering-to-match-datasets.html
genomic_idx <- match(rownames(meta_filtered), colnames(uniref_counts))
uniref_counts_ordered  <- uniref_counts[ , genomic_idx]

#--- Data Sparsity
sum(uniref_counts_ordered == 0)/(dim(uniref_counts_ordered)[1]*dim(uniref_counts_ordered)[2])
#0.359398

UniRef_names <- row.names(uniref_counts_ordered)
UniRef_names <- data.frame(UniRef_names)
row.names(UniRef_names) <- UniRef_names$UniRef_names
colnames(UniRef_names) <- "gene_short_name"

UniRef_counts_Dataset <- new_cell_data_set(as.matrix(uniref_counts_ordered),
                                           cell_metadata = meta_filtered,
                                           gene_metadata = UniRef_names)


UniRef_counts_Dataset <- preprocess_cds(UniRef_counts_Dataset, num_dim = 50)
#UniRef_counts_Dataset <- preprocess_cds(UniRef_counts_Dataset, num_dim = 50, norm_method = "none") # When VSD data was Used

#UniRef_counts_Dataset <- align_cds(UniRef_counts_Dataset, alignment_group = "Participant.ID") # To remove batch effect

## Step 3: Reduce the dimensions using UMAP
UniRef_counts_Dataset <- reduce_dimension(UniRef_counts_Dataset)
plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE,  color_cells_by = "Participant.ID")
plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE,  color_cells_by = "diagnosis")

#Cluster the cells
UniRef_counts_Dataset <- cluster_cells(UniRef_counts_Dataset)
plot_cells(UniRef_counts_Dataset, color_cells_by = "partition")

## Step 5: Learn a graph
UniRef_counts_Dataset <- learn_graph(UniRef_counts_Dataset)

## Step 6: Order cells
# UniRef_counts_Dataset <- order_cells(UniRef_counts_Dataset) Don't run this command if you are Running this on Cluster

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(UniRef_counts_Dataset, time_bin="CD"){
  cell_ids <- which(colData(UniRef_counts_Dataset)[, "diagnosis"] == time_bin)
  
  closest_vertex <-
    UniRef_counts_Dataset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(UniRef_counts_Dataset), ])
  root_pr_nodes <-
    igraph::V(principal_graph(UniRef_counts_Dataset)[["UMAP"]])$name[as.numeric(names
                                                                                (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
UniRef_counts_Dataset <- order_cells(UniRef_counts_Dataset, root_pr_nodes=get_earliest_principal_node(UniRef_counts_Dataset))

jpeg("Trajectory_diagnosis.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE, color_cells_by = "diagnosis")
dev.off ()
jpeg("Trajectory_pseudotime.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE, color_cells_by = "pseudotime")
dev.off ()
jpeg("Trajectory_PID.jpg", height = 7, width = 7, units = 'in', res = 600)
plot_cells(UniRef_counts_Dataset, label_groups_by_cluster=FALSE, color_cells_by = "Participant.ID")
dev.off ()
