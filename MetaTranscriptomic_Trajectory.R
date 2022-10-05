#-- Get UniRef Counts (MetaTranscitpomic Data)
uniref_counts <- read.csv (file="~/Library/CloudStorage/Box-Box/David_Casero_Lab/HMP_IBD_Project/Count_Trajectory/New_count_Trajectory_Mapped_Counts/MetaTranscriptomic_Data/Final_filtered_UniRef20k_and_SamplesLowdepth_proportionLeft.csv", sep=",", row.names =1, header = T)

#-- Get metadata
meta <- read.csv(file="~/Library/CloudStorage/Box-Box/David_Casero_Lab/HMP_IBD_Project/Count_Trajectory/New_count_Trajectory_Mapped_Counts/MetaTranscriptomic_Data/metaTrascriptomic.csv", header = T, row.names = 1)
Samples_to_be_retained <- colnames(uniref_counts)
meta_filtered <- meta[rownames(meta) %in% Samples_to_be_retained, ]

#--to get same samples are were in the metadata file
samples_in_meta <- rownames(meta_filtered)
uniref_counts <- uniref_counts[samples_in_meta]

# To reorder: https://hbctraining.github.io/Intro-to-R-flipped/lessons/09_reordering-to-match-datasets.html
genomic_idx <- match(rownames(meta_filtered), colnames(uniref_counts))
uniref_counts_ordered  <- uniref_counts[ , genomic_idx]
#--- Data Sparsity
sum(uniref_counts_ordered == 0)/(dim(uniref_counts_ordered)[1]*dim(uniref_counts_ordered)[2])
#0.359398


#- Get gene name File (here it is taxa names)
UniRef_names <- row.names(uniref_counts_ordered)
UniRef_names <- data.frame(UniRef_names)
row.names(UniRef_names) <- UniRef_names$UniRef_names
colnames(UniRef_names) <- "gene_short_name"
fd <- new("AnnotatedDataFrame", data = UniRef_names)

#--- Format metadata of Phenodata
pd <- new("AnnotatedDataFrame", data = meta_filtered)

#uniref_counts_ordered_format <- lapply(uniref_counts_ordered, as.numeric)
#uniref_counts_ordered_format <- data.frame(uniref_counts_ordered_format)
#rownames(uniref_counts_ordered_format) <- rownames(uniref_counts_ordered)
#-- Monocle Object
#uniref_counts_ordered <- as.matrix(uniref_counts_ordered)
if(class(uniref_counts_ordered) != 'matrix')
  uniref_counts_ordered <- as.matrix(uniref_counts_ordered)

UniRef_counts_Dataset <- newCellDataSet(as(uniref_counts_ordered, "sparseMatrix"),
                                     phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
UniRef_counts_Dataset <- newCellDataSet(uniref_counts_ordered,
                                        phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
#if(class(UniRef_counts_Dataset) != 'matrix')
 # UniRef_counts_Dataset <- as.matrix(UniRef_counts_Dataset)

#Estimate size factors and dispersions
UniRef_counts_Dataset <- estimateSizeFactors(UniRef_counts_Dataset)
UniRef_counts_Dataset <- estimateDispersions(UniRef_counts_Dataset)


#Filtering low-quality cells
UniRef_counts_Dataset <- detectGenes(UniRef_counts_Dataset, min_expr = 0.1)
print(head(fData(UniRef_counts_Dataset)))
expressed_genes <- row.names(subset(fData(UniRef_counts_Dataset),
                                    num_cells_expressed >= 5)) #The vector expressed_genes now holds the identifiers for taxa Counts in at least 5 Samples of the data set.


#-- Phenodata with number of genes expressed and Size factor
print(head(pData(UniRef_counts_Dataset)))

##########################
#-- Distribution of Counts
##########################
pData(UniRef_counts_Dataset)$Total_counts <- Matrix::colSums(exprs(UniRef_counts_Dataset))
#Ref_counts_Dataset <- Ref_counts_Dataset[,pData(Ref_counts_Dataset)$Total_counts < 1e6]
upper_bound <- 10^(mean(log10(pData(UniRef_counts_Dataset)$Total_counts)) +
                     2*sd(log10(pData(UniRef_counts_Dataset)$Total_counts)))
lower_bound <- 10^(mean(log10(pData(UniRef_counts_Dataset)$Total_counts)) -
                     2*sd(log10(pData(UniRef_counts_Dataset)$Total_counts)))
qplot(Total_counts, data = pData(UniRef_counts_Dataset), color = diagnosis, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)


##########################
#Clustering cells without marker genes
##########################
disp_table <- dispersionTable(UniRef_counts_Dataset)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
UniRef_counts_Dataset <- setOrderingFilter(UniRef_counts_Dataset, unsup_clustering_genes$gene_id)
plot_ordering_genes(UniRef_counts_Dataset)

plot_pc_variance_explained(UniRef_counts_Dataset, return_all = F) # norm_method='log'


############################################################################################
#--- Trajectory Analysis -- Odering Genes based on Differentially Expressed among Diagnosis
############################################################################################
#Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(UniRef_counts_Dataset[expressed_genes,],
                                      fullModelFormulaStr = "~diagnosis")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
UniRef_counts_Dataset <- setOrderingFilter(UniRef_counts_Dataset, ordering_genes)
plot_ordering_genes(UniRef_counts_Dataset)

#Trajectory step 2: reduce data dimensionality
UniRef_counts_Dataset <- reduceDimension(UniRef_counts_Dataset, max_components = 2,
                                      method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory
library(igraph)
#library(SingleCellExperiment)
source("~/Library/CloudStorage/Box-Box/David_Casero_Lab/HMP_IBD_Project/Count_Trajectory/supervised_monocle_uniRef/ReducedDim.R")
source("~/Library/CloudStorage/Box-Box/David_Casero_Lab/HMP_IBD_Project/Count_Trajectory/supervised_monocle_uniRef/order_cells.R")
UniRef_counts_Dataset <- orderCells(UniRef_counts_Dataset, root_state = NULL, num_paths = NULL, reverse = NULL)

jpeg("Trajectory_selTaxa_diagnosis.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(UniRef_counts_Dataset, color_by = "diagnosis")
dev.off ()
jpeg("Trajectory_selTaxa_States.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(UniRef_counts_Dataset, color_by = "State")
dev.off ()

############################################################################################################
#------------------------------------- Individual Plots for CD, nonIBD and UC ------ **** Not mendatory
#"State" is just Monocle's term for the segment of the tree.
#The function below is handy for identifying the State which contains most of the cells from  *****CD*****
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$diagnosis)[,"CD"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
UniRef_counts_Dataset <- orderCells(UniRef_counts_Dataset, root_state = GM_state(UniRef_counts_Dataset))
jpeg("Trajectory_selUniref_pseudoTimeCD.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(UniRef_counts_Dataset, color_by = "Pseudotime")
dev.off ()

#The function below is handy for identifying the State which contains most of the cells from *****NonIBD*******
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$diagnosis)[,"nonIBD"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
UniRef_counts_Dataset <- orderCells(UniRef_counts_Dataset, root_state = GM_state(UniRef_counts_Dataset))
#The function below is handy for identifying the State which contains most of the cells from CD
jpeg("Trajectory_selUniRef_pseudoTime_nonIBD.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(UniRef_counts_Dataset, color_by = "Pseudotime")
dev.off ()

#The function below is handy for identifying the State which contains most of the cells from *****UC******
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$diagnosis)[,"UC"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
UniRef_counts_Dataset <- orderCells(UniRef_counts_Dataset, root_state = GM_state(UniRef_counts_Dataset))
#The function below is handy for identifying the State which contains most of the cells from CD
jpeg("Trajectory_selUniref_pseudoTime_UC.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(UniRef_counts_Dataset, color_by = "Pseudotime")
dev.off ()
############################################################################################################

jpeg("Trajectory_Stats.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(UniRef_counts_Dataset, color_by = "State") +
  facet_wrap(~State, nrow = 1)
dev.off ()

###################################################################################################
#-----*Few Selected Genes Based on Lowest Q values (Should be preesent in more than 100 Samples)
#--- Plot Selected Genes
###################################################################################################
DE_Genes_sample100 <- subset(diff_test_res, diff_test_res$num_cells_expressed > 100)

sel_genes <- row.names(subset(fData(UniRef_counts_Dataset),
                                gene_short_name %in% c("UniRef90_L1NLB8", "UniRef90_R6WWF8", "UniRef90_A0A1M7KQW2")))
jpeg("State_Genes_expression.jpg", height = 6, width = 5, units = 'in', res = 600)
plot_genes_jitter(UniRef_counts_Dataset[sel_genes,],
                  grouping = "State",
                  min_expr = 0.1)
dev.off ()

Ref_expressed_genes <- row.names(subset(fData(UniRef_counts_Dataset),
                                        num_cells_expressed >= 5))
UniRef_counts_Dataset_filtered <- UniRef_counts_Dataset[Ref_expressed_genes,]
my_genes <- row.names(subset(fData(UniRef_counts_Dataset_filtered),
                             gene_short_name %in% c("UniRef90_L1NLB8", "UniRef90_R6WWF8", "UniRef90_A0A1M7KQW2")))
cds_subset <- UniRef_counts_Dataset_filtered[my_genes,]
jpeg("Pseudotime_Genes_expression.jpg", height = 6, width = 5, units = 'in', res = 600)
plot_genes_in_pseudotime(cds_subset, color_by = "diagnosis")
dev.off ()

############################################################################################
#Alternative choices for ordering genes
#Ordering based on genes that differ between clusters
############################################################################################
#We recommend users a order samples using taxa selected by an unsupervised procedure called "dpFeature".
#To use dpFeature, we first select superset of feature taxa as taxa expressed in at least 5% of all the samples
UniRef_counts_Dataset_AC <- detectGenes(UniRef_counts_Dataset, min_expr = 0.1)
fData(UniRef_counts_Dataset_AC)$use_for_ordering <-
  fData(UniRef_counts_Dataset_AC)$num_cells_expressed > 0.05 * ncol(UniRef_counts_Dataset_AC)
#- PCA Analysis
plot_pc_variance_explained(UniRef_counts_Dataset_AC, return_all = F)
#- Reduce dimension
UniRef_counts_Dataset_AC <- reduceDimension(UniRef_counts_Dataset_AC,
                                      max_components = 2,
                                      norm_method = 'log',
                                      num_dim = 3,
                                      reduction_method = 'tSNE',
                                      verbose = T)

UniRef_counts_Dataset_AC <- clusterCells(UniRef_counts_Dataset_AC, verbose = F)

plot_cell_clusters(UniRef_counts_Dataset_AC, color_by = 'as.factor(Cluster)')
plot_cell_clusters(UniRef_counts_Dataset_AC, color_by = 'as.factor(diagnosis)')



#----- Work with Differentially Expressed Genes Based on State
diff_test_res_State <- differentialGeneTest(UniRef_counts_Dataset[expressed_genes,],
                                      fullModelFormulaStr = "~State")
sig_genes <- subset(diff_test_res_State, qval < 0.01)
sig_genes[,c("gene_short_name", "pval", "qval")]
check_marker <- UniRef_counts_Dataset[row.names(subset(fData(UniRef_counts_Dataset),
                                      gene_short_name %in% c("UniRef90_A0A014AUH4", "UniRef90_A0A015NVK9"))),]
plot_genes_jitter(check_marker, grouping = "State", ncol= 2)


#----- Work with Differentially Expressed Genes Based on PseudoTime
diff_test_res_Pseudotime <- differentialGeneTest(UniRef_counts_Dataset[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_genes_pseudotime <- subset(diff_test_res_Pseudotime, qval < 0.01)
sig_genes_pseudotime[,c("gene_short_name", "pval", "qval")]

check_marker <- UniRef_counts_Dataset[row.names(subset(fData(UniRef_counts_Dataset),
                                                       gene_short_name %in% c("UniRef90_A0A014AUH4", "UniRef90_A0A015NVK9"))),]
plot_genes_in_pseudotime(check_marker, color_by = "diagnosis")

sig_genes_pseudotime_rownames <- row.names(subset(diff_test_res_Pseudotime, qval < 0.01))
plot_pseudotime_heatmap(UniRef_counts_Dataset[sig_genes_pseudotime_rownames,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
