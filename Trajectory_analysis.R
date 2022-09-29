#--> Trajectory analyis on Counts
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("monocle")

library (monocle)

#do not run ---> As per the tutorial we need these three files
# http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
#HSMM_expr_matrix <- read.table("fpkm_matrix.txt")
#HSMM_sample_sheet <- read.delim("cell_sample_sheet.txt")
#HSMM_gene_annotation <- read.delim("gene_annotations.txt")

#- Get metadata
PTR_meta <- read.csv (file="~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.txt", sep="\t", row.names =1, header = T)
sample_sheet <- PTR_meta[,71:90]
pd <- new("AnnotatedDataFrame", data = sample_sheet)

#- Get Ref Gene Counts
Ref_counts <- read.csv(file="RefGenome_Counts_output.csv", row.names = 1, header = T, sep = ",")
# Fix names
rownames(Ref_counts) <- gsub(".hgm", "", rownames(Ref_counts), fixed=TRUE)
rownames(Ref_counts) <- gsub(".patric", "", rownames(Ref_counts), fixed=TRUE)
rownames(Ref_counts) <- gsub(".img", "", rownames(Ref_counts), fixed=TRUE)
rownames(Ref_counts) <- gsub("X", "", rownames(Ref_counts), fixed=TRUE)
samples <- row.names(PTR_meta)
Ref_counts_Sel <- Ref_counts[,samples] #-- to keep same samples as we have in the metadata file

#-- Check Sparsity(i.e. the percentage of 0 values) of a Dataframe
sum(Ref_counts_Sel == 0)/(dim(Ref_counts_Sel)[1]*dim(Ref_counts_Sel)[2])
#[1] 0.9082249: 90%

#- Get gene name File (here it is taxa names)
Refgid_gname_Org <- read.csv("../PTR_analysis_on_IGG_Ref/RefGid_Gnames.txt", sep = "\t", header = T, row.names = 1)
Refgid_gname_Ref_counts <- merge(Refgid_gname_Org, Ref_counts_Sel, by=0, all=F)
rownames(Refgid_gname_Ref_counts) <- Refgid_gname_Ref_counts$Row.names; Refgid_gname_Ref_counts$Row.names <- NULL
Refgid_gname <- Refgid_gname_Ref_counts[,1] #-- Just to get same names as we have in Expression matrix
Refgid_gname <- data.frame(Refgid_gname)
row.names(Refgid_gname) <- row.names(Refgid_gname_Ref_counts)
colnames(Refgid_gname) <- "gene_short_name"
#write.table(Refgid_gname, file="Refgid_gname.txt", sep="\t")
#_ These steps to reanme duplicate Taxa Names
Refgid_gname_uniq <- make.unique(Refgid_gname$gene_short_name)
Refgid_gname_uniq <- data.frame(Refgid_gname_uniq)
Refgid_gname_updated <- cbind(Refgid_gname, Refgid_gname_uniq)
Refgid_gname_updated <- data.frame(Refgid_gname_updated[,2])
row.names(Refgid_gname_updated) <- row.names(Refgid_gname)
colnames(Refgid_gname_updated) <- "gene_short_name"

fd <- new("AnnotatedDataFrame", data = Refgid_gname_updated)

#To work with count data, specify the negative binomial distribution as the expressionFamily argument to newCellDataSet:
Ref_counts_Dataset <- newCellDataSet(as.matrix(Ref_counts_Sel),
                       phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

#Estimate size factors and dispersions
Ref_counts_Dataset <- estimateSizeFactors(Ref_counts_Dataset)
Ref_counts_Dataset <- estimateDispersions(Ref_counts_Dataset)

#Filtering low-quality cells
Ref_counts_Dataset <- detectGenes(Ref_counts_Dataset, min_expr = 0.1)
print(head(fData(Ref_counts_Dataset)))
expressed_genes <- row.names(subset(fData(Ref_counts_Dataset),
                                    num_cells_expressed >= 5)) #The vector expressed_genes now holds the identifiers for taxa Counts in at least 5 Samples of the data set.

#-- Phenodata with number of genes expressed and Size factor
print(head(pData(Ref_counts_Dataset)))

##########################
#-- Distribution of Counts
##########################
pData(Ref_counts_Dataset)$Total_counts <- Matrix::colSums(exprs(Ref_counts_Dataset))
#Ref_counts_Dataset <- Ref_counts_Dataset[,pData(Ref_counts_Dataset)$Total_counts < 1e6]
upper_bound <- 10^(mean(log10(pData(Ref_counts_Dataset)$Total_counts)) +
                     2*sd(log10(pData(Ref_counts_Dataset)$Total_counts)))
lower_bound <- 10^(mean(log10(pData(Ref_counts_Dataset)$Total_counts)) -
                     2*sd(log10(pData(Ref_counts_Dataset)$Total_counts)))
qplot(Total_counts, data = pData(Ref_counts_Dataset), color = diagnosis, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

##########################
#-- Log Transformation
##########################
# Log-transform each value in the expression matrix.
#L <- log(exprs(Ref_counts_Dataset[expressed_genes,]))
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
#melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of the standardized gene expression values.
#qplot(value, geom = "density", data = melted_dens_df) +
#  stat_function(fun = dnorm, size = 0.5, color = 'red') +
#  xlab("Standardized log(FPKM)") +
#  ylab("Density")

##########################
#Clustering cells without marker genes
##########################
disp_table <- dispersionTable(Ref_counts_Dataset)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Ref_counts_Dataset <- setOrderingFilter(Ref_counts_Dataset, unsup_clustering_genes$gene_id)
plot_ordering_genes(Ref_counts_Dataset)

plot_pc_variance_explained(Ref_counts_Dataset, return_all = F) # norm_method='log'

############################################################################################
#--- Trajectory Analysis -- Odering Genes based on Differentially Expressed among Diagnosis
############################################################################################
#Trajectory step 1: choose genes that define a cell's progress
diff_test_res <- differentialGeneTest(Ref_counts_Dataset[expressed_genes,],
                                      fullModelFormulaStr = "~diagnosis")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
Ref_counts_Dataset <- setOrderingFilter(Ref_counts_Dataset, ordering_genes)
plot_ordering_genes(Ref_counts_Dataset)

#Trajectory step 2: reduce data dimensionality
Ref_counts_Dataset <- reduceDimension(Ref_counts_Dataset, max_components = 2,
                            method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory
Ref_counts_Dataset <- orderCells(Ref_counts_Dataset)
jpeg("Trajectory_selTaxa_diagnosis.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(Ref_counts_Dataset, color_by = "diagnosis")
dev.off ()
jpeg("Trajectory_selTaxa_States.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(Ref_counts_Dataset, color_by = "State")
dev.off ()

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
Ref_counts_Dataset <- orderCells(Ref_counts_Dataset, root_state = GM_state(Ref_counts_Dataset))
jpeg("Trajectory_selTaxa_pseudoTimeCD.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(Ref_counts_Dataset, color_by = "Pseudotime")
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
Ref_counts_Dataset <- orderCells(Ref_counts_Dataset, root_state = GM_state(Ref_counts_Dataset))
#The function below is handy for identifying the State which contains most of the cells from CD
jpeg("Trajectory_selTaxa_pseudoTime_nonIBD.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(Ref_counts_Dataset, color_by = "Pseudotime")
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
Ref_counts_Dataset <- orderCells(Ref_counts_Dataset, root_state = GM_state(Ref_counts_Dataset))
#The function below is handy for identifying the State which contains most of the cells from CD
jpeg("Trajectory_selTaxa_pseudoTime_UC.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(Ref_counts_Dataset, color_by = "Pseudotime")
dev.off ()

jpeg("Trajectory_Stats.jpg", height = 5, width = 5, units = 'in', res = 600)
plot_cell_trajectory(Ref_counts_Dataset, color_by = "State") +
  facet_wrap(~State, nrow = 1)
dev.off ()

####################
#--- Plot few Taxa
###################
blast_genes <- row.names(subset(fData(Ref_counts_Dataset),
                                gene_short_name %in% c("Roseburia intestinalis", "Subdoligranulum sp.", "Flavonifractor plautii")))
jpeg("State_Taxa_expression.jpg", height = 6, width = 5, units = 'in', res = 600)
plot_genes_jitter(Ref_counts_Dataset[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
dev.off ()

Ref_expressed_genes <- row.names(subset(fData(Ref_counts_Dataset),
                                    num_cells_expressed >= 5))
Ref_counts_Dataset_filtered <- Ref_counts_Dataset[Ref_expressed_genes,]
my_genes <- row.names(subset(fData(Ref_counts_Dataset_filtered),
                             gene_short_name %in% c("Roseburia intestinalis", "Subdoligranulum sp.", "Flavonifractor plautii")))
cds_subset <- Ref_counts_Dataset_filtered[my_genes,]
jpeg("Pseudotime_Taxa_expression.jpg", height = 6, width = 5, units = 'in', res = 600)
plot_genes_in_pseudotime(cds_subset, color_by = "diagnosis")
dev.off ()

############################################################################################
#Alternative choices for ordering genes
#Ordering based on genes that differ between clusters
############################################################################################
#We recommend users a order samples using taxa selected by an unsupervised procedure called "dpFeature".
#To use dpFeature, we first select superset of feature taxa as taxa expressed in at least 5% of all the samples
Ref_counts_Dataset <- detectGenes(Ref_counts_Dataset, min_expr = 0.1)
fData(Ref_counts_Dataset)$use_for_ordering <-
  fData(Ref_counts_Dataset)$num_cells_expressed > 0.05 * ncol(Ref_counts_Dataset)
#- PCA Analysis
plot_pc_variance_explained(Ref_counts_Dataset, return_all = F)
#- Reduce dimension
Ref_counts_Dataset <- reduceDimension(Ref_counts_Dataset,
                            max_components = 2,
                            norm_method = 'log',
                            num_dim = 3,
                            reduction_method = 'tSNE',
                            verbose = T, check_duplicates = FALSE)

Ref_counts_Dataset <- clusterCells(Ref_counts_Dataset, verbose = F)

jpeg("PlotClusters.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, color_by = 'as.factor(Cluster)')
dev.off ()
jpeg("PlotClusters_diagnosis.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, color_by = 'as.factor(diagnosis)')
dev.off ()
jpeg("PlotClusters_Visits.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, 1, 2, color = "visit")
dev.off ()
jpeg("PlotClusters_BowelSurgery.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, 1, 2, color = "bowel_surgery")
dev.off ()
jpeg("PlotCluster_Taxa_expression.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, 1, 2, color = "diagnosis",
                   markers = c("Roseburia intestinalis", "Subdoligranulum sp.", "Flavonifractor plautii"))
dev.off ()

#--- Check the Thresolds for Clustering and then Re-run the Analysis
plot_rho_delta(Ref_counts_Dataset, rho_threshold = 2, delta_threshold = 4 )

#---- Clustering Again
Ref_counts_Dataset <- clusterCells(Ref_counts_Dataset,
                         rho_threshold = 15,
                         delta_threshold = 10,
                         skip_rho_sigma = T,
                         verbose = F)
plot_cell_clusters(Ref_counts_Dataset, color_by = 'as.factor(Cluster)')
plot_cell_clusters(Ref_counts_Dataset, color_by = 'as.factor(diagnosis)')

jpeg("PlotClusters.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, color_by = 'as.factor(Cluster)')
dev.off ()
jpeg("PlotClusters_diagnosis.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, color_by = 'as.factor(diagnosis)')
dev.off ()
jpeg("PlotClusters_Visits.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, 1, 2, color = "visit")
dev.off ()
jpeg("PlotClusters_BowelSurgery.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, 1, 2, color = "bowel_surgery")
dev.off ()
jpeg("PlotCluster_Taxa_expression.jpg", height = 6, width = 7, units = 'in', res = 600)
plot_cell_clusters(Ref_counts_Dataset, 1, 2, color = "diagnosis",
                   markers = c("Roseburia intestinalis", "Subdoligranulum sp.", "Flavonifractor plautii"))
dev.off ()


#After clustering makes sense, 
#we can then perform differential gene expression test as a way to extract the genes that distinguish them.

clustering_DEG_genes <-
  differentialGeneTest(Ref_counts_Dataset[expressed_genes,],
                       fullModelFormulaStr = '~Cluster',
                       cores = 1)
#We will then select the top 1000 significant genes as the ordering genes.

Ref_counts_Dataset_ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:500]
Ref_counts_Dataset <-
  setOrderingFilter(Ref_counts_Dataset,
                    ordering_genes = Ref_counts_Dataset_ordering_genes)
Ref_counts_Dataset <- reduceDimension(Ref_counts_Dataset, method = 'DDRTree')
Ref_counts_Dataset <- orderCells(Ref_counts_Dataset)
#-- For CD
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$diagnosis)[,"CD"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
Ref_counts_Dataset <- orderCells(Ref_counts_Dataset, root_state = GM_state(Ref_counts_Dataset))

plot_cell_trajectory(Ref_counts_Dataset, color_by = "diagnosis")


