meta <- read.csv (file="meta_pid.csv", row.names = 1)
UniRef_counts <- read.csv(file="Bowtie2_map_uniRefGeneCounts.csv", sep = ",", row.names = 1)
UniRef_counts <- data.frame(t(UniRef_counts))
library (labdsv)
UniRef_counts_meta <- merge(meta, UniRef_counts, by=0, all = F)
rownames(UniRef_counts_meta) <- UniRef_counts_meta$Row.names; UniRef_counts_meta$Row.names <- NULL

#for (i in 1:length(pid)){
#  count <- subset(UniRef_counts_meta, pid == pid[i])
#  count <- count[,3:ncol(count)]
#  Uniref_present <- data.frame(colnames(count_filtered))
#  name <- rep(1,nrow(Uniref_present))
#  Uniref_present_df <- cbind(Uniref_present, name)
#  rownames(Uniref_present_df) <- Uniref_present_df$colnames.count_filtered.; Uniref_present_df$colnames.count_filtered. <- NULL
#  colnames(Uniref_present_df) <- pid[i]
  
#  remainingUniref = count[,!(colnames(count) %in% rownames(Uniref_present_df))]
#  Uniref_absent <- data.frame(colnames(remainingUniref))
#  name<-rep(0,nrow(Uniref_absent)) #-- Assign Zeros to this patient
# Uniref_absent_df <- cbind(Uniref_absent, name)
#  rownames(Uniref_absent_df) <- Uniref_absent_df$colnames.remainingUniref.; Uniref_absent_df$colnames.remainingUniref. <- NULL
# colnames(Uniref_absent_df) <- pid[i]
  
#  Presence_absence <- rbind(Uniref_absent_df, Uniref_present_df)
#  write.csv2(Presence_absence, paste0("~/Box/David_Casero_Lab/HMP_IBD_Project/Count_Trajectory/New_count_Trajectory_Mapped_Counts/Check_outputs/",pid[i],"_presence_absence.csv"),row.names = TRUE)
#}

UniRef_counts_meta_New <- UniRef_counts_meta[,c(1,3:ncol(UniRef_counts_meta))]
#aggregate(cbind(UniRef90_A0A009F602, UniRef90_A0A010Z266) ~ pid, data = UniRef_counts_meta, FUN = mean, na.rm = TRUE)
UniRef_aggregate <- aggregate(x = UniRef_counts_meta_New[ , colnames(UniRef_counts_meta_New) != "pid"],             # Mean by group
          by = list(UniRef_counts_meta_New$pid),
          FUN = sum)
rownames(UniRef_aggregate) <- UniRef_aggregate$Group.1; UniRef_aggregate$Group.1 <- NULL
UniRef_aggregate_CheckCopy <- UniRef_aggregate
UniRef_aggregate[UniRef_aggregate > 0] <- 1 

UniRef_aggregate_Transpose <- data.frame(t(UniRef_aggregate))
UniRef_psesncein_EachPID <- data.frame(rowSums(UniRef_aggregate_Transpose))
write.csv(UniRef_psesncein_EachPID, file = "UniRef_psesncein_EachPID.csv")
