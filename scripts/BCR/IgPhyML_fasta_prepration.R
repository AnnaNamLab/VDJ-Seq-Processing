library(dplyr)

sample="HL10"
DominantChain="IGL"
model_name="GY"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample))
dir.create(paste0('igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined'))
dir.create(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final'))

filtered_BCR_processed = read.csv(paste0('no_doublet_VDJfixed_clean_BCR_data_HRS_Bcells_',DominantChain,'.csv'))
filtered_BCR_processed= filtered_BCR_processed %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count))  


### read singlet

pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/FINAL_doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))
singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
filtered_BCR_processed= filtered_BCR_processed[filtered_BCR_processed$cell_id1 %in% singlet,]

####reading the cell states##
new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_May2025/2025-07-03_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)

##### filtering the cells to keep HLsample HRS and Bcells after removing VDJ and ; 
HLsample_HRS=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("HRS"),] 
HLsample_Bcells=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("Bcells"),] 
HLsample_Tcells=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("Tcells"),] 

write.table(paste0(HLsample_HRS$cell_id1,'-1'),'HLsample_HRS.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)
write.table(paste0(HLsample_Bcells$cell_id1,'-1'),'HLsample_Bcells.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)
write.table(paste0(HLsample_Tcells$cell_id1,'-1'),'HLsample_Tcells.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)

###import the output file of data with fixed VDJ and removed doublets, call it filtered_BCR_processed # this file is free of doublets and fixed VDJ

filtered_BCR_processed=filtered_BCR_processed[filtered_BCR_processed$cell_id1 %in% HLsample_HRS$cell_id1,]# getting only HRS cells
ind=match(filtered_BCR_processed$cell_id1,HLsample_HRS$cell_id1)# adding the cell state
filtered_BCR_processed$new_cell_state=HLsample_HRS$SubtypeName[ind]
copy_filtered_BCR_processed=filtered_BCR_processed


##### keeping only this group of HRS cells ("Cycling",   "Inflammatory",      "Metabolic",    "Bcell-HIGH", "Neural-Glial")
filtered_BCR_processed=filtered_BCR_processed[filtered_BCR_processed$new_cell_state %in% c("MetaboHIGH", "Inflammatory", "Cycling","Neural" ,"BcellHIGH", "UPR_Secretory"),]
filtered_BCR_processed$new_cell_state = gsub("UPR_Secretory","Secretory" , filtered_BCR_processed$new_cell_state)

## reading the airr.tsv file output of TRUST4 to extract the full sequence
TRUST4_airr = read.table(Sys.glob('*_barcode_airr.tsv'), sep='\t', header=TRUE)
filtered_BCR_processed_sequence_ind= match(filtered_BCR_processed$contig_id, TRUST4_airr$sequence_id)
filtered_BCR_processed$sequence = TRUST4_airr$sequence[filtered_BCR_processed_sequence_ind]
filtered_BCR_processed$productive = TRUST4_airr$productive[filtered_BCR_processed_sequence_ind]
filtered_BCR_processed$junction_aa = TRUST4_airr$junction_aa[filtered_BCR_processed_sequence_ind]

#filtered_BCR_processed = filtered_BCR_processed[which(filtered_BCR_processed$productive==TRUE & is.na(filtered_BCR_processed$CDR3)==FALSE), ]
filtered_BCR_processed = filtered_BCR_processed[which(is.na(filtered_BCR_processed$CDR3)==FALSE), ]


# Step 1: Create initial status
filtered_BCR_processed$cdr3_status <- ifelse(
  is.na(filtered_BCR_processed$junction_aa) | filtered_BCR_processed$junction_aa == "",
  "missing",
  ifelse(
    substr(filtered_BCR_processed$junction_aa, 1, 1) != "C" &
      substr(filtered_BCR_processed$junction_aa, nchar(filtered_BCR_processed$junction_aa),
             nchar(filtered_BCR_processed$junction_aa)) %in% c("W", "F"),
    "invalid_start",
    ifelse(
      substr(filtered_BCR_processed$junction_aa, 1, 1) == "C" &
        !(substr(filtered_BCR_processed$junction_aa, nchar(filtered_BCR_processed$junction_aa),
                 nchar(filtered_BCR_processed$junction_aa)) %in% c("W", "F")),
      "invalid_stop",
      ifelse(
        substr(filtered_BCR_processed$junction_aa, 1, 1) == "C" &
          substr(filtered_BCR_processed$junction_aa, nchar(filtered_BCR_processed$junction_aa),
                 nchar(filtered_BCR_processed$junction_aa)) %in% c("W", "F"),
        "start and stop valid",
        "invalid"
      )
    )
  )
)


# Step 2: Append CDR3 similarity for problematic sequences
problematic <- c("missing", "invalid_start", "invalid_stop", "invalid")

filtered_BCR_processed$cdr3_status <- ifelse(
  filtered_BCR_processed$cdr3_status %in% problematic & 
    !is.na(filtered_BCR_processed$CDR3_germline_similarity) &
    filtered_BCR_processed$CDR3_germline_similarity >= 0.7,
  paste0(filtered_BCR_processed$cdr3_status, "_CDR3_similar_germline"),
  filtered_BCR_processed$cdr3_status
)

# Check results
table(filtered_BCR_processed$cdr3_status)

# Define valid start and end codons
valid_start <- c("TGT", "TGC")       # Start C codons
valid_end   <- c("TTT", "TTC", "TGG") # End F/W codons
# Function to check and filter sequences
df_filtered <- filtered_BCR_processed %>%
  rowwise() %>%
  filter(
    substr(CDR3, 1, 3) %in% valid_start &
      substr(CDR3, nchar(CDR3)-2, nchar(CDR3)) %in% valid_end
  ) %>%
  ungroup()


# Create a bar plot of CDR3 status
pdf(paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/",sample,"_CDR3_stop_end_validation.pdf"))
ggplot(filtered_BCR_processed, aes(x = cdr3_status)) +
  geom_bar(fill = "steelblue4") +
  theme_minimal() +
  labs(
    title = "Distribution of CDR3 Status",
    x = "CDR3 Status",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggplot(df_filtered, aes(x = cdr3_status)) +
  geom_bar(fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Distribution of CDR3 Status",
    x = "CDR3 Status",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

dev.off()
write.table(filtered_BCR_processed,  paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/filtered_BCR_processed_CDR3_stop_end_validationInvalidIncluded.csv"))

invalid_contig= setdiff(filtered_BCR_processed$contig_id , df_filtered$contig_id)
write.table(filtered_BCR_processed[filtered_BCR_processed$contig_id %in% invalid_contig,],paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/InvalidContigs.csv") )

filtered_BCR_processed = df_filtered
#generating list of seqs based on unique CDR3
seqs <- readDNAStringSet(paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/MAFFTaligned_StopCodonNNN_CDR3_sequences_",sample,"_new_corrected_validCodon_fixed.fasta"))

df <- data.frame(
  header = names(seqs),
  sequence = as.character(seqs),
  stringsAsFactors = FALSE
)
split_list <- split(df, df$sequence)

sampled_list <- lapply(split_list, function(x) {
  x[sample(nrow(x), 1), , drop = FALSE]
})



sampled_df <- do.call(rbind, sampled_list)

# Convert back to DNAStringSet
sampled_seqs <- DNAStringSet(sampled_df$sequence)
names(sampled_seqs) <- sampled_df$header

write.csv(filtered_BCR_processed , 'filtered_BCR_processed.csv')
# Write to new FASTA
writeXStringSet(sampled_seqs, paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/sampled_unique_sequences_",sample,".fasta"))



