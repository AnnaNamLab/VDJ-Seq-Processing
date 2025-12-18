
# this step is for prepraring the required fasta file for running IgPhyML

library(dplyr)
library(ggplot2)
library(Biostrings)

sample="sample1"     # sample name
DominantChain="IGL"  # clone dominant chain
model_name="GY"      # model used for generating tree in IgPhyMl

setwd(paste0('set_working_directory/',sample))      
dir.create(paste0('igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined')) # output fodler
dir.create(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final')) # outputfolder

filtered_BCR_processed = read.csv(paste0('no_doublet_VDJfixed_clean_BCR_data_HRS_Bcells_',DominantChain,'.csv'))  # this is the output of expandedClone_refinement.R
filtered_BCR_processed= filtered_BCR_processed %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count))  # selecting the contig with highest number of reads


## read singlet

pilot_cellStatus_file= read.csv(Sys.glob(paste0('./doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))
singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
filtered_BCR_processed= filtered_BCR_processed[filtered_BCR_processed$cell_id1 %in% singlet,]

## reading the cell states
new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_May2025/2025-07-03_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)

## filtering the cells to keep HLsample HRS and Bcells after removing VDJ  

HLsample_HRS=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("HRS"),] 
HLsample_Bcells=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("Bcells"),] 
HLsample_Tcells=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("Tcells"),] 

write.table(paste0(HLsample_HRS$cell_id1,'-1'),'HLsample_HRS.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)
write.table(paste0(HLsample_Bcells$cell_id1,'-1'),'HLsample_Bcells.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)
write.table(paste0(HLsample_Tcells$cell_id1,'-1'),'HLsample_Tcells.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)

## import the output file of data with fixed VDJ and removed doublets, call it filtered_BCR_processed # this file is free of doublets and fixed VDJ

filtered_BCR_processed=filtered_BCR_processed[filtered_BCR_processed$cell_id1 %in% HLsample_HRS$cell_id1,]# getting only HRS cells
ind=match(filtered_BCR_processed$cell_id1,HLsample_HRS$cell_id1)# adding the cell state
filtered_BCR_processed$new_cell_state=HLsample_HRS$SubtypeName[ind]
copy_filtered_BCR_processed=filtered_BCR_processed


## keeping only this group of HRS cells ("Cycling",   "Inflammatory",      "Metabolic",    "Bcell-HIGH", "Neural-Glial")

filtered_BCR_processed=filtered_BCR_processed[filtered_BCR_processed$new_cell_state %in% c("MetaboHIGH", "Inflammatory", "Cycling","Neural" ,"BcellHIGH", "UPR_Secretory"),]
filtered_BCR_processed$new_cell_state = gsub("UPR_Secretory","Secretory" , filtered_BCR_processed$new_cell_state)

## reading the airr.tsv file output of TRUST4 to extract the full sequence

TRUST4_airr = read.table(Sys.glob('*_barcode_airr.tsv'), sep='\t', header=TRUE)
filtered_BCR_processed_sequence_ind= match(filtered_BCR_processed$contig_id, TRUST4_airr$sequence_id)
filtered_BCR_processed$sequence = TRUST4_airr$sequence[filtered_BCR_processed_sequence_ind]
filtered_BCR_processed$productive = TRUST4_airr$productive[filtered_BCR_processed_sequence_ind]
filtered_BCR_processed$junction_aa = TRUST4_airr$junction_aa[filtered_BCR_processed_sequence_ind]
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


## if there is stop codon in the seq replaced with NNN
mask_stops_if_present <- function(seq) {
  # Convert to character to allow manipulation
  s <- as.character(seq)
  len <- nchar(s)
  
  # Prepare output sequence as character vector
  seq_chars <- strsplit(s, "")[[1]]
  
  # Only replace codons that start at 1, 4, 7, ...
  for (i in seq(1, len - 2, by = 3)) {
    codon <- substr(s, i, i + 2)
    
    # Skip codons that contain gaps (-)
    if (grepl("-", codon)) next
    
    # Replace stop codons with NNN
    if (codon %in% c("TAA", "TAG", "TGA")) {
      seq_chars[i:(i+2)] <- "N"
    }
  }
  
  # Return a DNAString with gaps preserved
  return(DNAString(paste(seq_chars, collapse = "")))
}


fasta_seqs <- DNAStringSet(filtered_BCR_processed$CDR3)
names(fasta_seqs) <- paste0(filtered_BCR_processed$contig_id, " ",filtered_BCR_processed$cdr3_status)

# Apply the function to all sequences

writeXStringSet(fasta_seqs, paste0('/Users/saramoein/Documents/new_run_HL_May2025/', sample, '/igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/raw_fasta.fasta'))

# -------------------------------
# Save masked sequences
# -------------------------------


Sys.setenv(PATH = paste("/opt/homebrew/bin", Sys.getenv("PATH"), sep = ":"))
system("which mafft")

fasta_mafft_input <- paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/raw_fasta.fasta')
fasta_mafft_output <- paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/mafft_fasta_aligned.fasta')


fasta_mafft <- readDNAStringSet(fasta_mafft_input)
fasta_mafft_clean <- DNAStringSet(gsub("\\s+", "", toupper(as.character(fasta_mafft))))
writeXStringSet(fasta_mafft_clean, fasta_mafft_input, format = "fasta")
system(paste("mafft --auto", fasta_mafft_input, ">", fasta_mafft_output))

full_aligned_sequences <- readDNAStringSet(fasta_mafft_output)
full_aligned_sequences <- DNAStringSet(toupper(as.character(full_aligned_sequences)))
fasta_mafft_output_uppercase <- paste0(dirname(fasta_mafft_output), "/mafft_fasta_aligned_uppercase.fasta")

full_aligned_sequences <-  DNAStringSet(lapply(full_aligned_sequences, mask_stops_if_present))

writeXStringSet(full_aligned_sequences, fasta_mafft_output_uppercase)

pad_to_multiple_of_3 <- function(seq) {
  s <- as.character(seq)
  remainder <- nchar(s) %% 3
  if (remainder != 0) {
    pad_len <- 3 - remainder
    s <- paste0(s, paste(rep("-", pad_len), collapse = ""))
  }
  return(s)
}

full_aligned_sequences_fixed <- DNAStringSet(sapply(full_aligned_sequences, pad_to_multiple_of_3))
names(full_aligned_sequences_fixed) <- names(full_aligned_sequences)
# Save the aligned sequences
writeXStringSet(full_aligned_sequences_fixed, paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/semiRefined_ExpandedCloneChange/MAFFTaligned_StopCodonNNN_CDR3_sequences_",sample,"_new_corrected_validCodon.fasta"))
BrowseSeqs(full_aligned_sequences_fixed , htmlFile = paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/semiRefined_ExpandedCloneChange/MAFFTaligned_StopCodonNNN_CDR3_sequences_",sample,"_new_corrected_validCodon.html"))


#generating list of seqs based on unique CDR3
seqs <- full_aligned_sequences_fixed  

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
sampled_df <- do.call(rbind, sampled_list)


write.csv(filtered_BCR_processed , 'filtered_BCR_processed.csv')

# Write to new FASTA
writeXStringSet(sampled_seqs, paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/sampled_unique_sequences_",sample,".fasta"))



