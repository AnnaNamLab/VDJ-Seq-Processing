library(Biostrings)
library(stringr)
library(dplyr)
library(colorRamp2)
library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(ape)
library(ggalt)
library(ggmsa)
library(phytools)


alignment_score <- function(CDR3_data, alignment_file) {
  # -----------------------------
  # CDR3_data: vector of original sequences (for reference only)
  # alignment_file: MAFFT-aligned FASTA file
  # contig_ids: vector of contig IDs matching CDR3_data
  # -----------------------------
  
  #  Read aligned sequences
  contig_ids= data.frame(CDR3_data)$contig_id
  
  contig_ids <- sub("startandstopvalid", " start and stop valid", contig_ids)
  contig_ids <- sub("missing", " missing", contig_ids)
  contig_ids <- sub("invalid", " invalid", contig_ids)
  aligned <- readDNAStringSet(alignment_file)
  aligned_names <- names(aligned)
  #aligned_contigs <- sub(" .*", "", aligned_names)  # remove trailing info if present
  
  #  Reorder aligned sequences to match input contig_ids
  
  aligned <- aligned[match(contig_ids, aligned_names)]
  
  # check contigs
  if (any(is.na(aligned))) {
    stop("Some contig_ids not found in alignment file.")
  }
  
  #  Initialize results
  n <- nrow(CDR3_data)
  dist1_matrix <- matrix(NA, n, n)
  main_string_score_align <- numeric(n)
  
  #  Loop over all sequences
  for (i in seq_len(n)) {
    print(i)
    main_string <- aligned[i]
    dist1 <- numeric(n)
    
    for (j in seq_len(n)) {
      
      # Compute Levenshtein (edit) distance between aligned sequences
      dist1[j] <- as.numeric(Biostrings::stringDist(DNAStringSet(c(main_string, aligned[j])), method = "levenshtein"))
    }
    
    # Mean distance for sequence i
    main_string_score_align[i] <- mean(dist1, na.rm = TRUE)
    dist1_matrix[i, ] <- dist1
  }
  
  # Return results (same structure as your original)
  res_alignment_score <- list(
    scores = main_string_score_align,
    dist1_matrix = dist1_matrix
  )
  
  names(res_alignment_score) <- c("scores", "dist1_matrix")
  
  return(res_alignment_score)
}





sample="HL10"
DominantChain="IGL"
model_name="GY"
working_dir= '/Users/saramoein/Documents/new_run_HL_May2025/'
setwd(paste0(working_dir,sample))
dir.create(paste0('igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined'))
dir.create(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final'))


####reading the cell states##
new_clusters= read.csv(paste0(working_dir, '/2025-07-03_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv'))
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)


### read singlet

pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/FINAL_doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))
singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])


filtered_BCR_processed = read.csv(paste0('filtered_BCR_processed.csv'))

########## reading the tree from IgPhyML

tree= read.tree(paste0("/Users/saramoein/Documents/new_run_HL_May2025/",sample,'/igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/sampled_unique_sequences_',sample,'.fasta_igphyml_tree_singleCDR3_',sample,'.txt'))
tree <- ladderize(tree)
##########NEED TO EXPAND THE TREE

filtered_BCR_processed$status_singlet=NA
filtered_BCR_processed$status_singlet[filtered_BCR_processed$cell_id1 %in% singlet]="singlet"
filtered_BCR_processed$new_VDJ= paste0(filtered_BCR_processed$corrected_Vgene," ",filtered_BCR_processed$corrected_Dgene," ",filtered_BCR_processed$corrected_Jgene)

summarized_filtered_BCR_processed_filtered= filtered_BCR_processed  %>% summarise(n = n(), contig_id,new_VDJ,malig_status,sequence,new_cell_state,CDR3,status_singlet, cdr3_status)



### adding the TipLabels.filtered to tree data
summarized_filtered_BCR_processed_filtered <- as.data.frame(summarized_filtered_BCR_processed_filtered)
summarized_filtered_BCR_processed_filtered$contig_id= paste0(summarized_filtered_BCR_processed_filtered$contig_id,summarized_filtered_BCR_processed_filtered$cdr3_status)
summarized_filtered_BCR_processed_filtered$contig_id = gsub("start and stop valid", "startandstopvalid", summarized_filtered_BCR_processed_filtered$contig_id)
#colnames(summarized_filtered_BCR_processed_filtered)[which(colnames(summarized_filtered_BCR_processed_filtered )=="contig_id")] = "full_contig_id"
#summarized_filtered_BCR_processed_filtered$contig_id= sub(" .*", "", summarized_filtered_BCR_processed_filtered$full_contig_id)

####################### EXPANDING THE TREE ##############
########################################################

# ----------------------------
# 1. Ensure tree has branch lengths
# ----------------------------
if(is.null(tree$edge.length)) tree$edge.length <- rep(1e-6, nrow(tree$edge))
tree$edge.length[is.na(tree$edge.length)] <- 1e-6

# ----------------------------
# 2. Prepare mapping table
# ----------------------------
# Table: summarized_filtered_BCR_processed_filtered
# Columns: contig_id, CDR3
contig_to_cdr3 <- summarized_filtered_BCR_processed_filtered %>%
  select(contig_id, CDR3)

# CDR3 -> all contigs
cdr3_map <- split(contig_to_cdr3$contig_id, contig_to_cdr3$CDR3)

# ----------------------------
# 3. Expand tree by adding duplicates as star clades
# ----------------------------



expanded_tree_original <- tree  # start from your original tree
tiny_length <- 1e-6             # very small branch length

for (tip in tree$tip.label) {
  
  # Get CDR3 for this tip
  cdr3seq <- contig_to_cdr3$CDR3[match(tip, contig_to_cdr3$contig_id)]
  if (is.na(cdr3seq)) next
  
  # All contigs with the same CDR3
  all_contigs <- cdr3_map[[cdr3seq]]
  
  # Skip if only one contig (already in tree)
  if (length(all_contigs) <= 1) next
  
  # New tips to add (exclude the tip already in tree)
  new_tips <- setdiff(all_contigs, tip)
  if (length(new_tips) == 0) next
  
  # Find tip index in current tree
  tip_idx <- match(tip, expanded_tree_original$tip.label)
  if (is.na(tip_idx)) next
  
  # Find the branch length of this tip
  branch_len <- expanded_tree_original$edge.length[
    which(expanded_tree_original$edge[,2] == tip_idx)
  ]
  
  # Safe tiny position: cannot exceed current branch
  safe_pos <- min(tiny_length, branch_len / 2)
  
  # Add each new tip as a sister
  for (nt in new_tips) {
    single_tip <- stree(1, tip.label = nt)
    single_tip$edge.length <- safe_pos
    
    expanded_tree_original <- ape::bind.tree(
      expanded_tree_original,
      single_tip,
      where = tip_idx,
      position = safe_pos
    )
  }
}


##############################REMOVING DUPLICATE TIPS

# Step 1: Identify duplicates
tip_labels <- expanded_tree_original$tip.label
dup_flags <- duplicated(tip_labels) | duplicated(tip_labels, fromLast = TRUE)

# Get the list of duplicate tip names
dup_names <- unique(tip_labels[dup_flags])

# Step 2: Process each duplicate group
for(tip_name in dup_names){
  
  # Find all positions of this tip
  positions <- which(expanded_tree_original$tip.label == tip_name)
  
  # If more than 1, keep the first, modify it temporarily
  if(length(positions) > 1){
    
    first_pos <- positions[1]
    
    # Temporarily rename the first occurrence to avoid accidental deletion
    expanded_tree_original$tip.label[first_pos] <- paste0(expanded_tree_original$tip.label[first_pos], "_KEEP")
    
    # Remove all other duplicates
    to_drop <- expanded_tree_original$tip.label[positions[-1]]
    expanded_tree_original <- drop.tip(expanded_tree_original, to_drop)
    
    # Restore the first occurrence label back to original
    expanded_tree_original$tip.label[expanded_tree_original$tip.label == paste0(tip_name, "_KEEP")] <- tip_name
    
    
  }
}

tips_to_drop <- which(is.na(expanded_tree_original$tip.label))
expanded_tree_original <- ape::drop.tip(expanded_tree_original, tips_to_drop)




# ----------------------------
# 4. Plot with ggtree
# ----------------------------
ig_expanded <- ggtree(expanded_tree_original , branch.length='branch.length',color="black", size=0.5, linetype=1, ladderize=TRUE) 
#geom_tippoint(shape=21, size=2, fill="white", color="black")

#####################################

summarized_filtered_BCR_processed_filtered = summarized_filtered_BCR_processed_filtered

TipLabels.filtered <- data.frame(tree.label = ig_expanded$data$label, contig_id = ig_expanded$data$label)
summarized_filtered_BCR_processed_filtered <- full_join(x=summarized_filtered_BCR_processed_filtered, y=TipLabels.filtered, by= "contig_id")


### this line probably will be removed
summarized_filtered_BCR_processed_filtered <- summarized_filtered_BCR_processed_filtered[is.na(summarized_filtered_BCR_processed_filtered$contig_id)==FALSE,]
#summarized_filtered_BCR_processed_filtered <- summarized_filtered_BCR_processed_filtered[is.na(summarized_filtered_BCR_processed_filtered$tree.label)==FALSE,]
summarized_filtered_BCR_processed_filtered <- summarized_filtered_BCR_processed_filtered %>%
  distinct(contig_id, .keep_all = TRUE)


rownames(summarized_filtered_BCR_processed_filtered) <- summarized_filtered_BCR_processed_filtered$contig_id




TipLabels.filtered.WithoutNAs <- TipLabels.filtered[is.na(TipLabels.filtered$tree.label)==FALSE ,]### this line probably will be removed
TipLabels.filtered.WithoutNAs <- TipLabels.filtered.WithoutNAs %>%
  distinct(contig_id, .keep_all = TRUE)

rownames(TipLabels.filtered.WithoutNAs) <- TipLabels.filtered.WithoutNAs$tree.label
summarized_filtered_BCR_processed_filtered <-summarized_filtered_BCR_processed_filtered[rownames(TipLabels.filtered.WithoutNAs),]

ig_expanded$data$label <- as.character(ig_expanded$data$label)
summarized_filtered_BCR_processed_filtered$contig_id <- as.character(summarized_filtered_BCR_processed_filtered$contig_id)

clusters.nodeNetwork.filtered =summarized_filtered_BCR_processed_filtered#   filtered_BCR_processed
## creating an empty dataframe and add tree data
emptyRows <- data.frame(matrix(ncol = ncol(clusters.nodeNetwork.filtered), nrow = (nrow(TipLabels.filtered) - nrow(clusters.nodeNetwork.filtered))))

colnames(emptyRows) <- colnames(clusters.nodeNetwork.filtered)
clusters.nodeNetwork.filtered <- rbind(clusters.nodeNetwork.filtered, emptyRows)


### assigning the tips to tree
TipLabels.filtered$status <- clusters.nodeNetwork.filtered$malig_status
TipLabels.filtered$id <-TipLabels.filtered$contig_id
TipLabels.filtered$vdj_status <- clusters.nodeNetwork.filtered$new_VDJ
TipLabels.filtered$new_cell_state <- clusters.nodeNetwork.filtered$new_cell_state
TipLabels.filtered$cdr3_status <- clusters.nodeNetwork.filtered$cdr3_status
TipLabels.filtered$NodeSize <- clusters.nodeNetwork.filtered$n
TipLabels.filtered$CDR3 <- clusters.nodeNetwork.filtered$CDR3
TipLabels.filtered$status_singlet <- clusters.nodeNetwork.filtered$status_singlet
ig_expanded$data$cdr3_status <- as.factor(TipLabels.filtered$cdr3_status)# malignancy status
ig_expanded$data$vdj_status <- as.factor(TipLabels.filtered$vdj_status) # VDJ states 
ig_expanded$data$new_cell_state <- as.factor(TipLabels.filtered$new_cell_state) # cell states
ig_expanded$data$NodeSize <- as.numeric(as.character(TipLabels.filtered$NodeSize)) # The number of cells related to the node associated with a clade
ig_expanded$data$CDR3 <- (TipLabels.filtered$CDR3) 
ig_expanded$data$status_singlet <- (TipLabels.filtered$status_singlet) 
ig_expanded$data$id<-TipLabels.filtered$id
ig_expanded$data$cdr3_status<-TipLabels.filtered$cdr3_status

####PLOTING THE TREE WITH ASSIGNMENTS

r <- ig_expanded+ geom_text2(aes(label=id), hjust=-.3) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

s <- ig_expanded + geom_tiplab(size=3)  + geom_tippoint(aes(col = vdj_status), size= 5, alpha=.75) # + ggtitle(paste0(Variables$TumorId, " - VDJ status"))
q <-  ig_expanded + geom_text2(aes(label=id), hjust=-.3)
cdr3_tree <- ig_expanded + geom_text2(aes(label=CDR3), hjust=-.3) #geom_tippoint(aes(col = cdr3), size= 3, alpha=.75)

cdr3_condition <- ig_expanded + geom_tiplab(size=3)  + geom_tippoint(aes(col = cdr3_status), size= 5, alpha=.75)

pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/VDJ_ig_',sample,'.pdf'),width = 50 ,height = 740)
plot(s)
dev.off()

pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/cdr3_condition_ig_',sample,'.pdf'),width = 50 ,height = 740)
plot(cdr3_condition)
dev.off()

pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/CDR3_ig_',sample,'.pdf'),width = 50 ,height = 740)
plot(cdr3_tree)
dev.off()
brewer.pal(n=5, name = "Set3")
color.tree <- c("BcellHIGH" = "#8DD3C7", "Secretory" = "#FFFFB3", "Inflammatory" = "#BEBADA", "Neural" = "#FB8072", 
                "MetaboHIGH" = "#80B1D3", "Cycling" = "#FDB462")

c <- ig_expanded + geom_tippoint(aes(col = new_cell_state), size= 40, alpha=.75)+ scale_color_manual(values=color.tree)


pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/cell_state_ig_',sample,'_2.pdf'),width = 50, height = 740 )
plot(c)
dev.off()



pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/label_ig_',sample,'.pdf'),width = 50, height = 740 )
plot(r)
dev.off()

long_clade_data=data.frame(table(summarized_filtered_BCR_processed_filtered$CDR3))

long_clade_data=long_clade_data[order(long_clade_data$Freq, decreasing = TRUE),]
summarized_filtered_BCR_processed_filtered$cluster="NA"


ALL_long_clade_data= long_clade_data[which(long_clade_data$Freq > 10),]

long_clade_contigs= summarized_filtered_BCR_processed_filtered$contig_id[which(summarized_filtered_BCR_processed_filtered$CDR3 %in% ALL_long_clade_data$Var1)]


######REMOVING THE CONTIG_ID ON THE MAIN CLADE THAT ARE NOT THE SAME AS DOMINANT CDR3

seq= readDNAStringSet(paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/MAFFTaligned_StopCodonNNN_CDR3_sequences_",sample,"_new_corrected_validCodon_fixed.fasta"))
fasta_names <- names(seq)
fasta_names= gsub(" missing","missing",fasta_names)
fasta_names= gsub(" start and stop valid","startandstopvalid", fasta_names)
fasta_names= gsub(" invalid","invalid", fasta_names)
names(seq) = fasta_names
tree_ladderized <- ladderize(expanded_tree_original)
ordered_tips <- tree_ladderized$tip.label[tree_ladderized$edge[tree_ladderized$edge[, 2] <= length(tree_ladderized$tip.label), 2]]

seq_ordered <- seq[match(ordered_tips, names(seq))]
writeXStringSet(seq_ordered, paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample,'/igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/msa_seq_validCodon_fixed.fasta'))
BrowseSeqs(seq_ordered,paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample,'/igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/msa_seq_validCodon_fixed.html') )

library(pagedown)

outfile_html <- paste0('/Users/saramoein/Documents/new_run_HL_May2025/',
                       sample,
                       '/igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/msa_seq_validCodon_fixed.html')

outfile_pdf <- sub("\\.html$", ".pdf", outfile_html)

pagedown::chrome_print(input = outfile_html, output = outfile_pdf)

summarized_filtered_BCR_processed_filtered$clone_id <- as.integer(factor(summarized_filtered_BCR_processed_filtered$CDR3))

freq_df = ALL_long_clade_data
freq_df$clone_name <- paste0("Clone_", seq_len(nrow(freq_df)), "_size",freq_df$Freq)
df = summarized_filtered_BCR_processed_filtered
# Merge: left join by CDR3
df2 <- df %>%
  dplyr::left_join(freq_df[, c("Var1", "clone_name")],
                   by = c("CDR3" = "Var1"))

df2$clone_name[is.na(df2$clone_name)] <- "non-Expanded"
summarized_filtered_BCR_processed_filtered = df2

summarized_filtered_BCR_processed_filtered$cluster = ifelse(summarized_filtered_BCR_processed_filtered$CDR3 %in% ALL_long_clade_data$Var1 , "Expanded" , 'non-Expanded')

summarized_filtered_BCR_processed_filtered_exclude_mainClade= summarized_filtered_BCR_processed_filtered[summarized_filtered_BCR_processed_filtered$cluster=='non-Expanded' , ]
summarized_filtered_BCR_processed_filtered_mainClade= summarized_filtered_BCR_processed_filtered[summarized_filtered_BCR_processed_filtered$cluster=='Expanded' , ]




clone_list <- unique(df2$clone_name)


for (cl in clone_list) {
  
  
  data= summarized_filtered_BCR_processed_filtered %>% filter(clone_name == cl)  #summarized_filtered_BCR_processed_filtered_exclude_mainClade
  
  df =data.frame(cell_state= c("BcellHIGH", "Secretory" ,"Cycling","MetaboHIGH", "Inflammatory" ,"Neural"),
                 value = c(0,0,0,0,0,0),
                 color= color.tree)
  
  con= 0
  pie_colors_cell_state=vector()
  for (cs in c("BcellHIGH", "Secretory" ,"Cycling","MetaboHIGH", "Inflammatory" ,"Neural")){
    con = con +1
    ind1= which(data$new_cell_state == cs)
    
    df$value[con] = length(ind1)
    
  }
  
  ind_df = which(df$value != 0 )
  df = df[ind_df, ]
  df$color <- color.tree[df$cell_state]
  pie_colors_cell_state <- df$color
  
  pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/',sample,'_PIE_',DominantChain,'_cellState_',cl,'.pdf'))
  pie(df$value , df$cell_state ,col=pie_colors_cell_state)
  dev.off()
  
}


data= summarized_filtered_BCR_processed_filtered_mainClade
df =data.frame(cell_state= c("BcellHIGH", "Secretory" ,"Cycling","MetaboHIGH", "Inflammatory" ,"Neural"),
               value = c(0,0,0,0,0,0),
               color= color.tree)

con= 0
pie_colors_cell_state=vector()
for (cs in c("BcellHIGH", "Secretory" ,"Cycling","MetaboHIGH", "Inflammatory" ,"Neural")){
  con = con +1
  ind1= which(data$new_cell_state == cs)
  
  df$value[con] = length(ind1)
  
}

ind_df = which(df$value != 0 )
df = df[ind_df, ]
df$color <- color.tree[df$cell_state]
pie_colors_cell_state <- df$color

pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/',sample,'_PIE_',DominantChain,'_cellState_Expanded','.pdf'))
pie(df$value , df$cell_state ,col=pie_colors_cell_state)
dev.off()



node_depths <- node.depth.edgelength(tree)
# 
# # Get tip labels (these are usually your contig_ids)
tips <- tree$tip.label
# 
# # Extract timing for just the tips (contigs)
tip_depths <- node_depths[1:length(tips)]

timing_df <- data.frame(
  contig_id = tips,
  relative_timing_ape = tip_depths
)


summarized_filtered_BCR_processed_filtered = summarized_filtered_BCR_processed_filtered %>%
  left_join(timing_df, by = "contig_id")



library(tidyr)

summarized_filtered_BCR_processed_filtered <- 
  summarized_filtered_BCR_processed_filtered %>%
  group_by(clone_id) %>%
  fill(relative_timing_ape, .direction = "downup") %>%  # fill NA both down and up
  ungroup()

summarized_filtered_BCR_processed_filtered_exclude_mainClade= summarized_filtered_BCR_processed_filtered[summarized_filtered_BCR_processed_filtered$cluster=='non-Expanded' , ]


pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/',sample,'_NEW_colorClades_CELLSTATE_barplot_phylo_ExcludeLongerCluster_timing_fixed.pdf'))  
print(ggplot(summarized_filtered_BCR_processed_filtered_exclude_mainClade, aes(x=cluster, y=round(as.numeric(relative_timing_ape),4), fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree))

dev.off()
write.csv(summarized_filtered_BCR_processed_filtered,paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/','summarized_filtered_BCR_processed_igphyml_',sample,'.csv'))
summarized_filtered_BCR_processed_filtered=summarized_filtered_BCR_processed_filtered[order(summarized_filtered_BCR_processed_filtered$cluster),]

gc()


CDR3_data <- data.frame(
  contig_id = summarized_filtered_BCR_processed_filtered$contig_id,
  CDR3 = summarized_filtered_BCR_processed_filtered$CDR3,
  stringsAsFactors = FALSE
)
summarized_filtered_BCR_processed_filtered$cdr3_diversity = alignment_score(CDR3_data,paste0("./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/MAFFTaligned_StopCodonNNN_CDR3_sequences_",sample,"_new_corrected_validCodon_fixed.fasta"))$scores
gc()


summarized_filtered_BCR_processed_filtered_exclude_mainClade = summarized_filtered_BCR_processed_filtered[summarized_filtered_BCR_processed_filtered$cluster=='non-Expanded',]
pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/',sample,'_NEW_colorClades_CELLSTATE_barplot_phylo_ExcludeLongerCluster_CDR3diversity.pdf'))  
print(ggplot(summarized_filtered_BCR_processed_filtered_exclude_mainClade, aes(x=cluster, y=cdr3_diversity, fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree))

dev.off()

write.csv(summarized_filtered_BCR_processed_filtered,paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/','summarized_filtered_BCR_processed_igphyml_',sample,'.csv'))

### test expansion
node_depths <- node.depth.edgelength(expanded_tree_original)
# 
# # Get tip labels (these are usually your contig_ids)
tips <- expanded_tree_original$tip.label
# 
# # Extract timing for just the tips (contigs)
tip_depths <- node_depths[1:length(tips)]

timing_df <- data.frame(
  contig_id = tips,
  relative_timing_ape_original = tip_depths
)


summarized_filtered_BCR_processed_filtered = summarized_filtered_BCR_processed_filtered %>%
  left_join(timing_df, by = "contig_id")


summarized_filtered_BCR_processed_filtered <- 
  summarized_filtered_BCR_processed_filtered %>%
  group_by(clone_id) %>%
  fill(relative_timing_ape_original, .direction = "downup") %>%  # fill NA both down and up
  ungroup()

summarized_filtered_BCR_processed_filtered_exclude_mainClade= summarized_filtered_BCR_processed_filtered[summarized_filtered_BCR_processed_filtered$cluster=='non-Expanded' , ]


pdf(paste0('./igphyml_modelGY_NNNstopCodon_MAFFT_validCodon_semiRefined/singleCDR3_final/',sample,'_NEW_colorClades_CELLSTATE_barplot_phylo_ExcludeLongerCluster_timing_fixed_fullTreeTIming.pdf'))  
print(ggplot(summarized_filtered_BCR_processed_filtered_exclude_mainClade, aes(x=cluster, y=round(as.numeric(relative_timing_ape_original),4), fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree))

dev.off()

