library(dplyr)
library(ggplot2)
library(Biostrings)
library(circlize)
library(ComplexHeatmap)

### we first run trust4, then getting cdr3_raw.out file
### doublets and noise cells are removed. Only Bcells and HRS cells are included.
###setting the directory
sample="HL20"
DominantChain="IGK"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_August2024/',sample))


### reading the cell types from GEX
# new_clusters=read.csv('/users/saramoein/downloads/2023-09-14_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
# new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)
new_clusters= read.csv('/users/saramoein/downloads/2024-08-13_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)



#### extracting the Bcells nad HRS from GEX data per sample
sample_HRS_Bcells=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS", "Bcells"),]

### after running TRUST4 clustering, we pull from HPC the outcome; this file contains the corrected V/D/J genes onlu for HRS cells. So we need to merge it with the full file output from trust4 clustering
sample_trust4_per_chain= read.csv(paste0(sample,'_FILTERED_out_clone_dominantChain_',DominantChain,'.csv'))
sample_trust4_per_chain= sample_trust4_per_chain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
colnames(sample_trust4_per_chain)[1:16]=c("X","consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )


#### we have VJ assignment and also getting the MainCelltype status from GEX data 
if (DominantChain=="IGL" | DominantChain=="IGK"){
  
  singlet_low_confidential= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_low_confident_singlet_IGL_IGK_read_thre1.txt'))
  singlet= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_singlet_IGL_IGK_read_thre1.txt'))
  doublet= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_doublet_IGL_IGK_read_thre1.txt'))
  
  sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$corrected_Vgene," ",sample_trust4_per_chain$corrected_Jgene)
  ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
  sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind]
  sample_trust4_per_chain=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% union(singlet$V1, singlet_low_confidential$V1),]
  
}else{
  
  singlet_low_confidential= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_low_confident_singlet_IGH_read_thre1.txt'))
  singlet= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_singlet_IGH_read_thre1.txt'))
  doublet= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_doublet_IGH_read_thre1.txt'))
  
  sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$corrected_Vgene," ",sample_trust4_per_chain$corrected_Dgene," ",sample_trust4_per_chain$corrected_Jgene)
  ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
  sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind] 
  sample_trust4_per_chain=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% union(singlet$V1, singlet_low_confidential$V1),]
}


merge_data= sample_trust4_per_chain

merge_data=data.frame(merge_data)
merge_data= merge_data[order(merge_data$VDJ, decreasing=TRUE),]

write.csv(merge_data,paste0('no_doublet_VDJfixed_clean_BCR_data_HRS_Bcells_',DominantChain,'.csv'))

  #### V9 is the CDR3 
  #calculate Levenshtein distance between two strings
  dist1=as.matrix(stringDist(merge_data$CDR3, method = "levenshtein"))
  
  col_runif = colorRamp2(c(0, max(dist1)), c("gray", "red"))
  pdf(paste0(sample,'BCR_CLONE_heatmap_lv_distance_CDR3_',DominantChain,'_afterVJcorrection_dend.pdf'), width=23, height=23)
  column_ha = HeatmapAnnotation(cell_status = merge_data$malig_status, corrected_Vgene= merge_data$corrected_Vgene)
  row_ha = rowAnnotation(cell_status= merge_data$malig_status, corrected_Jgene= merge_data$corrected_Jgene)
  print(Heatmap(dist1,col= col_runif,top_annotation = column_ha, right_annotation = row_ha))
  dev.off()
  graphics.off()
  
  
  
  