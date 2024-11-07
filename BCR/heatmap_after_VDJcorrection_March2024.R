
library(dplyr)
library(ggplot2)
library(Biostrings)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(insect)
library(aphid)
library(ape)
library(phangorn)
### we first run trust4, then getting cdr3_raw.out file
### doublets and noise cells are removed. Only Bcells and HRS cells are included.
###setting the directory

sample="HL20"
DominantChain="IGL"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_August2024/',sample))


### reading the cell types from GEX
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
  
  singlet_low_confidential= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_low_confident_singlet_IGL_IGK_read_thre1.txt'))[1])
  singlet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_singlet_IGL_IGK_read_thre1.txt'))[2])
  doublet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_doublet_IGL_IGK_read_thre1.txt')))
  
  sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$corrected_Vgene," ",sample_trust4_per_chain$corrected_Jgene)
  ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
  sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind]
  sample_trust4_per_chain=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% union(singlet$V1, singlet_low_confidential$V1),]
  
}else{
  
  singlet_low_confidential= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_low_confident_singlet_IGH_read_thre1.txt'))[1])
  singlet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_singlet_IGH_read_thre1.txt'))[2])
  doublet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_doublet_IGH_read_thre1.txt')))
  
  sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$corrected_Vgene," ",sample_trust4_per_chain$corrected_Dgene," ",sample_trust4_per_chain$corrected_Jgene)
  ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
  sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind] 
  sample_trust4_per_chain=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% union(singlet$V1, singlet_low_confidential$V1),]
}

  
  data=sample_trust4_per_chain
  data=data.frame(data)
  data= data[order(data$VDJ, decreasing=TRUE),]
  
  write.csv(data,paste0('no_doublet_VDJfixed_clean_BCR_data_HRS_Bcells_',DominantChain,'.csv'))
  dna_sample=char2dna(data$CDR3, simplify = FALSE) ## simplify= FALSE is the default
  ##### to align the sequence
  class(dna_sample)
  X11= align(dna_sample)
  rownames(X11)<- data$contig_id
  ## we need to change the phylo class
  X111.phydat<-as.phyDat(X11)
  HLsample_X11= X11
  model_name="F81"
  #genrating the distance matrix from aligned sequence
  dist.X111.phydat = dist.dna(X11, model = model_name, variance = FALSE,
                              gamma = FALSE, pairwise.deletion = TRUE,
                              base.freq = NULL, as.matrix = TRUE)# , indel=TRUE
  
  
  giveNAs = which(is.nan(as.matrix(dist.X111.phydat)),arr.ind=TRUE)
  
  if (nrow(giveNAs)>0){
    fix_dist.X111.phydat = dist.X111.phydat[-c(giveNAs[,1]),-c(giveNAs[,2])]
  }else if (nrow(giveNAs)==0){
    
    fix_dist.X111.phydat=dist.X111.phydat
    
  }
  
  
  
  
  
  # V9 is the CDR3
  #calculate Levenshtein distance between two strings
  
  dist1= as.matrix(fix_dist.X111.phydat)
  row_dist1= rownames(dist1)
  new_dat1_ind=match(row_dist1,data$contig_id)
  new_data= data[new_dat1_ind,]
  
  col_runif = colorRamp2(c(0,max(dist1)), c("gray","red"))
  
  column_ha = HeatmapAnnotation(cell_status = new_data$MainCelltype, V_gene= new_data$V_gene)
  row_ha = rowAnnotation(cell_status= new_data$MainCelltype, j_gene= new_data$J_gene)
  

  pdf(paste0(sample,'BCR_CLONE_heatmap_lv_distance_CDR3_',DominantChain,'_afterVJcorrection_dend.pdf'), width=23, height=23)
  print(Heatmap(dist1,col= col_runif,top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,show_row_names = FALSE,cluster_columns=TRUE,show_row_dend = TRUE))
  dev.off() 
  