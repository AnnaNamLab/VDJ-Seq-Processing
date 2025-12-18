
# Load required libraries
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

sample= "sample1"
input_directory = paste0(set_working_dir,'/',sample)
metadata_file = "metadata.cell.csv"
trust4_cluster_clone = Sys.glob('*clone*.tsv')[1] 
light_chain_cellStatus = Sys.glob(paste0('./doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv'))
heavy_chain_cellStatus = Sys.glob(paste0('./doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGH.csv'))

# -----------------------
# Set working directory
# -----------------------
dir.create(input_directory, showWarnings = FALSE)
setwd(input_directory)

### reading the cell types from GEX
new_clusters= read.csv(metadata_file)
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)


### extracting the Bcells nad HRS from GEX data per sample
sample_HRS_Bcells=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS", "Bcells"),]
sample_HRS=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS"),]

### after running TRUST4 clustering, we pull from HPC the outcome 
FILTERED_out_clone= read.table(trust4_cluster_clone,header= FALSE, sep='\t')
colnames(FILTERED_out_clone)=c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )


### from TRUST4 clustering we need to define the chain; V_gene is corresponded to V_gene
FILTERED_out_clone$chain=substr(FILTERED_out_clone$V_gene,1,3)
FILTERED_out_clone$cell_id1=substr(FILTERED_out_clone$contig_id,1,16)
FILTERED_out_clone= FILTERED_out_clone[FILTERED_out_clone$cell_id1 %in% sample_HRS_Bcells$cell_id1,]

### separating the locus per IGK, IGL, IGH

for (chain in c("IGL","IGK", "IGH")){
  sample_trust4_per_chain= FILTERED_out_clone[FILTERED_out_clone$chain==chain,]
    
  ### we an have VJ assignment and also getting the MainCelltype status from GEX data 
  if (chain=="IGL" | chain=="IGK"){
  
    pilot_cellStatus_file= read.csv(light_chain_cellStatus) 
    singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
    
    sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$V_gene," ",sample_trust4_per_chain$J_gene)
    ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
    sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind]
    
    sample_trust4_per_chain= sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% singlet,]
    
    sample_trust4_per_chain= sample_trust4_per_chain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
    if (chain=="IGL"){
      data_IGL=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS$cell_id1,]
      full_data_IGL= sample_trust4_per_chain
    }else if (chain=="IGK") {
      data_IGK=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS$cell_id1,]
      full_data_IGK= sample_trust4_per_chain
    }
    
  }else{
    
    pilot_cellStatus_file= read.csv(heavy_chain_cellStatus) 
    singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
    
    sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$V_gene," ",sample_trust4_per_chain$D_gene," ",sample_trust4_per_chain$J_gene)
    ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
    sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind] 
    sample_trust4_per_chain= sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% singlet,]
    sample_trust4_per_chain= sample_trust4_per_chain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
    data_IGH=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS$cell_id1,]
    
  }
  
  
  ### generating the disnatnce matrix
  data=sample_trust4_per_chain
  data=data.frame(data)
  data= data[order(data$VDJ, decreasing=TRUE),]
  dna_sample=char2dna(data$CDR3, simplify = FALSE) ## simplify= FALSE is the default
  ### to align the sequence
  class(dna_sample)
  X11= align(dna_sample)
  rownames(X11)<- data$contig_id
  
  ### we need to change the phylo class
  X111.phydat<-as.phyDat(X11)
  HLsample_X11= X11
  model_name="F81"
  
  ### genrating the distance matrix from aligned sequence
  dist.X111.phydat = dist.dna(X11, model = "F81", variance = FALSE,
                              gamma = FALSE, pairwise.deletion = TRUE,
                              base.freq = NULL, as.matrix = TRUE)# , indel=TRUE
  
  
  giveNAs = which(is.nan(as.matrix(dist.X111.phydat)),arr.ind=TRUE)
  
  if (nrow(giveNAs)>0){
    fix_dist.X111.phydat = dist.X111.phydat[-c(giveNAs[,1]),-c(giveNAs[,2])]
  }else if (nrow(giveNAs)==0){
    
    fix_dist.X111.phydat=dist.X111.phydat
    
  }
  
 
  ### calculate Levenshtein distance between two strings
  
    dist1= as.matrix(dist.X111.phydat)
    row_dist1= rownames(dist1)
    new_dat1_ind=match(row_dist1,data$contig_id)
    new_data= data[new_dat1_ind,]
    
    col_runif = colorRamp2(c(0,max(dist1,na.rm = TRUE)), c("gray","red"))
    
    column_ha = HeatmapAnnotation(cell_status = new_data$MainCelltype, V_gene= new_data$V_gene, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
    row_ha = rowAnnotation(cell_status= new_data$MainCelltype, j_gene= new_data$J_gene, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
    
 
    pdf(paste0(sample,'_heatmap_lv_distance_CDR3_',chain,'_beforeVJcorrection.pdf'),width=15,heigh=15)
    print(Heatmap(dist1,col= col_runif,top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,
                  show_row_names = FALSE,cluster_columns=TRUE,show_row_dend = TRUE, raster_resize_mat = TRUE))
    dev.off()
  
}

