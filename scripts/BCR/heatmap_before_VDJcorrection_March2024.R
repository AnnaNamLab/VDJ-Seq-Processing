
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

### generating the entropy of the sequences. IF the entropy is high, then the biology model for very different number gives NaN or zeros.
calculate_entropy <- function(alignment_matrix) {
  alignment_length <- ncol(alignment_matrix)  # Number of positions
  entropies <- numeric(alignment_length)  # Initialize entropy vector
  
  for (pos in 1:alignment_length) {
    column <- alignment_matrix[, pos]  # Extract column
    counts <- table(column)  # Count occurrences of each base or gap
    freqs <- counts / sum(counts)  # Calculate frequencies
    entropies[pos] <- -sum(freqs * log2(freqs), na.rm = TRUE)  # Compute entropy
  }
  
  return(entropies)
}

sample="HL10"
dir.create(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample))
setwd(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample))
### reading the cell types from GEX


new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_Jan2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)


#### extracting the Bcells nad HRS from GEX data per sample
sample_HRS_Bcells=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS", "Bcells"),]
sample_HRS=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS"),]
### after running TRUST4 clustering, we pull from HPC the outcome 
FILTERED_out_clone= read.table(Sys.glob('*clone*.tsv')[1],header= FALSE, sep='\t')
colnames(FILTERED_out_clone)=c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )



### from TRUST4 clustering we need to define the chain; V_gene is corresponded to V_gene
FILTERED_out_clone$chain=substr(FILTERED_out_clone$V_gene,1,3)
FILTERED_out_clone$cell_id1=substr(FILTERED_out_clone$contig_id,1,16)
FILTERED_out_clone= FILTERED_out_clone[FILTERED_out_clone$cell_id1 %in% sample_HRS_Bcells$cell_id1,]



### separating the locus per IGK, IGL, IGH
#IGK
for (chain in c("IGL","IGK", "IGH")){
  sample_trust4_per_chain= FILTERED_out_clone[FILTERED_out_clone$chain==chain,]
  

  #### we an have VJ assignment and also getting the MainCelltype status from GEX data 
  if (chain=="IGL" | chain=="IGK"){
    
     # # 
     pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/REPEAT_other_doublets_BCR_thre07_ent08//',sample,'_Rawdata_IGK_IGL.csv')))
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
      pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/REPEAT_other_doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGH.csv')))
      
      singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
   
     sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$V_gene," ",sample_trust4_per_chain$D_gene," ",sample_trust4_per_chain$J_gene)
     ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
     sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind] 
     sample_trust4_per_chain= sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% singlet,]
     sample_trust4_per_chain= sample_trust4_per_chain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
     data_IGH=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS$cell_id1,]
     
   }

  
   
  
## generating the disnatnce matrix
data=sample_trust4_per_chain
data=data.frame(data)
data= data[order(data$VDJ, decreasing=TRUE),]



dna_sample=char2dna(data$CDR3, simplify = FALSE) ## simplify= FALSE is the default
##### to align the sequence
class(dna_sample)
X11= align(dna_sample)
rownames(X11)<- data$contig_id

#entropy_seq = mean(calculate_entropy(X11), na.rm=TRUE)
## we need to change the phylo class
X111.phydat<-as.phyDat(X11)
HLsample_X11= X11
model_name="F81"

#genrating the distance matrix from aligned sequence
dist.X111.phydat = dist.dna(X11, model = "F81", variance = FALSE,
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

if(sd(dist.X111.phydat, na.rm=TRUE)){
        dist1= as.matrix(dist.X111.phydat)
        row_dist1= rownames(dist1)
        new_dat1_ind=match(row_dist1,data$contig_id)
        new_data= data[new_dat1_ind,]
        
        col_runif = colorRamp2(c(0,max(dist1,na.rm = TRUE)), c("gray","red"))
        
        column_ha = HeatmapAnnotation(cell_status = new_data$MainCelltype, V_gene= new_data$V_gene, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
        row_ha = rowAnnotation(cell_status= new_data$MainCelltype, j_gene= new_data$J_gene, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
        
        # print(Heatmap(dist1, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        #       cluster_columns=TRUE,cluster_rows=TRUE, name = "mat", row_split = NULL,
        #       show_row_dend = TRUE, 
        pdf(paste0(sample,'_heatmap_lv_distance_CDR3_',chain,'_beforeVJcorrection.pdf'),width=15,heigh=15)
          print(Heatmap(dist1,col= col_runif,top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,
                      show_row_names = FALSE,cluster_columns=TRUE,show_row_dend = TRUE, raster_resize_mat = TRUE))
        dev.off()
} else {
  pdf(paste0(sample,'_heatmap_lv_distance_CDR3_',chain,'_beforeVJcorrection.pdf'),width=15,heigh=15)
    heatmap(dist.X111.phydat, main = paste0("entropy= ", entropy_seq, " and std distance matrix = ", sd(dist.X111.phydat, na.rm=TRUE)))
  dev.off()
}

}

### WE HAVE A VJ CORRECTION STEP....
### THIS STEP ADDS A COLUMN NEW_VJ ASSIGNMENT


join_IGL_IGK= full_join(data_IGL, data_IGK, by="cell_id1")

df1 <- join_IGL_IGK %>% mutate(IGK_read = if_else(is.na(read_fragment_count.y), 0,read_fragment_count.y ))

df1 <- df1 %>% mutate(IGL_read = if_else(is.na(read_fragment_count.x), 0,read_fragment_count.x ))
pdf('join_IGL_IGK_read_HRScells_zoomed.pdf')
lim= max(IQR(df1$IGL_read),IQR(df1$IGK_read))
plot(df1$IGL_read,df1$IGK_read, xlim=c(0,500) , ylim=c(0,500))
dev.off()
pdf('join_IGL_IGK_read_HRScells.pdf')

plot(df1$IGL_read,df1$IGK_read)
dev.off()

write.csv(table(data_IGK$MainCelltype),'data_IGK_HRS.csv')

write.csv(table(data_IGL$MainCelltype),'data_IGL_HRS.csv')

write.csv(table(data_IGH$MainCelltype),'data_IGH_HRS.csv')

table(data_IGK$MainCelltype)
table(data_IGL$MainCelltype)
table(data_IGH$MainCelltype)

