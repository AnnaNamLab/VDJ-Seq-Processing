
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

sample="HL8"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_August2024/',sample))
### reading the cell types from GEX


new_clusters= read.csv('/users/saramoein/downloads/2024-08-13_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
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

for (chain in c("IGL","IGK","IGH")){
  sample_trust4_per_chain= FILTERED_out_clone[FILTERED_out_clone$chain==chain,]
  

  #### we an have VJ assignment and also getting the MainCelltype status from GEX data 
  if (chain=="IGL" | chain=="IGK"){
    
    singlet_low_confidential= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_singlet_IGL_IGK_read_thre1.txt'))[1])
    singlet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_singlet_IGL_IGK_read_thre1.txt'))[2])
    doublet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_doublet_IGL_IGK_read_thre1.txt')))
    
     sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$V_gene," ",sample_trust4_per_chain$J_gene)
     ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
     sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind]
     sample_trust4_per_chain= sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% union(singlet_low_confidential$V1,singlet$V1),]
     sample_trust4_per_chain= sample_trust4_per_chain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
     if (chain=="IGL"){
       data_IGL=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS$cell_id1,]
       full_data_IGL= sample_trust4_per_chain
     }else if (chain=="IGK") {
       data_IGK=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS$cell_id1,]
       full_data_IGK= sample_trust4_per_chain
     }
     
    }else{
     
      singlet_low_confidential= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_singlet_IGH_read_thre1.txt'))[1])
      singlet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_singlet_IGH_read_thre1.txt'))[2])
      doublet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','*_doublet_IGH_read_thre1.txt')))
  
     sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$V_gene," ",sample_trust4_per_chain$D_gene," ",sample_trust4_per_chain$J_gene)
     ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
     sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind] 
     sample_trust4_per_chain= sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% union(singlet_low_confidential$V1,singlet$V1),]
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

# print(Heatmap(dist1, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
#       cluster_columns=TRUE,cluster_rows=TRUE, name = "mat", row_split = NULL,
#       show_row_dend = TRUE, 
pdf(paste0(sample,'_heatmap_lv_distance_CDR3_',chain,'_beforeVJcorrection_test2.pdf'), height=14, width=14)
print(Heatmap(dist1,col= col_runif,top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,show_row_names = FALSE,cluster_columns=TRUE,show_row_dend = TRUE))
dev.off()



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

