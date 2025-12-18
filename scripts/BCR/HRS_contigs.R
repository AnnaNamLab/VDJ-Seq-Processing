

library(dplyr)


# this step is for generating the HRS_contigs required for IgBlast run

sample="sample1"
setwd(paste0('set_working_directory/',sample))
FILTERED_out_clone=read.table('cluster_clone.tsv',sep='\t') # output from trust-cluster.py

## extracting the dominant clone
View(FILTERED_out_clone) # 
dominant_V="IGLV4-69*01" 
dominant_D="*"
dominant_J="IGLJ1*01"
chain=substr(dominant_V,1,3)

## reading the GEX cell type annotation 
cell_type_annotation=read.csv('cell_metadata.csv')
cell_type_annotation$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", cell_type_annotation$Full.cell_id)
sample_cell_type_annotation=cell_type_annotation[cell_type_annotation$Patient== sample,]

## extracting the clustering results from TRUST4; we use this file for all analysis since it has all the complete CDR3s. This file contains TRUST4 V/D/J

colnames(FILTERED_out_clone)=c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )

FILTERED_out_clone$chain=substr(FILTERED_out_clone$V_gene,1,3)
FILTERED_out_clone_dominantChain=FILTERED_out_clone[which(FILTERED_out_clone$chain==chain),]


## the cuurent VDJ gene columns are renamed to TRUST4_V/D/J_calles
colnames(FILTERED_out_clone_dominantChain)[3]<-"trust4_v_call"
colnames(FILTERED_out_clone_dominantChain)[4]<-"trust4_d_call"
colnames(FILTERED_out_clone_dominantChain)[5]<-"trust4_j_call"
colnames(FILTERED_out_clone_dominantChain)[14]<-"sequence_id"


## extracting the HRS contigs from data 
FILTERED_out_clone_dominantChain$cell_id1= substr(FILTERED_out_clone_dominantChain$sequence_id,1,16)
FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[FILTERED_out_clone_dominantChain$cell_id1 %in% sample_cell_type_annotation$cell_id1,]
FILTERED_out_clone_dominantChain_index= match(FILTERED_out_clone_dominantChain$cell_id1, sample_cell_type_annotation$cell_id1)
FILTERED_out_clone_dominantChain$malig_status=sample_cell_type_annotation$MainCelltype[FILTERED_out_clone_dominantChain_index]
FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[which(FILTERED_out_clone_dominantChain$malig_status %in% c("HRS","Bcells")),]
HRS_FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[which(FILTERED_out_clone_dominantChain$malig_status=="HRS"),]
HRS_FILTERED_out_clone_dominantChain= HRS_FILTERED_out_clone_dominantChain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 


pilot_cellStatus_file= read.csv(Sys.glob(paste0('./doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))

singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])

HRS_FILTERED_out_clone_dominantChain= HRS_FILTERED_out_clone_dominantChain[HRS_FILTERED_out_clone_dominantChain$cell_id1 %in% (singlet),]

## HRS contigs are used for the next step for igblasting
write.table(HRS_FILTERED_out_clone_dominantChain$sequence_id,paste0('contigs_',sample,'_',chain,'.txt'),row.names= FALSE, col.names= FALSE, quote= FALSE)

