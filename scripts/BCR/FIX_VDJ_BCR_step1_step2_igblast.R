### this code has 2 main steps for fixing V/D/J assignment obtained from TRUST4. Step 1 is based on igblast, and step2 is based on fasta file to find any contig that has similar V/D/J as 
### dominant V/D/J (clonal V/D/J).

library(dplyr)
#########STEP1#######


### step1 receives all the HRS_contigs of a sample and generates the igblast VDJs from trust4 fasta file and then merges the generated igblast V/D/Js as new columns 
### to the trust4_clustered file. Then Any of the contigs that their igblast V/D/J is the same as dominant VDJ, then it is replaced. Else, the trust4 V/D/J is returened.

### input dominant V_D_J ; this is obtained from the heatmap in previous steps



sample="HL10"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample))

## reading the GEX cell type annotation 
cell_type_annotation=new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_Jan2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
cell_type_annotation$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", cell_type_annotation$Full.cell_id)
sample_cell_type_annotation=cell_type_annotation[cell_type_annotation$Patient== sample,]

## extracting the clustering results from TRUST4; we use this file for all analysis since it has all the complete CDR3s. This file contains TRUST4 V/D/J

FILTERED_out_clone=read.table('out_clone_clustered.tsv',sep='\t')
## extracting the dominant clone
View(FILTERED_out_clone) # 
dominant_V="IGKV2-28*01"
dominant_D="*"
dominant_J="IGKJ5*01"
chain=substr(dominant_V,1,3)
colnames(FILTERED_out_clone)=c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )

FILTERED_out_clone$chain=substr(FILTERED_out_clone$V_gene,1,3)
FILTERED_out_clone_dominantChain=FILTERED_out_clone[which(FILTERED_out_clone$chain==chain),]


### the cuurent VDJ gene columns are renamed to TRUST4_V/D/J_calles
colnames(FILTERED_out_clone_dominantChain)[3]<-"trust4_v_call"
colnames(FILTERED_out_clone_dominantChain)[4]<-"trust4_d_call"
colnames(FILTERED_out_clone_dominantChain)[5]<-"trust4_j_call"
colnames(FILTERED_out_clone_dominantChain)[14]<-"sequence_id"


### extracting the HRS contigs from data 
FILTERED_out_clone_dominantChain$cell_id1= substr(FILTERED_out_clone_dominantChain$sequence_id,1,16)
FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[FILTERED_out_clone_dominantChain$cell_id1 %in% sample_cell_type_annotation$cell_id1,]
FILTERED_out_clone_dominantChain_index= match(FILTERED_out_clone_dominantChain$cell_id1, sample_cell_type_annotation$cell_id1)
FILTERED_out_clone_dominantChain$malig_status=sample_cell_type_annotation$MainCelltype[FILTERED_out_clone_dominantChain_index]
FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[which(FILTERED_out_clone_dominantChain$malig_status %in% c("HRS","Bcells")),]
HRS_FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[which(FILTERED_out_clone_dominantChain$malig_status=="HRS"),]
HRS_FILTERED_out_clone_dominantChain= HRS_FILTERED_out_clone_dominantChain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 


##### keeping only singlet and low_confident_singlet
# singlet_low_confidential= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_low_confident_singlet_IGL_IGK_read_thre1.txt'))
# singlet= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_singlet_IGL_IGK_read_thre1.txt'))
# doublet= read.table(paste0('./new_new_doublets_fixed/combined_light_chains_automated/','2024-09-23','_doublet_IGL_IGK_read_thre1.txt'))

pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/REPEAT_other_doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))

singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])

HRS_FILTERED_out_clone_dominantChain= HRS_FILTERED_out_clone_dominantChain[HRS_FILTERED_out_clone_dominantChain$cell_id1 %in% (singlet),]

### HRS contigs are used for the next step for igblasting
write.table(HRS_FILTERED_out_clone_dominantChain$sequence_id,paste0('contigs_',sample,'_',chain,'.txt'),row.names= FALSE, col.names= FALSE, quote= FALSE)


### Bash script running on HPC
### filtering the airr.tsv file (from TRUST4) based on cells_HRS. That gives the contigs of HRS cells.
### Then the contigs are extracted in the text file. For example, for HL10, we will have contigs_HL10.txt


#To fix the VDJs first we add a new column to data based on igblast tool. To find the clone VDJ, for those non-clonal HRS cell's contigs, we try to ge
# #!/bin/bash
# 
# #SBATCH --partition=scu-cpu   # cluster-specific
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --job-name=igblast
# #SBATCH --time=48:00:00   # HH/MM/SS
# #SBATCH --mem=512G   # memory requested, units available: K,M,G,T
# 
# fastafile=./test/out_FR2_annot.fa
## we also convertthe fasta file to single line fatsa for using in the second step
# awk '/^>/&&NR>1{print "";}{printf "%s",/^>/ ? $0" " : $0}' ./test/out_FR2_annot.fa > ./test/singleline_HL10_s1s2_annot.fa
# 
# cd /athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin
# 
# sed '/^>/ s/ .*//' ./test/out_FR2_annot.fa > ./test/fix_my_sample.fa ### this line fixes the fasta file headers and only keeps the contig id
# grep -w -A 1 -f  ./test/contigs_HL10.txt ./test/fix_my_sample.fa --no-group-separator > ./test/sub_fix_my_sample.fa.  ### this line subset the fasta file to keep only HRS_contigs
# 
# igblastn -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGHV_clean.fasta -num_alignments_V 3 -germline_db_D /home/sam4032/share/igblast/database/imgt_human_IGHD_clean.fasta -num_alignments_D 3 -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGHJ_clean.fasta -num_alignments_J 3 -query /athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin/test/sub_fix_my_sample.fa -organism human -outfmt 19 -out /athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin/test/igblast_HL10_output2.txt
# igblastn -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGHV_clean.fasta -num_alignments_V 3 -germline_db_D /home/sam4032/share/igblast/database/imgt_human_IGHD_clean.fasta -num_alignments_D 3 -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGHJ_clean.fasta -num_alignments_J 3 -query /athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin/test/sub_fix_my_sample.fa -organism human -outfmt 7 -out /athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin/test/igblast_HL10_output2_7.txt
# 



# reading the igbalst table; the output of igblast is airr
igblast_output=read.table(paste0(sample,'_igblast_output2_',chain,'.txt'),header=T,sep='\t')

length(HRS_FILTERED_out_clone_dominantChain$cell_id1)
common=intersect(substr(igblast_output$sequence_id ,1,16) , HRS_FILTERED_out_clone_dominantChain$cell_id1)
length(common)
## renames in VDJ gene columns to igblast_V/D/J_calls
colnames(igblast_output)[10]="igblast_v_call"
colnames(igblast_output)[11]="igblast_d_call"
colnames(igblast_output)[12]="igblast_j_call"

## merging TRUST4_cluster file with igblast table

merge_igblast_trust4_all = merge(HRS_FILTERED_out_clone_dominantChain,igblast_output,by="sequence_id")
merge_igblast_trust4_all$igblast_d_call[is.na(merge_igblast_trust4_all$igblast_d_call) == TRUE] <- "*"
merge_igblast_trust4=merge_igblast_trust4_all#cbind(merge_igblast_trust4_all[,c(1:30)])



### the idea is to generate a new column that select the between igblast ot trust4 
merge_igblast_trust4$trust4_igblast_v_call=NA
merge_igblast_trust4$trust4_igblast_d_call=NA
merge_igblast_trust4$trust4_igblast_j_call=NA

for (i in 1:nrow(merge_igblast_trust4)){
  
  if (grepl(dominant_V,merge_igblast_trust4$igblast_v_call[i],fixed =TRUE) & grepl(dominant_D , merge_igblast_trust4$igblast_d_call[i],fixed =TRUE) & grepl(dominant_J,merge_igblast_trust4$igblast_j_call[i] ,fixed =TRUE)){
    merge_igblast_trust4$trust4_igblast_v_call[i]= dominant_V
    merge_igblast_trust4$trust4_igblast_j_call[i]= dominant_J 
    merge_igblast_trust4$trust4_igblast_d_call[i]= dominant_D  
  }else{
    merge_igblast_trust4$trust4_igblast_v_call[i]= merge_igblast_trust4$trust4_v_call[i] 
    merge_igblast_trust4$trust4_igblast_j_call[i]= merge_igblast_trust4$trust4_j_call[i] 
    merge_igblast_trust4$trust4_igblast_d_call[i]= merge_igblast_trust4$trust4_d_call[i]
  }
  
  
}

#########STEP2#######

# Step2 aims to select any rows in contig in TRUST4 fasta file that has common V/D/J with dominant V/D/J. We should extract top 3 V/D/J gebes from igbalst.
## any rows in igblast results that contains the dominant V/D/J  can give us top 3 V/D/J genes
## all those contigs are extracted from fasta file


fasta_subset=read.table('singleline_TRUST_all_BCR_R2_annot.fa')

### extract one full sequence from the dominant VDJ and extract from igblast v1,v2,v3,j1,j2,j3,d1,d2,d3

v1="IGKV2-28*01"
v2="IGKV2D-28*01"
v3="IGKV2-40*01"
d1=NA
d2=NA
d3=NA
j1="IGKJ5*01"
j2="IGKJ1*01"
j3="IGKJ3*01"



merge_igblast_trust4$fix_trust4_fasta_v_call="NA"
merge_igblast_trust4$fix_trust4_fasta_j_call="NA"
merge_igblast_trust4$fix_trust4_fasta_d_call="NA"
sel_dat =data.frame()

if ((!is.na(v1) & !is.na(j1) & chain!="IGH") | (!is.na(v1) & !is.na(j1)  & !is.na(d1) & chain=="IGH")){
      
      
  
      fasta_subset1= fasta_subset[grepl(v1, as.character(fasta_subset$V4),fix=TRUE),]
      fasta_subset2= fasta_subset[grepl(v2, as.character(fasta_subset$V4),fix=TRUE),]
      fasta_subset3= fasta_subset[grepl(v3, as.character(fasta_subset$V4),fix=TRUE),]
      data= rbind(fasta_subset1,fasta_subset2,fasta_subset3)
      
 
      
      fasta_subset4= data[grepl(j1, as.character(data$V6),fix=TRUE),]
      fasta_subset5= data[grepl(j2, as.character(data$V6),fix=TRUE),]
      fasta_subset6= data[grepl(j3, as.character(data$V6),fix=TRUE),]
      
      data= rbind(fasta_subset4,fasta_subset5,fasta_subset6)
      
      if (chain=="IGH"){
       fasta_subset7= data[grepl(d1, as.character(data$V5),fix=TRUE),]
       fasta_subset8= data[grepl(d2, as.character(data$V5),fix=TRUE),]
       fasta_subset9= data[grepl(d3, as.character(data$V5),fix=TRUE),]
       data= rbind(fasta_subset7,fasta_subset8,fasta_subset9)
      }
      sel_data= data
      sel_data$cell_id1= substr(sel_data$V1,2,17)
      sel_data$sequence_id= substr(sel_data$V1,2,50)
      
      sel_dat=distinct(sel_data)
      
      # Here we replace dominant V/D/J for any of the contigs that are available in the fasta_subset with equal V/D/J as in dominant V/D/J
      ind=match(merge_igblast_trust4$sequence_id , sel_data$sequence_id)
      
      merge_igblast_trust4$fix_trust4_fasta_v_call[which(ind>0)]= dominant_V
      merge_igblast_trust4$fix_trust4_fasta_d_call[which(ind>0)]= dominant_D
      merge_igblast_trust4$fix_trust4_fasta_j_call[which(ind>0)]= dominant_J
      
}

#column "final_fixed_v_call" is the final V/D/J calls. We initiate with the final V/D/J calls from STEP1 and then continue final assignment
## then selecting between fixed_trust4 and igblast_trust4 columns from step1

merge_igblast_trust4$final_fixed_v_call=merge_igblast_trust4$trust4_igblast_v_call
merge_igblast_trust4$final_fixed_j_call=merge_igblast_trust4$trust4_igblast_j_call
merge_igblast_trust4$final_fixed_d_call=merge_igblast_trust4$trust4_igblast_d_call

for (i in 1:nrow(merge_igblast_trust4)){ #  & (merge_igblast_trust4$fix_trust4_fasta_v_call[i]!=NA)
  if ((merge_igblast_trust4$trust4_igblast_v_call[i]==dominant_V) | (merge_igblast_trust4$fix_trust4_fasta_v_call[i]==dominant_V) & (is.na(merge_igblast_trust4$fix_trust4_fasta_v_call[i])==FALSE)){
    merge_igblast_trust4$final_fixed_v_call[i]=dominant_V
  }
  if ((merge_igblast_trust4$trust4_igblast_d_call[i]==dominant_D) | (merge_igblast_trust4$fix_trust4_fasta_d_call[i]==dominant_D) & (is.na(merge_igblast_trust4$fix_trust4_fasta_d_call[i])==FALSE)){
    merge_igblast_trust4$final_fixed_d_call[i]=dominant_D
    
  }
  if ((merge_igblast_trust4$trust4_igblast_j_call[i]==dominant_J) | (merge_igblast_trust4$fix_trust4_fasta_j_call[i]==dominant_J) & (is.na(merge_igblast_trust4$fix_trust4_fasta_j_call[i])==FALSE)){
    merge_igblast_trust4$final_fixed_j_call[i]=dominant_J
    
  }  
  
  
}


## final file with fixed V/D/J in column "final_fixed_v/d/j_call"


write.csv(merge_igblast_trust4,paste0(sample,'_merge_igblast_trust4',chain,'.csv'))


FILTERED_out_clone_dominantChain$corrected_Vgene=FILTERED_out_clone_dominantChain$trust4_v_call
FILTERED_out_clone_dominantChain$corrected_Dgene=FILTERED_out_clone_dominantChain$trust4_d_call
FILTERED_out_clone_dominantChain$corrected_Jgene=FILTERED_out_clone_dominantChain$trust4_j_call

for (i in 1:nrow(FILTERED_out_clone_dominantChain)){
  ind= which(FILTERED_out_clone_dominantChain$sequence_id[i] == merge_igblast_trust4$sequence_id)
  
  if (length(ind)>0) {
    FILTERED_out_clone_dominantChain$corrected_Vgene[i]= merge_igblast_trust4$final_fixed_v_call[ind]
    FILTERED_out_clone_dominantChain$corrected_Dgene[i]= merge_igblast_trust4$final_fixed_d_call[ind]
    FILTERED_out_clone_dominantChain$corrected_Jgene[i]= merge_igblast_trust4$final_fixed_j_call[ind]
    
  }else {
    FILTERED_out_clone_dominantChain$corrected_Vgene[i]= FILTERED_out_clone_dominantChain$trust4_v_call[i]
    FILTERED_out_clone_dominantChain$corrected_Dgene[i]= FILTERED_out_clone_dominantChain$trust4_d_call[i]
    FILTERED_out_clone_dominantChain$corrected_Jgene[i]= FILTERED_out_clone_dominantChain$trust4_j_call[i]
    
    
  } 
  
  if (FILTERED_out_clone_dominantChain$trust4_v_call[i] %in% c(v1,v2,v3) & FILTERED_out_clone_dominantChain$malig_status[i]=="HRS"){
    FILTERED_out_clone_dominantChain$corrected_Vgene[i]=dominant_V
    
  }
  if (FILTERED_out_clone_dominantChain$trust4_d_call[i] %in% c(d1,d2,d3) & FILTERED_out_clone_dominantChain$malig_status[i]=="HRS"){
    FILTERED_out_clone_dominantChain$corrected_Dgene[i]=dominant_D
    
  }
  
  if (FILTERED_out_clone_dominantChain$trust4_j_call[i] %in% c(j1,j2,j3) & FILTERED_out_clone_dominantChain$malig_status[i]=="HRS"){
    FILTERED_out_clone_dominantChain$corrected_Jgene[i]=dominant_J
    
  }
  
  
  
}
    
write.csv(FILTERED_out_clone_dominantChain,paste0(sample,'_FILTERED_out_clone_dominantChain_',chain,'.csv'))

