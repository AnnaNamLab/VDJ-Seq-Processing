################################################################################
# Purpose:
# This script identifies and corrects BCR V/D/J gene calls for contigs belonging
# to a dominant (clonal) V/D/J using multiple evidence sources:
# 1) TRUST4 clustering output
# 2) IgBlast AIRR output
# 3) FASTA-level V/D/J annotation
# 4) Cell-type annotation (HRS vs B cells)
# 5) Singlet/doublet filtering
#
# The final goal is to enforce consistent dominant V/D/J calls in malignant HRS
################################################################################

library(dplyr)

################################################################################
# STEP 0: Sample setup and dominant clone definition
################################################################################

sample="sample1"
setwd(paste0('set_working_directory/',sample))

# Read TRUST4 clustering output
FILTERED_out_clone=read.table('cluster_clone.tsv',sep='\t')

# Define dominant (clonal) V/D/J genes
# D is set to "*" for light chain

dominant_V="IGLV4-69*01"
dominant_D="*"
dominant_J="IGLJ1*01"

# Infer chain type (IGH / IGK / IGL) from V gene prefix
chain=substr(dominant_V,1,3)

################################################################################
# STEP 1: Cell-type annotation and TRUST4 contig preprocessing
################################################################################

# Read cell-type annotation (GEX metadata)
cell_type_annotation=read.csv('cell_metadata.csv')
# Extract short cell barcode (used for matching TRUST4 contigs)
cell_type_annotation$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", cell_type_annotation$Full.cell_id)

# Keep only annotations for the current sample
sample_cell_type_annotation=cell_type_annotation[cell_type_annotation$Patient== sample,]

# Rename TRUST4 columns to meaningful names
colnames(FILTERED_out_clone)=c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )

# Assign chain based on V gene
FILTERED_out_clone$chain=substr(FILTERED_out_clone$V_gene,1,3)

# Keep only contigs from the dominant chain
FILTERED_out_clone_dominantChain=FILTERED_out_clone[which(FILTERED_out_clone$chain==chain),]


# Rename TRUST4 V/D/J columns
colnames(FILTERED_out_clone_dominantChain)[3]<-"trust4_v_call"
colnames(FILTERED_out_clone_dominantChain)[4]<-"trust4_d_call"
colnames(FILTERED_out_clone_dominantChain)[5]<-"trust4_j_call"

# Rename contig identifier
colnames(FILTERED_out_clone_dominantChain)[14]<-"sequence_id"


################################################################################
# STEP 1.1: Match contigs to cell types and keep HRS/B cells only
################################################################################

# Extract cell barcode from contig ID
FILTERED_out_clone_dominantChain$cell_id1= substr(FILTERED_out_clone_dominantChain$sequence_id,1,16)

# Keep contigs with annotated cell types
FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[FILTERED_out_clone_dominantChain$cell_id1 %in% sample_cell_type_annotation$cell_id1,]

# Map malignancy status to contigs
FILTERED_out_clone_dominantChain_index= match(FILTERED_out_clone_dominantChain$cell_id1, sample_cell_type_annotation$cell_id1)
FILTERED_out_clone_dominantChain$malig_status=sample_cell_type_annotation$MainCelltype[FILTERED_out_clone_dominantChain_index]

# Keep only HRS and B cells
FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[which(FILTERED_out_clone_dominantChain$malig_status %in% c("HRS","Bcells")),]

# For HRS cells: keep the contig with the highest read count per cell
HRS_FILTERED_out_clone_dominantChain=FILTERED_out_clone_dominantChain[which(FILTERED_out_clone_dominantChain$malig_status=="HRS"),]
HRS_FILTERED_out_clone_dominantChain= HRS_FILTERED_out_clone_dominantChain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 


################################################################################
# STEP 1.2: Remove doublets and keep singlets only
################################################################################

pilot_cellStatus_file= read.csv(Sys.glob(paste0('./doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))

# Extract singlet cell barcodes
singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])

# Keep only singlet HRS contigs
HRS_FILTERED_out_clone_dominantChain= HRS_FILTERED_out_clone_dominantChain[HRS_FILTERED_out_clone_dominantChain$cell_id1 %in% (singlet),]

# Export contig IDs for IgBlast input
write.table(HRS_FILTERED_out_clone_dominantChain$sequence_id,paste0('contigs_',sample,'_',chain,'.txt'),row.names= FALSE, col.names= FALSE, quote= FALSE)



################################################################################
# STEP 1.3: Read and merge IgBlast AIRR output
################################################################################

igblast_output= read.csv(paste0(sample,'_igblast_output2_IGL_correrct.csv'), 
                         header = TRUE,     # or FALSE if no header
                         fill = TRUE,       # fills missing columns with NA
                         blank.lines.skip = FALSE,  # do NOT skip empty lines
                         stringsAsFactors = FALSE,  # keep text as strings
                         na.strings = c("", "NA"))   # treat empty cells as NA

# Remove empty sequences
igblast_output= igblast_output[is.na(igblast_output$sequence)==FALSE, ]
length(HRS_FILTERED_out_clone_dominantChain$cell_id1)
common=intersect(substr(igblast_output$sequence_id ,1,16) , HRS_FILTERED_out_clone_dominantChain$cell_id1)
length(common)

# Rename IgBlast V/D/J columns
colnames(igblast_output)[10]="igblast_v_call"
colnames(igblast_output)[11]="igblast_d_call"
colnames(igblast_output)[12]="igblast_j_call"

# Merge TRUST4 and IgBlast by contig ID
merge_igblast_trust4_all = merge(HRS_FILTERED_out_clone_dominantChain,igblast_output,by="sequence_id")

# Replace missing D calls with "*"
merge_igblast_trust4_all$igblast_d_call[is.na(merge_igblast_trust4_all$igblast_d_call) == TRUE] <- "*"
merge_igblast_trust4=merge_igblast_trust4_all#cbind(merge_igblast_trust4_all[,c(1:30)])

################################################################################
# STEP 1.4: TRUST4 vs IgBlast reconciliation
################################################################################

# Initialize reconciled V/D/J columns
merge_igblast_trust4$trust4_igblast_v_call=NA
merge_igblast_trust4$trust4_igblast_d_call=NA
merge_igblast_trust4$trust4_igblast_j_call=NA

# Prefer IgBlast only if it matches the dominant clone
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

################################################################################
# STEP 2: FASTA-level rescue of dominant V/D/J
################################################################################

fasta_subset=read.table('singleline_TRUST_all_BCR_R2_annot.fa')

## extract one full sequence from the dominant VDJ and extract from igblast v1,v2,v3,j1,j2,j3,d1,d2,d3

v1="IGLV4-69*01"
v2="IGLV4-69*02"
v3="IGLV4-60*033"
d1=NA
d2=NA
d3=NA
j1="IGLJ1*01"
j2="IGLJ6*01"
j3="IGLJ2*01"

merge_igblast_trust4$fix_trust4_fasta_v_call="NA"
merge_igblast_trust4$fix_trust4_fasta_j_call="NA"
merge_igblast_trust4$fix_trust4_fasta_d_call="NA"
sel_dat =data.frame()

if ((!is.na(v1) & !is.na(j1) & chain!="IGL") | (!is.na(v1) & !is.na(j1)  & !is.na(d1) & chain=="IGL")){ 
  fasta_subset1= fasta_subset[grepl(v1, as.character(fasta_subset$V4),fix=TRUE),]
  fasta_subset2= fasta_subset[grepl(v2, as.character(fasta_subset$V4),fix=TRUE),]
  fasta_subset3= fasta_subset[grepl(v3, as.character(fasta_subset$V4),fix=TRUE),]
  data= rbind(fasta_subset1,fasta_subset2,fasta_subset3)
  fasta_subset4= data[grepl(j1, as.character(data$V6),fix=TRUE),]
  fasta_subset5= data[grepl(j2, as.character(data$V6),fix=TRUE),]
  fasta_subset6= data[grepl(j3, as.character(data$V6),fix=TRUE),]
  
  data= rbind(fasta_subset4,fasta_subset5,fasta_subset6)
  
  if (chain=="IGL"){
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


write.csv(merge_igblast_trust4,paste0(sample,'_merge_igblast_trust4_corrected',chain,'.csv'))


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

write.csv(FILTERED_out_clone_dominantChain,paste0(sample,'_FILTERED_out_clone_dominantChain_',chain,'_corrected.csv'))

