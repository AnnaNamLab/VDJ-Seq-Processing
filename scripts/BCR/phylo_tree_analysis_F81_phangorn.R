#### THIS CODE IS TO GENERATE THE PHYLOGENETIC TREE, AFTER FIXING THE VDJS. THE TREE IS BASED ON ALL HRS (only HRS). TREES ARE BASED ON NO GERMLINE OR REFERENCE. 
rm(list=ls())
library(dplyr)
library(phangorn)
library(insect)
library(aphid)
library(ggtree)
library(RColorBrewer)
library(tidytree)
library(stringdist)
library(colorRamp2)
library(ComplexHeatmap)
library(seqinr)
library(TreeTools)
library(Biostrings)
library(stringr)
library(tidyverse)
library(ape)
library(DECIPHER)


heatmap_prox <- function(all_distances_df,DominantChain,genetic_distance,struct){
  
  library(colorRamp2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  data_heatmap = all_distances_df %>% group_by(cell_state) %>% summarize(Mean_lv = mean(lv, na.rm=TRUE))
  tmp = extract(data_heatmap, cell_state, c('row_state', 'col_state'), "([^ ]+) (.*)")
  
  data_heatmap$row_state= tmp$row_state
  data_heatmap$col_state= tmp$col_state
  mat = matrix(NA,6,6)
  rownames(mat)= c("BcellHIGH", "Cycling", "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory" )
  colnames(mat)= c("BcellHIGH", "Cycling", "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory" )
  
  
  for (l in c(1:nrow(data_heatmap))){
    
    for (i in c(1:6)){
      
      for (j in c(1:6)){
        
        indi=which(data_heatmap$row_state[l]==rownames(mat))
        indj=which(data_heatmap$col_state[l]==colnames(mat))
        mat[indi,indj]= data_heatmap$Mean_lv[l]
        mat[indj,indi]= data_heatmap$Mean_lv[l]
        
        
        
      }
      
    }
  }
  
  ind=which(rowSums(mat , na.rm = TRUE)==0)
  if (length(ind>0)){
    fix_mat= mat[,-c(ind)]
    fix_mat= fix_mat[-c(ind),]
  }else {
    fix_mat = mat
  }
  
  col_runif = colorRamp2(c(0,max(fix_mat[which(is.na(fix_mat)== FALSE)])/4,max(fix_mat[which(is.na(fix_mat)== FALSE)])/3 ,max(fix_mat[which(is.na(fix_mat)== FALSE)])/2,max(fix_mat[which(is.na(fix_mat)== FALSE)])), c("gray","yellow","orange","brown","red"))
  
  column_ha = HeatmapAnnotation(cell_status = rownames(fix_mat), col = list(cell_status = c("MetaboHIGH" = "#8DD3C7", "Inflammatory" = "#FFFFB3", "Cycling"= "#BEBADA","Neural" = "#FB8072","BcellHIGH" ="#80B1D3", "UPR_Secretory"= "#111999")))
  
  row_ha = rowAnnotation(cell_status = rownames(fix_mat), col = list(cell_status = c("MetaboHIGH" = "#8DD3C7", "Inflammatory" = "#FFFFB3", "Cycling"= "#BEBADA","Neural" = "#FB8072","BcellHIGH" ="#80B1D3", "UPR_Secretory"= "#111999")))
  pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_HEATMAP_CELLSTATE_averageLV_',struct,'.pdf'),width=10,height=10)
  print(Heatmap(fix_mat ,top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,show_row_names = FALSE, na_col = "black", cluster_columns=FALSE,
                cluster_rows=FALSE, cell_fun = function(j, i, x, y, width, height, fill) 
                {
                  grid.text(sprintf("%.1f", fix_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                }
  ))
  dev.off()
  
  
  write.csv(fix_mat,paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_HEATMAP_CELLSTATE_averageLV_',struct,'.csv'))
  
  
  return()
  
}       



alignment <- function(CDR3_data){
  
  seqs= DNAStringSet(CDR3_data)
  
  aligned <- AlignSeqs(seqs)
  
  return(aligned)  
  
}

alignment_score <- function(CDR3_data){
  dist1=NA
  dist1_matrix= matrix(NA, length(CDR3_data),length(CDR3_data))
  main_string_score_align=NA
  for (i in c(1:length(CDR3_data))){
    main_string= CDR3_data[i]
    
    for (j in c(1:length(CDR3_data))){
      dist1[j]= Biostrings::stringDist(AlignSeqs(DNAStringSet(c(main_string, CDR3_data[j]))))
      #dist[j]= Biostrings::stringDist(DNAStringSet(c(main_string, CDR3_data[j])))
    }
    main_string_score_align[i]=mean(dist1)
    dist1_matrix[i,]=dist1
    
  }
  res_alignment_score=list(main_string_score_align,dist1_matrix)
  names(res_alignment_score)=c("scores","dist1_matrix")
  return(res_alignment_score)
  
}



sample="HL12"
DominantChain="IGL"
model_name="F81"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample))
dir.create(paste0('no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name))

filtered_BCR_processed = read.csv(paste0('no_doublet_VDJfixed_clean_BCR_data_HRS_Bcells_',DominantChain,'.csv'))
#colnames(filtered_BCR_processed) = c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )
#filtered_BCR_processed$cell_id1=substr(filtered_BCR_processed$contig_id,1,16)
filtered_BCR_processed= filtered_BCR_processed %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count))  


### read contigs with 1 read (noise) and singlet
# singlet_low_confidential= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/*_singlet_IGL_IGK_read_thre1.txt'))[1])
# singlet= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/*_singlet_IGL_IGK_read_thre1.txt'))[2])
# doublet= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/*_doublet_IGL_IGK_read_thre1.txt')))
# 
pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))

singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
filtered_BCR_processed= filtered_BCR_processed[filtered_BCR_processed$cell_id1 %in% singlet,]

####reading the cell states##
new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_May2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)

##### filtering the cells to keep HLsample HRS and Bcells after removing VDJ and ; 
HLsample_HRS=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("HRS"),] 
HLsample_Bcells=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("Bcells"),] 
HLsample_Tcells=new_clusters[new_clusters$Patient==sample & new_clusters$MainCelltype %in% c("Tcells"),] 

write.table(paste0(HLsample_HRS$cell_id1,'-1'),'HLsample_HRS.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)
write.table(paste0(HLsample_Bcells$cell_id1,'-1'),'HLsample_Bcells.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)
write.table(paste0(HLsample_Tcells$cell_id1,'-1'),'HLsample_Tcells.txt',col.names=FALSE, row.names= FALSE, quote= FALSE)

###import the output file of data with fixed VDJ and removed doublets, call it filtered_BCR_processed # this file is free of doublets and fixed VDJ


filtered_BCR_processed=filtered_BCR_processed[filtered_BCR_processed$cell_id1 %in% HLsample_HRS$cell_id1,]# getting only HRS cells
ind=match(filtered_BCR_processed$cell_id1,HLsample_HRS$cell_id1)# adding the cell state
filtered_BCR_processed$new_cell_state=HLsample_HRS$SubtypeName[ind]
copy_filtered_BCR_processed=filtered_BCR_processed


##### keeping only this group of HRS cells ("Cycling",   "Inflammatory",      "Metabolic",    "Bcell-HIGH", "Neural-Glial")
filtered_BCR_processed=filtered_BCR_processed[filtered_BCR_processed$new_cell_state %in% c("MetaboHIGH", "Inflammatory", "Cycling","Neural" ,"BcellHIGH", "UPR_Secretory"),]

## reading the airr.tsv file output of TRUST4 to extract the full sequence
TRUST4_airr = read.table(Sys.glob('*_barcode_airr.tsv'), sep='\t', header=TRUE)
#TRUST4_airr$sequence_id=sub('_[^_]*$', '', TRUST4_airr$sequence_id)
filtered_BCR_processed_sequence_ind= match(filtered_BCR_processed$contig_id, TRUST4_airr$sequence_id)
filtered_BCR_processed$sequence = TRUST4_airr$sequence[filtered_BCR_processed_sequence_ind]

###### GENERATING THE MAXMUM LIKLIHOOD TREE ABSED ON ALIGNED SEQUENCE BASED ON SUBSTITUTION MODEL F81
# ind_notNA_sequence= which(is.na(filtered_BCR_processed$sequence)=="FALSE")
# 
# filtered_BCR_processed= filtered_BCR_processed[ind_notNA_sequence,]

# converting the charatcer to dna 
dna_sample=char2dna(filtered_BCR_processed$CDR3, simplify = FALSE) ## simplify= FALSE is the default
##### to align the sequence
class(dna_sample)
X11= align(dna_sample)
rownames(X11)<- filtered_BCR_processed$contig_id
## we need to change the phylo class
X111.phydat<-as.phyDat(X11)
HLsample_X11= X11

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

### maximum liklihood tree
set.seed(10)
tre.ini<-upgma(fix_dist.X111.phydat)
fit.ini <- pml(tre.ini, X111.phydat, k=4)#### pml tree
##optimizing the tree
fit <- optim.pml(fit.ini)

#### plonew_summarized_filtered_BCR_processeding the main structure of the tree without assignments
##
ape::write.tree(fit$tree,paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/','phylo_full_tree_Jan_2025'))
p <- ggtree(fit$tree, branch.length='branch.length',color="black", size=0.5, linetype=1, ladderize=TRUE) #+ geom_tiplab(size=2)#+ geom_tiplab(size=5) #+ geom_nodepoint(color="blue", size=3)


#p<- read.tree(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/','phylo_full_tree_Jan_2025'))
#p <- ggtree(p, branch.length='branch.length',color="black", size=0.5, linetype=1, ladderize=TRUE)#+ geom_nodepoint(color="blue", size=3) #geom_tiplab(size=5) 
### WE ARE ADDING THE ASSIGNMENTS

filtered_BCR_processed$status_singlet=NA
filtered_BCR_processed$status_singlet[filtered_BCR_processed$cell_id1 %in% singlet]="singlet"
#filtered_BCR_processed$status_singlet[filtered_BCR_processed$cell_id1 %in% singlet_low_confidential$V1]="singlet_low_confident"

filtered_BCR_processed$branch.length=p$data$branch.length[match(filtered_BCR_processed$contig_id,p$data$label)]     
filtered_BCR_processed$branch=p$data$branch[match(filtered_BCR_processed$contig_id,p$data$label)] 
filtered_BCR_processed$new_VDJ= paste0(filtered_BCR_processed$corrected_Vgene," ",filtered_BCR_processed$corrected_Dgene," ",filtered_BCR_processed$corrected_Jgene)



summarized_filtered_BCR_processed= filtered_BCR_processed %>% group_by(branch) %>% summarise(n = n(), contig_id,new_VDJ,malig_status,sequence,branch.length,branch,new_cell_state,CDR3,status_singlet)


### generating the tiplabels
TipLabels <- data.frame(tree.label = p$data$label, contig_id = p$data$label)
### adding the tiplabels to tree data
summarized_filtered_BCR_processed <- as.data.frame(summarized_filtered_BCR_processed)

summarized_filtered_BCR_processed <- full_join(x=summarized_filtered_BCR_processed, y=TipLabels, by= "contig_id")




### this line probebly will be removed
summarized_filtered_BCR_processed <- summarized_filtered_BCR_processed[is.na(summarized_filtered_BCR_processed$contig)==FALSE,]
rownames(summarized_filtered_BCR_processed) <- summarized_filtered_BCR_processed$contig

TipLabels.WithoutNAs <- TipLabels[is.na(TipLabels$tree.label)==FALSE,]### this line probably will be removed
rownames(TipLabels.WithoutNAs) <- TipLabels.WithoutNAs$contig
summarized_filtered_BCR_processed <-summarized_filtered_BCR_processed[rownames(TipLabels.WithoutNAs),]



clusters.nodeNetwork =summarized_filtered_BCR_processed#   filtered_BCR_processed
## creating an empty dataframe and add tree data
emptyRows <- data.frame(matrix(ncol = ncol(clusters.nodeNetwork), nrow = (nrow(TipLabels) - nrow(clusters.nodeNetwork))))
colnames(emptyRows) <- colnames(clusters.nodeNetwork)
clusters.nodeNetwork1 <- rbind(clusters.nodeNetwork, emptyRows)



### assigning the tips to tree
TipLabels$status <- clusters.nodeNetwork1$malig_status
TipLabels$id <-TipLabels$contig_id
TipLabels$vdj_status <- clusters.nodeNetwork1$new_VDJ
TipLabels$new_cell_state <- clusters.nodeNetwork1$new_cell_state

TipLabels$NodeSize <- clusters.nodeNetwork1$n
TipLabels$CDR3 <- clusters.nodeNetwork1$CDR3
TipLabels$status_singlet <- clusters.nodeNetwork1$status_singlet
p$data$status <- as.factor(TipLabels$status)# malignancy status
p$data$vdj_status <- as.factor(TipLabels$vdj_status) # VDJ states 
p$data$new_cell_state <- as.factor(TipLabels$new_cell_state) # cell states
p$data$NodeSize <- as.numeric(as.character(TipLabels$NodeSize)) # The number of cells related to the node associated with a clade
p$data$CDR3 <- (TipLabels$CDR3) 
p$data$status_singlet <- (TipLabels$status_singlet) 
p$data$id<-TipLabels$id



fasta_file=cbind(p$data$label,p$data$CDR3)
colnames(fasta_file)=c("id","cdr3")

write.csv(fasta_file,paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/','full_fasta_file_cdr3.csv'))
write.dna(X11, file = paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/','my_gene.fasta'), format = 'fasta')


####PLOTING THE TREE WITH ASSIGNMENTS
# scale_color_manual(values=c(malignant="red", pre_malignant="black",germline="grey")) 
r <- p + geom_text2(aes(label=id), hjust=-.3) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

s <- p + geom_tippoint(aes(col = vdj_status), size= 3, alpha=.75) # + ggtitle(paste0(Variables$TumorId, " - VDJ status"))
q <-  p + geom_text2(aes(label=id), hjust=-.3)
cdr3_tree <- p + geom_text2(aes(label=CDR3), hjust=-.3) #geom_tippoint(aes(col = cdr3), size= 3, alpha=.75)

brewer.pal(n=5, name = "Set3")
color.tree <- c("MetaboHIGH" = "#8DD3C7", "Inflammatory" = "#FFFFB3", "Cycling"= "#BEBADA","Neural" = "#FB8072","BcellHIGH" ="#80B1D3", "UPR_Secretory"= "#111999")
c <- p + geom_tippoint(aes(col = new_cell_state), size= 3, alpha=.75)+ scale_color_manual(values=color.tree)
status_singlet_gg <- p + geom_tippoint(aes(col = status_singlet), size= 3, alpha=.75) # + ggtitle(paste0(Variables$TumorId, " - VDJ status"))

pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_upgma_GGTREE_clone_HLsample_',DominantChain,'_VDJ.pdf'),width =10, height = 80)
plot(s)
dev.off()
pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_umpga_GGTREE_clone_HLsample_',DominantChain,'_status.pdf'),width =10, height = 80)
plot(r)
dev.off()
pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_upmga_GGTREE_clone_HLsample_',DominantChain,'_NEW_CELL_STATE.pdf'),width =10, height = 80)
plot(c)
dev.off()
pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_upmga_GGTREE_clone_HLsample_',DominantChain,'_contigID.pdf'),width =10, height = 80)
plot(q)
dev.off()

pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_upgma_GGTREE_clone_HLsample_',DominantChain,'_status_singlet.pdf'),width =10, height = 80)
plot(status_singlet_gg)
dev.off()


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_upgma_GGTREE_clone_HLsample_',DominantChain,'_CDR3.pdf'),width =10, height = 80)
plot(cdr3_tree)
dev.off()


cell_color= p + geom_tippoint(aes(col = new_cell_state), size= 8, alpha=.75) + scale_color_manual(values=color.tree) # +geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/',Sys.Date(),'_F81_mpl_opt_upmga_GGTREE_clone_HLsample_',DominantChain,'_COLLAPSED_COLORFULE_NEW_preMalignant.pdf'),width = 20, height = 20)
plot(cell_color)
dev.off()

################# COLORING CLADES FOR HLsample; For this purpose we extract the node of interest number from the tree; we have 3 interesting clades in this tree

###

#### removing the longer clade and regenerating the clusters
long_clade_data=data.frame(table(summarized_filtered_BCR_processed$CDR3))

long_clade_data=long_clade_data[order(long_clade_data$Freq, decreasing = TRUE),]

#long_clade= read.table('/Users/saramoein/Documents/new_run_HL_May2025/HL10/no_germline_phylo_HLsample_CDR3_IGL_modelF81/bs09_threshold50_gretaer5member/my_geneta_phylo_full_tree_Jan_2025_IGL_clusterPicks_cluster22_sequenceList.txt')
summarized_filtered_BCR_processed$cluster="NA"


long_clade_contigs= summarized_filtered_BCR_processed$contig_id[which(summarized_filtered_BCR_processed$CDR3==long_clade_data$Var1[1])]
cut_tree <- TreeTools::DropTipPhylo(as.phylo(r), long_clade_contigs, preorder = TRUE, check = TRUE)


##tried mukltiple threhsold "c( 0.05, 0.1, 0.3, 0.5)"
for (genetic_distance in c(0.1)){
  dir.create(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/phylopart',genetic_distance))
  ape::write.tree(cut_tree,paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/phylopart',genetic_distance,'/','phylo_full_tree_Jan_2025_excludeLONGclade_new'))
}
########
#RUNNING THE PHYLOPART ON UNIX COMMAN LINE

# phylopart0.5 % for genetic_distance in 0.1; do echo $genetic_distance; cd /Users/saramoein/Documents/new_run_HL_May2025/HL12/no_germline_phylo_HLsample_CDR3_IGL_modelF81/phylopart$genetic_distance; java -Xmx32G -jar /users/saramoein/downloads/PhyloPart_v2.1/PhyloPart_v2.1.jar phylo_full_tree_Jan_2025_excludeLONGclade_new  $genetic_distance -oOUT$genetic_distance.txt ; done;
#########
dir.create(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/phylopart',genetic_distance,'/proximity_cell_state/'))
genetic_distance = 0.1#0.05 , 0.1, 0.3, 0.5
other_clusters=read.table(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/phylopart',genetic_distance,'/','OUT',genetic_distance,'.txt'), sep=',',header=TRUE)

ind = match(summarized_filtered_BCR_processed$contig_id, other_clusters$leafname)
summarized_filtered_BCR_processed$cluster=  other_clusters$clustername[ind]
summarized_filtered_BCR_processed$cluster[which(summarized_filtered_BCR_processed$CDR3==long_clade_data$Var1[1])]="LONG_cluster"

clusters= unique(summarized_filtered_BCR_processed$cluster)


library(Polychrome)
P36 = createPalette(50,  c("#ff0000", "#00ff00", "#0000ff", "#0ff0ff"))
k=0
r_copy= r
summarized_filtered_BCR_processed$col = NA
summarized_filtered_BCR_processed$cdr3_diversity_per_clade = NA


for (c in clusters){
  print(c)
  k = k +1
  #node_number = getMRCA(as.phylo(r), summarized_filtered_BCR_processed$contig_id[which(summarized_filtered_BCR_processed$cluster == c)])
  ind_data = which(summarized_filtered_BCR_processed$cluster==c)
  data= summarized_filtered_BCR_processed[ind_data,]
  summarized_filtered_BCR_processed$col[ind_data]= P36[k]
  
  if (!c %in% c("LONG_cluster","0")){
    summarized_filtered_BCR_processed$cdr3_diversity_per_clade[ind_data]= alignment_score(data$CDR3)$scores
  }  
  #r_copy= r_copy  + geom_hilight(node= node_number, fill= P36[k] ) 
  
  colors_cell_state<- color.tree
  library(TSP)
  rgb <- col2rgb(colors_cell_state)
  tsp <- as.TSP(dist(t(rgb)))
  sol <- solve_TSP(tsp, control = list(repetitions = 1e3))
  sorted_cols_cell_state=colors_cell_state[sol]
  
  
  df =data.frame(cell_state= c("BcellHIGH", "Cycling" , "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory"),
                 value = c(0,0,0,0,0,0),
                 color= c("#80B1D3" , "#BEBADA" ,  "#FFFFB3" ,"#8DD3C7","#FB8072" ,"#111999"))
  
  con= 0 
  pie_colors_cell_state=vector()
  for (cs in c("BcellHIGH", "Cycling" , "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory")){
    con = con +1
    ind1= which(data$new_cell_state == cs)
    
    df$value[con] = length(ind1)
    
  }
  
  ind_df = which(df$value != 0 )
  df = df[ind_df, ]
  pie_colors_cell_state <- df$color
  
  pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_pie_cell_state_cluster',c,'.pdf'))
  pie(df$value , df$cell_state ,col=pie_colors_cell_state)
  dev.off()
  
  pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_pie_VDJ_cluster',c,'.pdf'))
  pie(table(data$new_VDJ), cex=0.7)
  dev.off()
  
  
  
  
  
  
  ### start
  # length_cluster= length(which(summarized_filtered_BCR_processed$cluster==c))
  # 
  # if (length_cluster>=5 & c!=0 ){
  #   data_BcellHIGH= summarized_filtered_BCR_processed[which(summarized_filtered_BCR_processed$new_cell_state == "Bcell-HIGH" & summarized_filtered_BCR_processed$cluster==c & summarized_filtered_BCR_processed$cluster!= 0), ]
  #   data_Cycling= summarized_filtered_BCR_processed[which(summarized_filtered_BCR_processed$new_cell_state == "Cycling" & summarized_filtered_BCR_processed$cluster==c & summarized_filtered_BCR_processed$cluster!= 0),]
  #   data_Neural= summarized_filtered_BCR_processed[which(summarized_filtered_BCR_processed$new_cell_state == "Neural-Glial" & summarized_filtered_BCR_processed$cluster==c & summarized_filtered_BCR_processed$cluster!= 0),]
  #   data_inflammatory= summarized_filtered_BCR_processed[which(summarized_filtered_BCR_processed$new_cell_state == "Inflammatory" & summarized_filtered_BCR_processed$cluster==c & summarized_filtered_BCR_processed$cluster!= 0),]
  #   data_MetaboHIGH= summarized_filtered_BCR_processed[which(summarized_filtered_BCR_processed$new_cell_state == "Metabolic" & summarized_filtered_BCR_processed$cluster==c & summarized_filtered_BCR_processed$cluster!= 0) ,]
  #   
  #   
  #   
  #   data_states=list();
  #   data_states[[1]]= data_BcellHIGH$CDR3
  #   data_states[[2]]= data_Cycling$CDR3
  #   data_states[[3]]= data_Neural$CDR3
  #   data_states[[4]]= data_inflammatory$CDR3
  #   data_states[[5]]= data_MetaboHIGH$CDR3
  #   
  #   length_data_states= length(which(sapply( data_states, function(x) length(x) != 0 )==TRUE))
  #   if (length_data_states>1){
  #     all_distances = list()
  #     counter = 0 
  #     colnames_vector=vector()
  #     mat1 = matrix(NA, 500000, 10 )
  #     for (i in c(1:5)){
  #       print(i)
  #       for (j in c(1:5)){
  #         print(j)
  #         if (i<j){
  #           
  #           states = c("Bcell-HIGH", "Cycling" , "Neural-Glial",  "Inflammatory","Metabolic") 
  #           a = data.frame(data_states[[i]])
  #           b = data.frame(data_states[[j]])
  #           colnames(a)="CDR3"
  #           colnames(b)="CDR3"
  #           if (nrow(a > 0) & nrow(b > 0)){
  #             counter = counter + 1
  #             dist1= cbind(expand.grid(a = a$CDR3, b = b$CDR3), lv = c(adist(a$CDR3, b$CDR3)))
  #             cell_states = paste0(states[i], " ",states[j])
  #             lv_distance = dist1
  #             
  #             all_distances[[counter]] = list(cell_states,lv_distance$lv)
  #             
  #             # colnames(mat1)[con]= paste0(states[i], " ",states[j])
  #             mat1[1:nrow(dist1),counter] = dist1$lv
  #             colnames_vector = rbind(colnames_vector , paste0(states[i], " ",states[j]))
  #           }
  #         }
  #       }
  #     }
  #     
  #     mat1=data.frame(mat1)
  #     
  #     colnames(mat1) = colnames_vector
  #     
  #     pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_F81_mpl_opt_umpga_GGTREE_clone_HLsample_',DominantChain,'_cellstarte_DISTANCE',c,'.pdf'),width=25)
  #     boxplot(as.matrix(mat1[,1:10]), las=1, cex.axis=0.8)
  #     dev.off()
  #   } 
  # }
  # 
}


new_summarized_filtered_BCR_processed= summarized_filtered_BCR_processed[which(!is.na(summarized_filtered_BCR_processed$CDR3==TRUE)),]

new_summarized_filtered_BCR_processed=new_summarized_filtered_BCR_processed[which(new_summarized_filtered_BCR_processed$cluster!="NA"),]

new_summarized_filtered_BCR_processed= new_summarized_filtered_BCR_processed[which(new_summarized_filtered_BCR_processed$cluster !=0),]




p$data$cluster=NA
p$data$cluster_col = NA
ind = match(p$data$id,summarized_filtered_BCR_processed$contig_id)
p$data$cluster = summarized_filtered_BCR_processed$cluster[ind]
p$data$cluster_col = summarized_filtered_BCR_processed$col[ind]

names(P36) = unique(summarized_filtered_BCR_processed$cluster)
p$data$cluster_col[which(p$data$cluster==0)] = NA
p$data$cluster[which(p$data$cluster==0)] = NA
cluster_color <- p + geom_tippoint(aes(col = cluster), size= 4, alpha=.75) + scale_color_manual(values=P36) + geom_text2(aes(label=cluster, subset=isTip), hjust=-.1) 


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/phylopart',genetic_distance,'/',Sys.Date(),'_F81_mpl_opt_umpga_GGTREE_clone_HLsample_',DominantChain,'COLORcluster.pdf'),width =20, height = 80) 
plot(cluster_color)
dev.off()






pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_pie_VDJ_allHRS.pdf'))
pie(table(new_summarized_filtered_BCR_processed$new_VDJ), cex=0.7)
dev.off()



###### GENERATING THE BOX PLOTS BASED ON BRANCH LENGTH FOR ALL CELL STATES PER CLADES

pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_timing.pdf'))  
print(ggplot(new_summarized_filtered_BCR_processed, aes(x=cluster, y=round(as.numeric(branch),4), fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree) +
        facet_wrap(~cluster, scale="free"))
dev.off()



new_summarized_filtered_BCR_processed_excludeLONGERcluster= new_summarized_filtered_BCR_processed[which(new_summarized_filtered_BCR_processed$cluster!= "LONG_cluster" & new_summarized_filtered_BCR_processed$cluster!=0),]


new_summarized_filtered_BCR_processed_excludeLONGERcluster$allCluster = "exclude_longer_cluster"
pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_ExcludeLongerCluster_timing.pdf'))  
print(ggplot(new_summarized_filtered_BCR_processed_excludeLONGERcluster, aes(x=allCluster, y=round(as.numeric(branch),4), fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree))

dev.off()


new_summarized_filtered_BCR_processed_LONGERcluster=new_summarized_filtered_BCR_processed[which(new_summarized_filtered_BCR_processed$cluster=="LONG_cluster"),]


data_timing= cbind(new_summarized_filtered_BCR_processed_excludeLONGERcluster$branch,new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state,new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster)
tmp= cbind(new_summarized_filtered_BCR_processed_LONGERcluster$branch, new_summarized_filtered_BCR_processed_LONGERcluster$new_cell_state,new_summarized_filtered_BCR_processed_LONGERcluster$cluster)
data_timing=rbind(data_timing,tmp)
data_timing=data.frame(data_timing)
colnames(data_timing)=c("branch","cell_state","cluster")
data_timing$new_cluster=ifelse(data_timing$cluster!="LONG_cluster","not_long_cluster","LONG_cluster")


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_LongerCluster_BRANCHlength_FULLPLOT.pdf'))  
print(ggplot(data_timing, aes(x=cell_state, y=as.numeric(branch), fill= new_cluster)) + geom_boxplot() + scale_fill_manual(values = c(not_long_cluster = "#8DD3C7", LONG_cluster = "#FB8072")))
# scale_fill_manual(values = c(Metabolic = "#8DD3C7", Inflammatory = "#FFFFB3", Cycling= "#BEBADA",`Neural-Glial` = "#FB8072",`Bcell-HIGH` ="#80B1D3")))
dev.off()






new_summarized_filtered_BCR_processed=new_summarized_filtered_BCR_processed[order(new_summarized_filtered_BCR_processed$col),]
#data_to_align=new_summarized_filtered_BCR_processed$CDR3
gc()
new_summarized_filtered_BCR_processed$allCluster = ifelse(new_summarized_filtered_BCR_processed$cluster!="LONG_cluster","not_long_cluster","LONG_cluster")
tmp_align=rbind(new_summarized_filtered_BCR_processed[which(new_summarized_filtered_BCR_processed$allCluster=="not_long_cluster"),],new_summarized_filtered_BCR_processed[which(new_summarized_filtered_BCR_processed$allCluster=="LONG_cluster")[1],])
tmp_align$cdr3_diversity = alignment_score(tmp_align$CDR3)$scores
gc()
new_summarized_filtered_BCR_processed$cdr3_diversity=NA
new_summarized_filtered_BCR_processed$cdr3_diversity[which(new_summarized_filtered_BCR_processed$allCluster=="LONG_cluster")]= tmp_align$cdr3_diversity[nrow(tmp_align)]
new_summarized_filtered_BCR_processed$cdr3_diversity[which(new_summarized_filtered_BCR_processed$allCluster=="not_long_cluster")]=tmp_align$cdr3_diversity[1:(nrow(tmp_align)-1)]

gc()
tmp_align$aligned_cdr3 = alignment(tmp_align$CDR3)
gc()
new_summarized_filtered_BCR_processed$aligned_cdr3=NA
new_summarized_filtered_BCR_processed$aligned_cdr3[which(new_summarized_filtered_BCR_processed$allCluster=="LONG_cluster")]= tmp_align$aligned_cdr3[nrow(tmp_align)]
new_summarized_filtered_BCR_processed$aligned_cdr3[which(new_summarized_filtered_BCR_processed$allCluster=="not_long_cluster")]=tmp_align$aligned_cdr3[1:(nrow(tmp_align)-1)]
aligned_data=DNAStringSet(new_summarized_filtered_BCR_processed$aligned_cdr3)
names(aligned_data)= new_summarized_filtered_BCR_processed$cluster
BrowseSeqs(aligned_data,htmlFile = paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_alignment_fulltree.html'))

gc()
new_summarized_filtered_BCR_processed_excludeLONGERcluster=new_summarized_filtered_BCR_processed_excludeLONGERcluster[order(new_summarized_filtered_BCR_processed_excludeLONGERcluster$col),]
data_to_align= new_summarized_filtered_BCR_processed_excludeLONGERcluster$CDR3
names(data_to_align)= new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster
new_summarized_filtered_BCR_processed_excludeLONGERcluster$allCluster = "exclude_longer_cluster"
gc()
new_summarized_filtered_BCR_processed_excludeLONGERcluster$cdr3_diversity = alignment_score(new_summarized_filtered_BCR_processed_excludeLONGERcluster$CDR3)$scores
gc()
aligned_data_excludedlong = alignment(data_to_align)
new_summarized_filtered_BCR_processed_excludeLONGERcluster$aligned_cdr3=aligned_data_excludedlong
BrowseSeqs(aligned_data_excludedlong,htmlFile = paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_alignment_excludeLongTree.html'))


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_fulltree_cdr3Diversity.pdf'))  
print(ggplot(new_summarized_filtered_BCR_processed, aes(x=cluster, y=cdr3_diversity_per_clade, fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree) +
        facet_wrap(~cluster, scale="free"))
dev.off()




pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_NEW_colorClades_CELLSTATE_barplot_phylo_ExcludeLongerCluster_cdr3Diversity.pdf'))  
print(ggplot(new_summarized_filtered_BCR_processed_excludeLONGERcluster, aes(x=allCluster, y=cdr3_diversity, fill= new_cell_state)) + geom_boxplot() +
        scale_fill_manual(values = color.tree))

dev.off()



# gc()
# dna_cdr3=char2dna(new_summarized_filtered_BCR_processed_excludeLONGERcluster$CDR3, simplify = FALSE) ## simplify= FALSE is the default
# ##### to align the sequence
# class(dna_cdr3)
# X11= align(dna_cdr3)
# rownames(X11)<- new_summarized_filtered_BCR_processed_excludeLONGERcluster$contig_id
# ## we need to change the phylo class
# X111.phydat<-as.phyDat(X11)
# HLsample_X11= X11
# #genrating the distance matrix from aligned sequence
# dist.X111.phydat = dist.dna(X11, model = model_name, variance = FALSE,
#                             gamma = FALSE, pairwise.deletion = TRUE,
#                             base.freq = NULL, as.matrix = TRUE)# , indel=TRUE
# giveNAs = which(is.nan(as.matrix(dist.X111.phydat)),arr.ind=TRUE)
# 
# if (nrow(giveNAs)>0){
#   fix_dist.X111.phydat = dist.X111.phydat[-c(giveNAs[,1]),-c(giveNAs[,2])]
# }else if (nrow(giveNAs)==0){
#   
#   fix_dist.X111.phydat=dist.X111.phydat
#   
# }
# gc()
# dist1= as.matrix(fix_dist.X111.phydat)
# row_dist1= rownames(dist1)
# 
# 
# ind=match(rownames(dist1),new_summarized_filtered_BCR_processed_excludeLONGERcluster$contig_id)
# dist_dat=new_summarized_filtered_BCR_processed_excludeLONGERcluster[ind,]
# col_runif = colorRamp2(c(0,max(dist1)), c("gray","red"))
# 
# #column_ha = HeatmapAnnotation(cell_status =  new_summarized_filtered_BCR_processed$new_cell_state, clade =  new_summarized_filtered_BCR_processed$allCluster, col = list(cell_status = c("Metabolic" = "green", "Inflammatory" = "yellow", "Cycling"= "purple","Neural-Glial" = "red","Bcell-HIGH" ="blue")))
# column_ha = HeatmapAnnotation(cell_status =  dist_dat$new_cell_state, clade =  dist_dat$allCluster, col = list(cell_status = c("Metabolic" = "#8DD3C7", "Inflammatory" = "#FFFFB3", "Cycling"= "#BEBADA","Neural-Glial" = "#FB8072","Bcell-HIGH" ="#80B1D3")))
# 
# row_ha = rowAnnotation(cell_status =  dist_dat$new_cell_state, clade =  dist_dat$allCluster, col = list(cell_status = c("Metabolic" = "#8DD3C7", "Inflammatory" = "#FFFFB3", "Cycling"= "#BEBADA","Neural-Glial" = "#FB8072","Bcell-HIGH" ="#80B1D3")))
# pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81/','/phylopart',genetic_distance,'/',Sys.Date(),'_HEATMAP_CELLSTATE_barplot_phylo_excludeLonger_Cluster.pdf'),width=10,height=10)
# print(Heatmap(dist1,col= col_runif, top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,show_row_names = FALSE))
# dev.off()
# 
# 
# 



data_BcellHIGH= new_summarized_filtered_BCR_processed_excludeLONGERcluster[which((new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state == "BcellHIGH") & (new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster!=0) ), ]
data_Cycling= new_summarized_filtered_BCR_processed_excludeLONGERcluster[which((new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state == "Cycling") & (new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster!=0) ),]
data_Neural= new_summarized_filtered_BCR_processed_excludeLONGERcluster[which((new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state == "Neural") & (new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster!=0)),]
data_inflammatory= new_summarized_filtered_BCR_processed_excludeLONGERcluster[which((new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state == "Inflammatory") & (new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster!=0) ),]
data_MetaboHIGH= new_summarized_filtered_BCR_processed_excludeLONGERcluster[which((new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state == "MetaboHIGH") & (new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster!=0) ) ,]
data_UPR_Secretory= new_summarized_filtered_BCR_processed_excludeLONGERcluster[which((new_summarized_filtered_BCR_processed_excludeLONGERcluster$new_cell_state == "UPR_Secretory") & (new_summarized_filtered_BCR_processed_excludeLONGERcluster$cluster!=0) ) ,]
data_states=list();
data_states[[1]]= data_BcellHIGH
data_states[[2]]= data_Cycling
data_states[[3]]= data_Neural
data_states[[4]]= data_inflammatory
data_states[[5]]= data_MetaboHIGH
data_states[[6]]= data_UPR_Secretory

length_data_states= length(which(sapply( data_states, function(x) length(x) >1)==TRUE))


gc()
if (length_data_states>1){     
  all_distances_df=data.frame()
  all_cdr3divergence_df=data.frame()
  
  all_distances = list()
  con = 0 
  colnames_vector=vector()
  
  for (i in c(1:6)){
    print(i)
    for (j in c(1:6)){
      print(j)
      if (i<=j){
        
        states = c("BcellHIGH", "Cycling" , "Neural",  "Inflammatory","MetaboHIGH","UPR_Secretory") 
        a = data.frame(data_states[[i]])
        b = data.frame(data_states[[j]])
        # colnames(a)="CDR3_aligned"
        # colnames(b)="CDR3_aligned"
        if (nrow(a)> 1 & nrow(b) > 1){
          con = con + 1
          dist1= cbind(expand.grid(a = a$aligned_cdr3, b = b$aligned_cdr3), lv = c(adist(a$aligned_cdr3, b$aligned_cdr3)))
          cell_states = paste0(states[i], " ",states[j])
          lv_distance = dist1
          
          all_distances[[con]] = list(cell_states,lv_distance$lv)
          
          tmp=data.frame(lv_distance$lv)
          colnames(tmp)= "lv"
          tmp$cell_state[1:length(lv_distance$lv)]=rep(cell_states,length(lv_distance$lv))
          
          all_distances_df= rbind(tmp,all_distances_df)
          cell_similar_ind= which(!(all_distances_df$cell_state %in% c("Cycling Cycling","BcellHIGH BcellHIGH","Neural Neural","MetaboHIGH MetaboHIGH","Inflammatory Inflammatory", "UPR_Secretory UPR_Secretory") & (all_distances_df$lv==0)))
          all_distances_df=all_distances_df[cell_similar_ind,]
          
          
          dvg= c(a$cdr3_diversity,b$cdr3_diversity)
          tmp_dvg= data.frame(dvg)
          colnames(tmp)= "prox_cdr3divergence"
          tmp_dvg$cell_state[1:length(tmp_dvg)]=rep(cell_states,length(tmp_dvg))
          all_cdr3divergence_df= rbind(tmp_dvg,all_cdr3divergence_df)
          
          # colnames(mat1)[con]= paste0(states[i], " ",states[j])
          #mat1[1:nrow(dist1),con] = dist1$lv
          colnames_vector = rbind(colnames_vector , paste0(states[i], " ",states[j]))
          
        }
      }
    }
    
  }
  
  
  
  
  library(rstatix)
  library(ggpubr)
  
  res.aov <- all_distances_df %>% anova_test(lv ~ cell_state)
  pwc <- all_distances_df %>%
    pairwise_t_test(lv ~ cell_state, p.adjust.method = "bonferroni")
  
  write.csv(data.frame(pwc),paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_lvDistance_excludeLongCluster_anova.csv'))
  
  
  pwc <- pwc %>% add_xy_position(x = "cell_state")
  group.order <- c(unique(pwc$group1), unique(pwc$group2)[!(unique(pwc$group2) %in% unique(pwc$group1))])
  
  pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_lvDistance_excludeLongCluster_anova.pdf'))
  g <- (ggboxplot(all_distances_df, x = "cell_state", y = "lv",fill="pink") +
          stat_pvalue_manual(pwc, label = "p.signif", tip.length = 0, step.increase = 0.2) +
          scale_x_discrete(limits = group.order) +
          labs(
            caption = get_pwc_label(pwc)
          ))+ rotate_x_text(90)
  print(g+  font("title", size = 14, color = "red", face = "bold.italic")+
          font("subtitle", size = 10, color = "orange")+
          font("caption", size = 10, color = "orange")+
          font("xlab", size = 12, color = "blue")+
          font("ylab", size = 12, color = "blue")+
          font("xy.text", size = 7, color = "black", face = "bold")
  )
  dev.off()
  
  res.aov.cdr3dvg <- all_cdr3divergence_df %>% anova_test(dvg ~ cell_state)
  pwc.dvg <- all_cdr3divergence_df %>%
    pairwise_t_test(dvg ~ cell_state, p.adjust.method = "bonferroni")
  
  write.csv(data.frame(pwc.dvg),paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_CDR3divergence_excludeLongCluster_anova.csv'))
  
  pwc.dvg <- pwc.dvg %>% add_xy_position(x = "cell_state")
  group.order <- c(unique(pwc.dvg$group1), unique(pwc.dvg$group2)[!(unique(pwc.dvg$group2) %in% unique(pwc.dvg$group1))])
  pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_CDR3divergence_excludeLongCluster_anova.pdf'))
  gg<- ggboxplot(all_cdr3divergence_df, x = "cell_state", y = "dvg",fill="brown") +
    stat_pvalue_manual(pwc.dvg, label = "p.signif", tip.length = 0, step.increase = 0.2) +
    scale_x_discrete(limits = group.order) +
    labs(
      #subtitle = get_test_label(res.aov.cdr3dvg, detailed = TRUE),
      caption = get_pwc_label(pwc.dvg)
    )+ rotate_x_text(90)
  print(gg +  font("title", size = 14, color = "red", face = "bold.italic")+
          font("subtitle", size = 10, color = "orange")+
          font("caption", size = 10, color = "orange")+
          font("xlab", size = 12, color = "blue")+
          font("ylab", size = 12, color = "blue")+
          font("xy.text", size = 7, color = "black", face = "bold")
  )
  
  dev.off()
  
  
  struct="exclude_longer_clade"
  heatmap_prox(all_distances_df,DominantChain,genetic_distance,struct)
  
  
}


# data_BcellHIGH= new_summarized_filtered_BCR_processed[which((new_summarized_filtered_BCR_processed$new_cell_state == "Bcell-HIGH") & (new_summarized_filtered_BCR_processed$cluster!=0) ), ]
# data_Cycling= new_summarized_filtered_BCR_processed[which((new_summarized_filtered_BCR_processed$new_cell_state == "Cycling") & (new_summarized_filtered_BCR_processed$cluster!=0) ),]
# data_Neural= new_summarized_filtered_BCR_processed[which((new_summarized_filtered_BCR_processed$new_cell_state == "Neural-Glial") & (new_summarized_filtered_BCR_processed$cluster!=0)),]
# data_inflammatory= new_summarized_filtered_BCR_processed[which((new_summarized_filtered_BCR_processed$new_cell_state == "Inflammatory") & (new_summarized_filtered_BCR_processed$cluster!=0) ),]
# data_MetaboHIGH= new_summarized_filtered_BCR_processed[which((new_summarized_filtered_BCR_processed$new_cell_state == "Metabolic") & (new_summarized_filtered_BCR_processed$cluster!=0) ) ,]
# data_states=list();
# data_states[[1]]= data_BcellHIGH
# data_states[[2]]= data_Cycling
# data_states[[3]]= data_Neural
# data_states[[4]]= data_inflammatory
# data_states[[5]]= data_MetaboHIGH
# length_data_states= length(which(sapply( data_states, function(x) length(x) != 0 )==TRUE))
# 
# 
# gc()
# if (length_data_states>1){     
#   all_distances_df=data.frame()
#   all_cdr3divergence_df=data.frame()
#   
#   all_distances = list()
#   con = 0 
#   colnames_vector=vector()
#   mat1 = matrix(NA, 500000, 10 )
#   for (i in c(1:5)){
#     print(i)
#     for (j in c(1:5)){
#       print(j)
#       if (i<=j){
#         
#         states = c("Bcell-HIGH", "Cycling" , "Neural-Glial",  "Inflammatory","Metabolic") 
#         a = data.frame(data_states[[i]])
#         b = data.frame(data_states[[j]])
#         # colnames(a)="CDR3_aligned"
#         # colnames(b)="CDR3_aligned"
#         if (nrow(a)> 0 & nrow(b) > 0){
#           con = con + 1
#           dist1= cbind(expand.grid(a = a$aligned_cdr3, b = b$aligned_cdr3), lv = c(adist(a$aligned_cdr3, b$aligned_cdr3)))
#           cell_states = paste0(states[i], " ",states[j])
#           lv_distance = dist1
#           
#           all_distances[[con]] = list(cell_states,lv_distance$lv)
#           
#           tmp=data.frame(lv_distance$lv)
#           colnames(tmp)= "lv"
#           tmp$cell_state[1:length(lv_distance$lv)]=rep(cell_states,length(lv_distance$lv))
#           
#           all_distances_df= rbind(tmp,all_distances_df)
#           
#           cell_similar_ind= which(!(all_distances_df$cell_state %in% c("Cycling Cycling","Bcell-HIGH Bcell-HIGH","Neural-Glial Neural-Glial","Metabolic Metabolic","Inflammatory Inflammatory") & (all_distances_df$lv==0)))
#           all_distances_df=all_distances_df[cell_similar_ind,]
#           
#           dvg= c(a$cdr3_diversity,b$cdr3_diversity)
#           tmp_dvg= data.frame(dvg)
#           colnames(tmp)= "prox_cdr3divergence"
#           tmp_dvg$cell_state[1:length(tmp_dvg)]=rep(cell_states,length(tmp_dvg))
#           all_cdr3divergence_df= rbind(tmp_dvg,all_cdr3divergence_df)
#           
#           # colnames(mat1)[con]= paste0(states[i], " ",states[j])
#           #mat1[1:nrow(dist1),con] = dist1$lv
#           colnames_vector = rbind(colnames_vector , paste0(states[i], " ",states[j]))
#           
#         }
#       }
#     }
#     
#   }
#   
#   
#   
#   
#   library(rstatix)
#   library(ggpubr)
#   
#   res.aov <- all_distances_df %>% anova_test(lv ~ cell_state)
#   pwc <- all_distances_df %>%
#     pairwise_t_test(lv ~ cell_state, p.adjust.method = "bonferroni")
#   
#   write.csv(data.frame(pwc),paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_lvDistance_fulltree_anova.csv'))
#   
#   
#   pwc <- pwc %>% add_xy_position(x = "cell_state")
#   group.order <- c(unique(pwc$group1), unique(pwc$group2)[!(unique(pwc$group2) %in% unique(pwc$group1))])
#   pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_lvDistance_fulltree_anova.pdf'))
#   g <- ggboxplot(all_distances_df, x = "cell_state", y = "lv",fill="gray") +
#     stat_pvalue_manual(pwc, label = "p.signif", tip.length = 0, step.increase = 0.2) +
#     scale_x_discrete(limits = group.order) +
#     labs(
#       #subtitle = get_test_label(res.aov, detailed = TRUE),
#       caption = get_pwc_label(pwc)
#     )+ rotate_x_text(90)
#   print(g+  font("title", size = 14, color = "red", face = "bold.italic")+
#           font("subtitle", size = 10, color = "orange")+
#           font("caption", size = 10, color = "orange")+
#           font("xlab", size = 12, color = "blue")+
#           font("ylab", size = 12, color = "blue")+
#           font("xy.text", size = 7, color = "black", face = "bold")
#   )
#   dev.off()
#   
#   res.aov.cdr3dvg <- all_cdr3divergence_df %>% anova_test(dvg ~ cell_state)
#   pwc.dvg <- all_cdr3divergence_df %>%
#     pairwise_t_test(dvg ~ cell_state, p.adjust.method = "bonferroni")
#   
#   write.csv(data.frame(pwc.dvg),paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_CDR3divergence_fulltree_anova.csv'))
#   
#   pwc.dvg <- pwc.dvg %>% add_xy_position(x = "cell_state")
#   group.order <- c(unique(pwc.dvg$group1), unique(pwc.dvg$group2)[!(unique(pwc.dvg$group2) %in% unique(pwc.dvg$group1))])
#   pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_CDR3divergence_fulltree_anova.pdf'))
#   gg<- ggboxplot(all_cdr3divergence_df, x = "cell_state", y = "dvg",fill="blue") +
#     stat_pvalue_manual(pwc.dvg, label = "p.signif", tip.length = 0, step.increase = 0.2) +
#     scale_x_discrete(limits = group.order) +
#     labs(
#       #subtitle = get_test_label(res.aov.cdr3dvg, detailed = TRUE),
#       caption = get_pwc_label(pwc.dvg)
#     )+ rotate_x_text(90)
#   
#   print(gg+  font("title", size = 14, color = "red", face = "bold.italic")+
#           font("subtitle", size = 10, color = "orange")+
#           font("caption", size = 10, color = "orange")+
#           font("xlab", size = 12, color = "blue")+
#           font("ylab", size = 12, color = "blue")+
#           font("xy.text", size = 7, color = "black", face = "bold")
#   )
#   
#   
#   dev.off()
#   struct="full_tree"
#   heatmap_prox(all_distances_df,DominantChain,genetic_distance,struct)
#   
# }
# 






data= new_summarized_filtered_BCR_processed_excludeLONGERcluster


df =data.frame(cell_state= c("BcellHIGH", "Cycling" , "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory"),
               value = c(0,0,0,0,0,0),
               color= c("#80B1D3" , "#BEBADA" ,  "#FFFFB3" ,"#8DD3C7","#FB8072" ,"#111999"))
con= 0 
pie_colors_cell_state=vector()
for (cs in c("BcellHIGH", "Cycling" , "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory")){
  con = con +1
  ind1= which(data$new_cell_state == cs)
  
  df$value[con] = length(ind1)
  
}

ind_df = which(df$value != 0 )
df = df[ind_df, ]
pie_colors_cell_state <- df$color

pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_PIE_',DominantChain,'_cellState_excludeLONGERcluster.pdf'))
pie(df$value , df$cell_state ,col=pie_colors_cell_state)
dev.off()


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_PIE_',DominantChain,'_VDJ_excludeLONGERcluster.pdf'))
pie(table(data$new_VDJ), cex=0.7)
dev.off()






data= new_summarized_filtered_BCR_processed

df =data.frame(cell_state= c("BcellHIGH", "Cycling" , "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory"),
               value = c(0,0,0,0,0,0),
               color= c("#80B1D3" , "#BEBADA" ,  "#FFFFB3" ,"#8DD3C7","#FB8072" ,"#111999"))

con= 0 
pie_colors_cell_state=vector()
for (cs in c("BcellHIGH", "Cycling" , "Inflammatory" , "MetaboHIGH" , "Neural", "UPR_Secretory")){
  con = con +1
  ind1= which(data$new_cell_state == cs)
  
  df$value[con] = length(ind1)
  
}

ind_df = which(df$value != 0 )
df = df[ind_df, ]
pie_colors_cell_state <- df$color

pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_PIE_',DominantChain,'_cellState_allTREE.pdf'))
pie(df$value , df$cell_state ,col=pie_colors_cell_state)
dev.off()


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/',Sys.Date(),'_PIE_',DominantChain,'_VDJ_allTREE.pdf'))
pie(table(data$new_VDJ), cex=0.7)
dev.off()

write.csv(new_summarized_filtered_BCR_processed,paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_model',model_name,'/','summarized_filtered_BCR_processed.csv'))


pdf(paste0('./no_germline_phylo_HLsample_CDR3_',DominantChain,'_modelF81','/phylopart',genetic_distance,'/proximity_cell_state/',Sys.Date(),'_cellstate_proximity_CDR3divergence_excludeLongCluster_BARPLOT.pdf'))
ggplot(all_distances_df, aes(x = cell_state, y = lv)) +
  geom_boxplot() +
  labs(title = "Bar Plot cell states distance",
       x = "cell state",
       y = "aligned CDR3 distance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

