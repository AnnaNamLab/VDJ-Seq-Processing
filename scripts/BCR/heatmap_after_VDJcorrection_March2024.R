# we first run trust4, then getting cdr3_raw.out file
### doublets and noise cells are removed. Only Bcells and HRS cells are included.
###setting the directory

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

args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
library(optparse)
option_list <- list(
  make_option(c("-s", "--sample"), type = "character", help = "Sample name"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to cell metadata CSV file"),
  make_option(c("-i", "--input"), type = "character", help = "path to output of trust4"),
  make_option(c("-d", "--DominantChain"), type = "character", help = "dominant chian"),
  make_option(c("-l", "--light_doublet"), type = "character", help = "output from doublets on light chain (e.g., HL1_Rawdata_IGK_IGL.csv)"),
  make_option(c("-v", "--heavy_doublet"), type = "character", help = "output from doublets on light chain (e.g., HL1_Rawdata_IGH.csv)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$input) || is.null(opt$metadata) || is.null(opt$sample) || is.null(opt$DominantChain) || is.null(opt$light_doublet) ||is.null(opt$heavy_doublet)) {
  stop("Please provide --input, --metadata, --sample, and --file arguments.")
}



sample <- opt$sample
metadata_file <- opt$metadata
input_directory <- opt$input
DominantChain <- opt$DominantChain
light_chain_cellStatus <- opt$light_doublet
heavy_chain_cellStatus <- opt$heavy_doublet


# -----------------------
# Set working directory
# -----------------------
dir.create(input_directory, showWarnings = FALSE)
setwd(input_directory)


#for (sample in c("HL10")){
      #DominantChain="IGK"
      #setwd(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample))
      
      
      ### reading the cell types from GEX
      #new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_May2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
      new_clusters= read.csv(metadata_file)
      new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)
      
      
      
      #### extracting the Bcells nad HRS from GEX data per sample
      sample_HRS_Bcells=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS", "Bcells"),]
      
      ### after running TRUST4 clustering, we pull from HPC the outcome; this file contains the corrected V/D/J genes onlu for HRS cells. So we need to merge it with the full file output from trust4 clustering
     # sample_trust4_per_chain= read.csv(paste0(sample,'_FILTERED_out_clone_dominantChain_',DominantChain,'.csv'))
      sample_trust4_per_chain= read.csv(paste0(sample,'_FILTERED_out_clone_dominantChain_',DominantChain,'.csv'))
      
      sample_trust4_per_chain= sample_trust4_per_chain %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
      colnames(sample_trust4_per_chain)[1:16]=c("X","consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )
      
      sample_trust4_per_chain = sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% sample_HRS_Bcells$cell_id1, ]
      
      #### we have VJ assignment and also getting the MainCelltype status from GEX data 
      if (DominantChain=="IGL" | DominantChain=="IGK"){
        
       # singlet_low_confidential= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/','*_low_confident_singlet_IGL_IGK_read_thre1.txt'))[1])
       # singlet= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/','*_singlet_IGL_IGK_read_thre1.txt'))[2])
       # doublet= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/','*_doublet_IGL_IGK_read_thre1.txt')))
       # pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/REPEAT_other_doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGK_IGL.csv')))
        pilot_cellStatus_file= read.csv(light_chain_cellStatus) 
        singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
        
        sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$corrected_Vgene," ",sample_trust4_per_chain$corrected_Jgene)
        ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
        sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind]
        sample_trust4_per_chain=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% singlet,]
        
      }else{
        
        # singlet_low_confidential= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/','*_low_confident_singlet_IGH_read_thre1.txt'))[1])
        # singlet= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/','*_singlet_IGH_read_thre1.txt'))[2])
        # doublet= read.table(Sys.glob(paste0('./doublets/combined_light_chains_automated/','*_doublet_IGH_read_thre1.txt')))
        # 
       # pilot_cellStatus_file= read.csv(Sys.glob(paste0('/Users/saramoein/Documents/new_run_HL_May2025/REPEAT_other_doublets_BCR_thre07_ent08/',sample,'_Rawdata_IGH.csv')))
        pilot_cellStatus_file= read.csv(heavy_chain_cellStatus) 
        singlet= unique(pilot_cellStatus_file$cell_id1[pilot_cellStatus_file$final_ReadStatus=="singlet"])
        
        sample_trust4_per_chain$VDJ= paste0(sample_trust4_per_chain$corrected_Vgene," ",sample_trust4_per_chain$corrected_Dgene," ",sample_trust4_per_chain$corrected_Jgene)
        ind=match(sample_trust4_per_chain$cell_id1,sample_HRS_Bcells$cell_id1)
        sample_trust4_per_chain$MainCelltype= sample_HRS_Bcells$MainCelltype[ind] 
        sample_trust4_per_chain=sample_trust4_per_chain[sample_trust4_per_chain$cell_id1 %in% singlet,]
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
        
        #column_ha = HeatmapAnnotation(cell_status = new_data$MainCelltype, V_gene= new_data$V_gene , col= list(cell_status=c("HRS"="red","Bcells"= "green")))
        #row_ha = rowAnnotation(cell_status= new_data$MainCelltype, j_gene= new_data$J_gene, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
        column_ha = HeatmapAnnotation(cell_status = new_data$MainCelltype, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
        row_ha = rowAnnotation(cell_status= new_data$MainCelltype, col= list(cell_status=c("HRS"="red","Bcells"= "green")))
        
      
        pdf(paste0(sample,'BCR_CLONE_heatmap_lv_distance_CDR3_',DominantChain,'_afterVJcorrection_dend_for_paper.pdf'), width=13, height=13)
        print(Heatmap(dist1,col= col_runif,top_annotation = column_ha, right_annotation = row_ha,show_column_names = FALSE,show_row_names = FALSE,cluster_columns=TRUE,show_row_dend = TRUE, raster_resize_mat = TRUE))
        dev.off() 
        
        
#}

