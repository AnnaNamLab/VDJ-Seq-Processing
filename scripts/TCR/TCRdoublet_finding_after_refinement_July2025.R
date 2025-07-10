library(data.table)
require(graphics)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(Polychrome)
library(tidyr)
library(stringdist)
library(ggplot2)


get_overlap_by_type <- function(df1, df2, dataset_name) {
  df1_HRS <- df1 %>% filter(cellType == "HRS") %>% select(cell_id1) %>% distinct()
  df2_HRS <- df2 %>% filter(cellType == "HRS") %>% select(cell_id1) %>% distinct()
  
  df1_Tcells <- df1 %>% filter(cellType == "Tcells") %>% select(cell_id1) %>% distinct()
  df2_Tcells <- df2 %>% filter(cellType == "Tcells") %>% select(cell_id1) %>% distinct()
  
  # Type A
  HRS_common <- intersect(df1_HRS$cell_id1, df2_HRS$cell_id1)
  HRS_all <- union(df1_HRS$cell_id1, df2_HRS$cell_id1)
  HRS_uncommon <- setdiff(HRS_all, HRS_common)
  
  # Type B
  Tcells_common <- intersect(df1_Tcells$cell_id1, df2_Tcells$cell_id1)
  Tcells_all <- union(df1_Tcells$cell_id1, df2_Tcells$cell_id1)
  Tcells_uncommon <- setdiff(Tcells_all, Tcells_common)
  
  summary_df <- data.frame(
    cellType = rep(c("HRS", "Tcells"), each = 2),
    overlap = rep(c("Common", "Uncommon"), 2),
    count = c(length(HRS_common), length(HRS_uncommon),
              length(Tcells_common), length(Tcells_uncommon))
  )
  
  
  return(list(
    summary = summary_df,
    HRS_all = HRS_all,
    Tcells_common = Tcells_common
  ))
  
}

kurtosis <- function (x) {
  n <- length(x)
  K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
  K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
  return(K)
}

normalized_entropy <- function(x) {
  probs <- x / sum(x)
  -sum(probs * log(probs)) / log(length(probs))
}

args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
library(optparse)
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to TRUST4 output directory"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to cell metadata CSV file"),
  make_option(c("-s", "--sample"), type = "character", help = "Sample name"),
  make_option(c("-f", "--file"), type = "character", help = "cdr3.out file name (e.g., sample_cdr3.out)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$input) || is.null(opt$metadata) || is.null(opt$sample) || is.null(opt$file)) {
  stop("Please provide --input, --metadata, --sample, and --file arguments.")
}

# Assign input paths
input_dir <- opt$input
metadata_file <- opt$metadata
sample <- opt$sample
cdr3_file <- opt$file


thre=0.7
entropy_ther=0.8
sub_folder= "TEST_doublets_TCR_thre07_ent08"
dir.create(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder))
dir.create(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/percentage'))


#for (sample in c("HL1","HL10","HL21","HL8","HL8R","HL15","HL20","HL16","HL12")){
  
  
  #setwd(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sample,"_TCR"))
  setwd(input_dir)
  #dir.create('./doublets_TCR_08ther/combined_lTRBt_chains_automated')
  
  ### using the trust4 raw_cdr3.out file as the input
  #out_FR2_cdr3=read.table(Sys.glob("*cdr3.out"),sep='\t')
  out_FR2_cdr3 <- read.table(cdr3_file, sep = '\t')
  
  
  colnames(out_FR2_cdr3) =c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "complete_vdj_assembly")
  data_raw= out_FR2_cdr3
  # extracting the cell_id and chain// the analysis on BCR is per chain: IGL/TRB/IGK
  data_raw$cell_id1=substr(data_raw$consensus_id,1,16)
  data_raw$chain=substr(data_raw$V_gene,1,3)
  

  
  #### extracting the Bcell, Tcell and HRS from GEX data
  
  #new_clusters= read.csv('/Users/saramoein/Documents/new_run_HL_May2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
  new_clusters <- read.csv(metadata_file)
  
  new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)
  
  ALLcells = new_clusters[new_clusters$Patient ==sample ,]
  
  
  out_FR2_cdr3$cell_id1=substr(out_FR2_cdr3$consensus_id,1,16)
  complete_out_FR2_cdr3 = out_FR2_cdr3[which(out_FR2_cdr3$CDR3_score> 0 ),] #### selecting the contigs with complete CDR3
  
  
  ############################BCR_Bcells 
  
  
  chain_names=c("TRA", "TRB")
  
  for (j in 1:2){ ### THROUGH 3 CHAINS IGL_IGK/TRB
    
    if (j==2){
      chain= c("TRB")
      chain_name="TRB"
    }else {
      chain=c("TRA")
      chain_name="TRA"
    }
    
    ### FILTERING THE BCR DATA BASED ON BCELL, TCELL, HRS
    data_rawCDR3= data_raw [data_raw$chain %in% chain,]
    rawCDR3_ALLcells= data_rawCDR3[data_rawCDR3$cell_id1 %in% ALLcells$cell_id1,]
    
    
    ### VDJ GENES ARE COMBINED TO SHOW SIMILARITY OF VDJ IN CONTIGS// THIS IS ONE OF THE CONDITIONS... 
    ### the VDJ are based main group of vdj and the sub_part fater the star is not considered to define the VDJ name.
    rawCDR3_ALLcells$main_vcall_family= sub("(.*?)[\\.|*].*", "\\1", rawCDR3_ALLcells$V_gene)
    rawCDR3_ALLcells$main_dcall_family= sub("(.*?)[\\.|*].*", "\\1", rawCDR3_ALLcells$D_gene)
    rawCDR3_ALLcells$main_jcall_family= sub("(.*?)[\\.|*].*", "\\1", rawCDR3_ALLcells$J_gene)
    rawCDR3_ALLcells$main_vdjcall= paste0(rawCDR3_ALLcells$main_vcall_family," ", rawCDR3_ALLcells$main_dcall_family," ", rawCDR3_ALLcells$main_jcall_family)
    
    ## removing the contigs with incomplete CDR3 (removing the missing values)
    if (j==2){
      rawCDR3_ALLcells= rawCDR3_ALLcells[(rawCDR3_ALLcells$V_gene!="*" & rawCDR3_ALLcells$D_gene!="*" & rawCDR3_ALLcells$J_gene!="*"),]
    } else {
      rawCDR3_ALLcells= rawCDR3_ALLcells[(rawCDR3_ALLcells$V_gene!="*" & rawCDR3_ALLcells$J_gene!="*"),]
    }
    
    rawCDR3_ALLcells= rawCDR3_ALLcells[which(rawCDR3_ALLcells$CDR3_score > 0), ]
    
    ### OBTAINING THE TOTAL CELLS ON EACH CHAIN
    cells=data.frame(unique(rawCDR3_ALLcells$cell_id1))
    colnames(cells)="id"
    
    #### AT THE END OF RUNNING, CELL ARE EDIVIDED TO 3 GROUPS: NOISE, DOUBLETS, SIGLETS
    singlet=matrix(NA,1,3)
    doublet=matrix(NA,1,3)
    low_confident_singlet=matrix(NA,1,3)
    NA_list=matrix(NA,1,3)
    singlet_dominant=matrix(NA,1,3)
    
    #### FIRST CONDITION IS TO CHECK IF threE IS A DOMINANT CONTIG BASED ON THE THRESHOLD
    #### EACH CONTIG HAS A READ COUNT IN CDR3.OUT FILE, COLUMN read_fragment_count
    #### IF A CONTIG IS BASED ON READS MORE THAN THE DEFINED THRESHOLD (thre) , then it is a singlet, since it has the 
    ### dominant number of reads
    
    
    #### one step After any step is to group_by contigs based on similar VDJ
    
    read_thre=1
    
    inf_ka = c()
    NA_list= c()
    
    for (i in 1:length(cells$id)){ ## per cell we decide what is it. Each cell can be either singlet, or doublet or noise
      
      sub_raw_dat=rawCDR3_ALLcells[which(rawCDR3_ALLcells$cell_id1==cells$id[i]),] ### first all contigs for that cell are selected
      agg_tbl <- sub_raw_dat %>% group_by(main_vdjcall) %>% summarise(read_fragment_count_group = sum(read_fragment_count)) ## the VDJs of contigs per that cells are aggregated to know how many reads we have per each VDJ
      sum1= sum(agg_tbl$read_fragment_count_group)
      agg_tbl$readRate=agg_tbl$read_fragment_count_group/sum1 ## readRate shows what percentage of the total reads belong to each contig.
      dominant_singlet=agg_tbl[agg_tbl$readRate >= thre,]
      ka= kurtosis(agg_tbl$read_fragment_count_group)
      ent= normalized_entropy(agg_tbl$read_fragment_count_group)
      
      
      if ((max(agg_tbl$read_fragment_count_group) == read_thre) &(nrow(agg_tbl)==1)){#### if the cell's maximum read count is less than the read_thre (here is 2) then that cell is noise
        
        low_confident_singlet=rbind(low_confident_singlet,c(cells$id[i],ka,ent))
        
      } else if (length(agg_tbl$read_fragment_count_group)==1) {
        singlet=rbind(singlet,c(cells$id[i],ka,ent))
        
      } else if (nrow(dominant_singlet)>0){ ## if three is dominant contig, then the cell is singlet
        #print('singlet')
        singlet=rbind(singlet,c(cells$id[i],ka,ent))
        singlet_dominant=rbind(singlet_dominant,c(cells$id[i],ka,ent))
        
      } else if (is.na(ka) | (is.infinite(ka) & sum(agg_tbl$read_fragment_count_group)<=3)){ 
        #  
        NA_list=rbind(NA_list,c(cells$id[i],ka,ent))
        if (is.infinite(ka)){
          inf_ka=rbind(cells$id[i],inf_ka)
         # print(ka)
        }
        
      } else if (ka <= 0){
        
        doublet=rbind(doublet,c(cells$id[i],ka,ent))
        
      } else if (ka > 0){
        
        singlet=rbind(singlet,c(cells$id[i],ka,ent))
        
      }  
      
      
      
    }
    singlet = na.omit( singlet)
    doublet= na.omit( doublet)
    low_confident_singlet = na.omit(low_confident_singlet)
    singlet_dominant= na.omit(singlet_dominant)
    ### The singlets should be only among Bcell and HRS. Thoise singlets that are not among HRS and Bcells are added to doublets.
    Bcell = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("Bcells"),]
    HRS = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("HRS"),]
    Tcell = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("Tcells"),]
    Macro_Mono=new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("Macro_Mono"),]
    DC=new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("DC"),]
    
    #nonHRSnonBcell = new_clusters[new_clusters$Patient == sample & !(new_clusters$MainCelltype %in% c("HRS", "Bcells")),]
    
    
    all.singlet <- rbind(singlet,low_confident_singlet)
    
    other_doublet= all.singlet[which(!all.singlet[,1] %in% Tcell$cell_id1== TRUE),]
    
    Tcell_singlet= all.singlet[which(all.singlet[,1]  %in% Tcell$cell_id1==TRUE),]
    doublet= rbind(doublet,other_doublet)
    
    other_doublet = data.frame(other_doublet)
    #colnames(other_doublet) = c("cell_id1","ka","Entropy")
    colnames(singlet) = c("cell_id1","ka","Entropy")
    colnames(doublet) = c("cell_id1","ka","Entropy")
    colnames(low_confident_singlet) = c("cell_id1","ka","Entropy")
    colnames(NA_list) = c("cell_id1","ka","Entropy") 
    colnames(Tcell_singlet) = c("cell_id1","ka","Entropy") 
    colnames(singlet_dominant) = c("cell_id1","ka","Entropy")
    colnames(other_doublet) = c("cell_id1","ka","Entropy")
    ###############################
    ############### we extract the maximum read contig per each cell and filter the original file
    
    library(dplyr)
    
    # Step 1: Aggregate read_fragment_count by main_vdjcall within each cell_id1
    vdj_aggregated <- rawCDR3_ALLcells %>%
      group_by(cell_id1, main_vdjcall) %>%
      summarise(aggregated_read_count = sum(read_fragment_count, na.rm = TRUE)) %>%
      arrange(cell_id1, desc(aggregated_read_count))  # Sort by read count within each cell_id1
    
    # Step 2: Calculate the total number of reads and unique main_vdjcall per cell_id1
    cell_info <- rawCDR3_ALLcells %>%
      group_by(cell_id1) %>%
      summarise(
        total_reads = sum(read_fragment_count, na.rm = TRUE),  # Total reads per cell_id1
        unique_vdjcall_count = n_distinct(main_vdjcall)  # Count of unique main_vdjcall per cell_id1
      )
    
    # Step 3: Get the top 2 main_vdjcall for each cell_id1 based on aggregated read counts
    top_vdj_calls <- vdj_aggregated %>%
      group_by(cell_id1) %>%
      slice_head(n = 10) %>%  # Get top 2 by read count
      ungroup()
    
    # Step 4: Pivot the data so that top 2 read counts are in separate columns
    
    vdj_counts_with_reads <- top_vdj_calls %>%
      group_by(cell_id1) %>%
      mutate(rank = rank(-aggregated_read_count, ties.method = "first")) %>%
      filter(rank <= 10) %>%
      ungroup() %>%
      pivot_wider(
        names_from = rank,
        values_from = aggregated_read_count,
        names_prefix = "read_",
        values_fill = 0
      )
    
    # Step 5: Merge the top 2 read counts with the cell_info data
    final_result <- cell_info %>%
      left_join(vdj_counts_with_reads, by = "cell_id1")  # Merge with the top VDJ calls data
    
    read_cols <- paste0("read_", 1:10)
    available_read_cols <- intersect(read_cols, colnames(final_result))
    
    
    
    final_aggregated_result <- final_result %>%
      group_by(cell_id1) %>%
      summarise(
        total_reads = first(total_reads),
        unique_vdjcall_count = first(unique_vdjcall_count),
        across(all_of(available_read_cols), ~sum(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    
    #rawCDR3_ALLcells = rawCDR3_ALLcells
    ###################### we add two columns to the cleaned file "cellType" for Tcell/Bcell/HRS and readStatus for signlet/doublet/noise  
    
    vdj_counts = final_aggregated_result
    vdj_counts$ReadStatus=NA
    ind= match(low_confident_singlet, vdj_counts$cell_id1)
    vdj_counts$ReadStatus[ind]= "low_confident_singlet"
    
    ind= match(doublet, vdj_counts$cell_id1)
    vdj_counts$ReadStatus[ind]= "doublet"
    
    ind= match(Tcell_singlet, vdj_counts$cell_id1)
    vdj_counts$ReadStatus[ind]= "singlet"
    
    ind= match(NA_list, vdj_counts$cell_id1)
    vdj_counts$ReadStatus[ind]= "Not_Assigned"
    
    
    vdj_counts$cellType=NA
    
    ind= match(Tcell$cell_id1, vdj_counts$cell_id1)
    vdj_counts$cellType[ind]= "Tcells"
    
    ind= match(Bcell$cell_id1, vdj_counts$cell_id1)
    vdj_counts$cellType[ind]= "Bcells"
    
    ind= match(HRS$cell_id1, vdj_counts$cell_id1)
    vdj_counts$cellType[ind]= "HRS"
    
    ind= match(Macro_Mono$cell_id1, vdj_counts$cell_id1)
    vdj_counts$cellType[ind]= "Macro_Mono"
    
    ind= match(DC$cell_id1, vdj_counts$cell_id1)
    vdj_counts$cellType[ind]= "DC"
    
    
    vdj_counts= data.frame(vdj_counts)
    vdj_counts$kurtosis = NA
    vdj_counts$ent =NA
    
    vdj_counts=vdj_counts[is.na(vdj_counts$cell_id1)==FALSE,]
    
    Tcell_singlet= data.frame(Tcell_singlet)
    doublet = data.frame(doublet)
    NA_list = data.frame(NA_list)
    low_confident_singlet = data.frame(low_confident_singlet)
    
    
    
    ind_singlet= match(Tcell_singlet$cell_id1 , vdj_counts$cell_id1)
    valid_indices <- !is.na(ind_singlet)
    vdj_counts$kurtosis[ind_singlet[valid_indices]] <- Tcell_singlet$ka[valid_indices]
    vdj_counts$ent[ind_singlet[valid_indices]] = Tcell_singlet$Entropy[valid_indices]
    
    ind_doublet= match(doublet$cell_id1, vdj_counts$cell_id1)
    valid_indices <- !is.na(ind_doublet)
    vdj_counts$kurtosis[ind_doublet[valid_indices]] <- doublet$ka[valid_indices]
    vdj_counts$ent[ind_doublet[valid_indices]] = doublet$Entropy[valid_indices]
    
    ind_NA= match(NA_list$cell_id1 , vdj_counts$cell_id1)
    valid_indices <- !is.na(ind_NA)
    vdj_counts$kurtosis[ind_NA[valid_indices]] <- NA_list$ka[valid_indices]
    vdj_counts$ent[ind_NA[valid_indices]] = NA_list$Entropy[valid_indices]
    
    
    ind_low_confident_singlet = match(low_confident_singlet$cell_id1 , vdj_counts$cell_id1)
    valid_indices <- !is.na(ind_low_confident_singlet)
    vdj_counts$kurtosis[ind_low_confident_singlet[valid_indices]] <- low_confident_singlet$ka[valid_indices]
    vdj_counts$ent[ind_low_confident_singlet[valid_indices]] = low_confident_singlet$Entropy[valid_indices]
    
    
    
    #vdj_counts= vdj_counts[vdj_counts$cellType!= "Bcells",]
    
    ind = match(vdj_counts$cell_id1 , vdj_aggregated$cell_id1)
    vdj_counts$chain = substr(vdj_aggregated$main_vdjcall[ind], 1,3)
    # 
    ind= match(vdj_counts$cell_id1 , new_clusters$cell_id1)
    vdj_counts$subType = new_clusters$SubtypeName[ind]
    
    
    vdj_counts=vdj_counts[!is.na(vdj_counts$cell_id1),]
    #### to detect the doubtful data
    
    vdj_counts$ReadStatus[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN")) & (vdj_counts$kurtosis=="Inf" | vdj_counts$kurtosis=="-Inf" | vdj_counts$kurtosis=="NaN" ) & (vdj_counts$total_reads >7) & (!vdj_counts$cell_id1 %in% singlet_dominant[,1]))] = "doublet"
    #vdj_counts$ReadStatus[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN")) & (vdj_counts$kurtosis=="Inf" | vdj_counts$kurtosis=="-Inf" | vdj_counts$kurtosis=="NaN" ) & (vdj_counts$total_reads >7) & (vdj_counts$cell_id1 %in% singlet[,1]))] = "singlet"
    vdj_counts$ReadStatus[which((vdj_counts$ent =="NaN")  & (vdj_counts$kurtosis=="NaN" ) & (vdj_counts$unique_vdjcall_count==1) )] = "singlet"
    
    
    weak_data_vdj_count_doublet = vdj_counts[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN"))  & (vdj_counts$kurtosis=="Inf" | vdj_counts$kurtosis=="-Inf" ) & vdj_counts$total_reads >7 ),]
    weak_data_vdj_count_NA = vdj_counts[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN"))  & (vdj_counts$kurtosis=="Inf" | vdj_counts$kurtosis=="-Inf" | vdj_counts$kurtosis=="NaN") & vdj_counts$total_reads <=7 ),]
    
    # sub_n1 <- nrow(weak_data_vdj_count_doublet)
    # sub_n2 <- nrow(weak_data_vdj_count_NA)
    # total_n <- nrow(vdj_counts)
    # other_n <- total_n - sub_n1-sub_n2
    
    # if (sub_n2/total_n < 0.03){
    #   vdj_counts$ReadStatus[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN")) & ((vdj_counts$kurtosis=="Inf") | (vdj_counts$kurtosis=="NaN")) & (vdj_counts$total_reads<=7) )] = "Not_Assigned"
    #   vdj_counts$ReadStatus[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN")) & (vdj_counts$kurtosis=="-Inf" ) & (vdj_counts$total_reads<=7) )] = "Not_Assigned"
    #   print("smaller")
    # } else {
      vdj_counts$ReadStatus[which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN")) & (vdj_counts$kurtosis=="Inf" | vdj_counts$kurtosis=="-Inf" ) & (vdj_counts$total_reads<=7) )] = "if_Not_Assigned"
      vdj_counts$ReadStatus[which(((vdj_counts$ent<=entropy_ther) & (vdj_counts$ent!= "NaN")) & (vdj_counts$kurtosis=="Inf"  | vdj_counts$kurtosis=="-Inf"  ) & (vdj_counts$total_reads<=7) )] = "singlet"
      
      ind = which(((vdj_counts$ent>entropy_ther) & (vdj_counts$ent!= "NaN")) & (vdj_counts$kurtosis=="Inf" | vdj_counts$kurtosis=="-Inf" ) & (vdj_counts$total_reads<=7))
      if_Not_Assigned= vdj_counts[ind,]
      
      sub_if_Not_Assigned = if_Not_Assigned[, c(1, (ncol(if_Not_Assigned)-4):(ncol(if_Not_Assigned)-2),ncol(if_Not_Assigned))]
      data_sub_if_Not_Assigned <- rawCDR3_ALLcells[rawCDR3_ALLcells$cell_id1 %in% sub_if_Not_Assigned$cell_id1,]
      data_sub_if_Not_Assigned_with_lv <- data_sub_if_Not_Assigned %>%
        group_by(cell_id1) %>%
        # Select top 2 contigs per cell by reads
        slice_max(order_by = read_fragment_count, n = 2, with_ties = FALSE) %>%
        # Compute normalized LV distance between CDR3s
        mutate(CDR3_LV_dist_col = if(n() == 2) {
          # Extract the two CDR3s
          cdr3s <- CDR3
          lv <- stringdist(cdr3s[1], cdr3s[2], method = "lv")
          norm_lv <- lv / max(nchar(cdr3s[1]), nchar(cdr3s[2]))
          rep(norm_lv, 2)
        } else {
          NA_real_
        }) %>%
        ungroup()
      
      singlet_selected= unique(data_sub_if_Not_Assigned_with_lv$cell_id1[which(data_sub_if_Not_Assigned_with_lv$CDR3_LV_dist_col<=0.2)])
      ind = match(singlet_selected , vdj_counts$cell_id1)
      vdj_counts$ReadStatus[ind]="singlet"
      
      vdj_counts$ReadStatus = ifelse(vdj_counts$ReadStatus=="if_Not_Assigned","Not_Assigned",vdj_counts$ReadStatus)
      
    #}
    
    
    
    # This step is to define the confidence level of the doublets
    
    doublet= vdj_counts$cell_id1[which(vdj_counts$ReadStatus=="doublet")]
    singlet= vdj_counts$cell_id1[which(vdj_counts$ReadStatus=="singlet")]
    NAs= vdj_counts$cell_id1[which(vdj_counts$ReadStatus=="Not_Assigned")]
    rawCDR3_ALLcells$ReadStatus="NA"
    other_cells_NA= vdj_counts$cell_id1[which((vdj_counts$ReadStatus=="Not_Assigned") | (vdj_counts$unique_vdjcall_count==1 & vdj_counts$total_reads<=3 & vdj_counts$cellType %in% c("Macro_Mono","Bcells","DC")))]
    rawCDR3_ALLcells$ReadStatus[match(doublet , rawCDR3_ALLcells$cell_id1)] = "doublet"
    rawCDR3_ALLcells$ReadStatus[match(singlet , rawCDR3_ALLcells$cell_id1)] = "singlet"
    rawCDR3_ALLcells$ReadStatus[match(NAs , rawCDR3_ALLcells$cell_id1)] = "Not_Assigned"
    rawCDR3_ALLcells$ReadStatus[match(other_cells_NA , rawCDR3_ALLcells$cell_id1)] = "Not_Assigned"
    rawCDR3_ALLcells$CDR1_length = NA
    rawCDR3_ALLcells$CDR2_length = NA
    rawCDR3_ALLcells$CDR1_length = nchar(rawCDR3_ALLcells$CDR1)
    rawCDR3_ALLcells$CDR2_length = nchar(rawCDR3_ALLcells$CDR2)
    
    
    # This step is to impute the status of the each contig based on the available contig's status
    rawCDR3_ALLcells <- rawCDR3_ALLcells %>%
      group_by(cell_id1) %>%
      mutate(status_impute = ifelse(
        ReadStatus== "NA",
        unique(ReadStatus[ReadStatus!="NA"])[1],
        ReadStatus
      )) %>%
      ungroup()
    
    rawCDR3_ALLcells$doublet_confidence = rawCDR3_ALLcells$status_impute
    
    
    # if TRB, then the CDR1 and CDR2 has different length comparing to IGL/IGK
    if (j==2){
      write.csv(vdj_counts,paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_TRB_",'countVDJ_perCELL_AfterInfFixing_POLISH_addingNA','.csv'),row.names = FALSE )#sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE
      
      ind_high= which((rawCDR3_ALLcells$status_impute =="doublet") & (rawCDR3_ALLcells$complete_vdj_assembly==1))
      ind_med= which((rawCDR3_ALLcells$status_impute =="doublet") &  (rawCDR3_ALLcells$CDR1_length <= 23  | rawCDR3_ALLcells$CDR2_length <= 20 | (rawCDR3_ALLcells$complete_vdj_assembly==0)) & (rawCDR3_ALLcells$CDR1_length !=1  | rawCDR3_ALLcells$CDR2_length !=1))
      ind_low = which(rawCDR3_ALLcells$status_impute =="doublet" & (rawCDR3_ALLcells$CDR1_length ==1  & rawCDR3_ALLcells$CDR2_length ==1))
      rawCDR3_ALLcells$doublet_confidence[ind_high]="doublet_high"
      rawCDR3_ALLcells$doublet_confidence[ind_med]="doublet_medium"
      rawCDR3_ALLcells$doublet_confidence[ind_low]="doublet_low"
      
    }else{
      write.csv(vdj_counts,paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_TRA_",'countVDJ_perCELL_AfterInfFixing_POLISH_addingNA','.csv'), row.names = FALSE)#sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE
      
      ind_high= which(rawCDR3_ALLcells$status_impute =="doublet"  & (rawCDR3_ALLcells$complete_vdj_assembly==1))
      ind_med= which(rawCDR3_ALLcells$status_impute =="doublet" &  (rawCDR3_ALLcells$CDR1_length <= 17  | rawCDR3_ALLcells$CDR2_length <= 8 | (rawCDR3_ALLcells$complete_vdj_assembly==0)) & (rawCDR3_ALLcells$CDR1_length !=1  | rawCDR3_ALLcells$CDR2_length !=1))
      ind_low = which(rawCDR3_ALLcells$status_impute =="doublet" & (rawCDR3_ALLcells$CDR1_length ==1  & rawCDR3_ALLcells$CDR2_length ==1))
      rawCDR3_ALLcells$doublet_confidence[ind_high]="doublet_high"
      rawCDR3_ALLcells$doublet_confidence[ind_med]="doublet_medium"
      rawCDR3_ALLcells$doublet_confidence[ind_low]="doublet_low"
      
    }
    
    df = rawCDR3_ALLcells
    
    
    sub_vdj_counts = vdj_counts[, c(1, (ncol(vdj_counts)-4):(ncol(vdj_counts)-2),ncol(vdj_counts))]
    df <- df %>%
      left_join(sub_vdj_counts, by = "cell_id1")
    
    # Step 1: Extract only doublets
    df_doublets <- df %>%
      filter(status_impute == "doublet")
    
    # Step 2: Get top 2 contigs per cell based on reads
    top2 <- df_doublets %>%
      group_by(cell_id1) %>%
      arrange(desc(read_fragment_count)) %>%
      slice_head(n = 2) %>%
      mutate(rank = row_number()) %>%
      ungroup()
    
    # Step 3: Pivot to wide to evaluate confidence
    top2_wide <- top2 %>%
      select(cell_id1, rank, doublet_confidence) %>%
      pivot_wider(names_from = rank, values_from = doublet_confidence, names_prefix = "conf_")
    
    
    
    df_with_lv <- df_doublets %>%
      group_by(cell_id1) %>%
      # Select top 2 contigs per cell by reads
      slice_max(order_by = read_fragment_count, n = 2, with_ties = FALSE) %>%
      # Compute normalized LV distance between CDR3s
      mutate(CDR3_LV_dist = if(n() == 2) {
        # Extract the two CDR3s
        cdr3s <- CDR3
        lv <- stringdist(cdr3s[1], cdr3s[2], method = "lv")
        norm_lv <- lv / max(nchar(cdr3s[1]), nchar(cdr3s[2]))
        rep(norm_lv, 2)
      } else {
        NA_real_
      }) %>%
      ungroup()
    
    
    library(dplyr)
    
    df_with_lv <- df_with_lv %>%
      group_by(cell_id1) %>%
      mutate(
        conf1 = first(doublet_confidence),
        conf2 = last(doublet_confidence),
        lv_dist = unique(CDR3_LV_dist),
        last_status = case_when(
          (conf1 == "doublet_high" & conf2 == "doublet_low" | conf1 == "doublet_low" & conf2 == "doublet_high") & lv_dist <= 0.2  ~ "singlet",
          (conf1 == "doublet_high" & conf2 == "doublet_low" | conf1 == "doublet_low" & conf2 == "doublet_high") & lv_dist > 0.2 ~ "doublet",
          (conf1 == "doublet_high" & conf2 == "doublet_medium" | conf1 == "doublet_medium" & conf2 == "doublet_high") & lv_dist > 0.2 ~ "doublet",
          (conf1 == "doublet_high" & conf2 == "doublet_medium" | conf1 == "doublet_medium" & conf2 == "doublet_high")  & lv_dist <= 0.2 ~ "singlet",
          (conf1 == "doublet_low" & conf2 == "doublet_medium" | conf1 == "doublet_medium" & conf2 == "doublet_low") & lv_dist <= 0.2 ~ "singlet",
          (conf1 == "doublet_low" & conf2 == "doublet_medium" | conf1 == "doublet_medium" & conf2 == "doublet_low")  & lv_dist > 0.2 ~ "doublet",
          conf1 == "doublet_high" & conf2 == "doublet_high" ~ "doublet",
          conf1 == "doublet_low" & conf2 == "doublet_low"  & lv_dist > 0.2 ~ "doublet",
          conf1 == "doublet_low" & conf2 == "doublet_low" &  lv_dist <= 0.2 ~ "singlet",
          conf1 == "doublet_medium" & conf2 == "doublet_medium" & lv_dist > 0.2 ~ "doublet",
          conf1 == "doublet_medium" & conf2 == "doublet_medium" & lv_dist <= 0.2 ~ "singlet",
          TRUE ~ NA_character_
        )
      ) %>%
      ungroup()
    
    
    
    
    
    # Step 5: Merge classification back into the full dataframe
    df <- df %>%
      left_join(df_with_lv %>% select(cell_id1 ,CDR3_LV_dist,last_status), by = "cell_id1")
    df <- df %>% distinct()
    df$final_ReadStatus = ifelse(is.na(df$last_status), df$status_impute, df$last_status)
    rawCDR3_ALLcells = df
    
    
    
    colnames(rawCDR3_ALLcells)[1:13] =c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "complete_vdj_assembly") 
    
    
    ind_other_doublet=which(rawCDR3_ALLcells$cellType %in% c("Macro_Mono","DC","Bcells") & (rawCDR3_ALLcells$final_ReadStatus=="singlet"))
    rawCDR3_ALLcells$final_ReadStatus[ind_other_doublet]="doublet"
    ind_HRS= which(rawCDR3_ALLcells$cellType %in% c("HRS"))
    rawCDR3_ALLcells$final_ReadStatus[ind_HRS]="singlet"
    
    write.csv(rawCDR3_ALLcells ,paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_Rawdata_",chain_name,".csv"),row.names = FALSE )#sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE
    
    ind_status = match(vdj_counts$cell_id1, rawCDR3_ALLcells$cell_id1)
    vdj_counts$final_ReadStatus= rawCDR3_ALLcells$final_ReadStatus[ind_status]
    
    
    
    if (j==1){
      
      light_vdj_counts= vdj_counts
    }else{
      heavy_vdj_counts= vdj_counts
      
    }
    
    weak_data_vdj_count_NA= vdj_counts[which(vdj_counts$final_ReadStatus=="Not_Assigned"),]
    weak_data_vdj_count_doublet= vdj_counts[which(vdj_counts$final_ReadStatus=="doublet"),]
    sub_n1 <- nrow(weak_data_vdj_count_doublet)
    sub_n2 <- nrow(weak_data_vdj_count_NA)
    total_n <- nrow(vdj_counts)
    other_n <- total_n - sub_n1-sub_n2
    
    df <- data.frame(
      category = c("doublet", "not_assigned","singlet"),
      count = c(sub_n1, sub_n2 ,other_n)
    )
    
    # Calculate percentages
    df$percent <- round(100 * df$count / sum(df$count), 1)
    df$label <- paste0(df$category, " (", df$percent, "%)")
    
    
    a0 <- ggplot(df, aes(x = "", y = count, fill = category)) +
      geom_col(width = 1, color = "white") +
      coord_polar(theta = "y") +
      geom_text(aes(label = label),
                position = position_stack(vjust = 1),
                size = 3) +
      scale_fill_manual(values = c("skyblue", "pink","lightgray")) +
      labs(title = "Proportion of cells with insufficient reads", x = NULL, y = NULL) +
      theme_void()
    

    
    long_reads <- weak_data_vdj_count_doublet %>%
      select(cell_id1, starts_with("read_")) %>%
      pivot_longer(cols = starts_with("read_"), 
                   names_to = "Read", 
                   values_to = "Count")
    
    long_reads= long_reads[long_reads$Count!=0 , ]  
    long_reads$Read <- factor(long_reads$Read, levels = paste0("read_", 1:10))
    
    c0 <- ggplot(long_reads, aes(x = Read, y = Count)) +
      geom_bar(stat = "identity", fill = "skyblue4") +
      facet_wrap(~ cell_id1, scales = "free_y", ncol = 10) +  # 10 columns grid
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      ) +
      labs(title = "Read Distributions per Doublet Cell ")
    
    
    
    long_reads <- weak_data_vdj_count_NA %>%
      select(cell_id1, starts_with("read_")) %>%
      pivot_longer(cols = starts_with("read_"), 
                   names_to = "Read", 
                   values_to = "Count")
    
    long_reads= long_reads[long_reads$Count!=0 , ]  
    long_reads$Read <- factor(long_reads$Read, levels = paste0("read_", 1:10))
    
    c1 <- ggplot(long_reads, aes(x = Read, y = Count)) +
      geom_bar(stat = "identity", fill = "pink4") +
      facet_wrap(~ cell_id1, scales = "free_y", ncol = 10) +  # 10 columns grid
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      ) +
      labs(title = "Read Distributions per Not_Assigned Cell ")
    
    
    set.seed(123)  # for reproducibility
    vdj_counts_singlet = vdj_counts[which(vdj_counts$final_ReadStatus=="singlet"),]
    sub_vdj_counts <- vdj_counts_singlet[sample(nrow(vdj_counts_singlet), size = 0.1 * nrow(vdj_counts_singlet), replace = FALSE), ]
    long_reads <- sub_vdj_counts %>%
      select(cell_id1, starts_with("read_")) %>%
      pivot_longer(cols = starts_with("read_"), 
                   names_to = "Read", 
                   values_to = "Count")
    
    long_reads= long_reads[long_reads$Count!=0 , ]  
    long_reads$Read <- factor(long_reads$Read, levels = paste0("read_", 1:10))
    
    c2 <- ggplot(long_reads, aes(x = Read, y = Count)) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      facet_wrap(~ cell_id1, scales = "free_y", ncol = 10) +  # 10 columns grid
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      ) +
      labs(title = "Read Distributions per Singlet Cell ")
    
    
    
    #### saving the noise/singlet/doublets to output files
    pdf(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_insufficient_reads_cellType_status_AfterInfFixing_POLISH_addingNA_",chain_name,".pdf"), width=20, height=10)  
    # Assume p, q, z, z1 may or may not exist or may be NULL
      print(a0)
    
    dev.off()
    pdf(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_insufficient_reads_cellType_status_AfterInfFixing_POLISH_addingNA_distributions_",chain_name,".pdf"), width=40, height=100)  
      print(c0)
      print(c1)
      print(c2)
    dev.off()
    
    
  }
  sub_heavy_vdj_counts= cbind(heavy_vdj_counts$cell_id1, heavy_vdj_counts$total_reads, heavy_vdj_counts$unique_vdjcall_count,
                              heavy_vdj_counts$final_ReadStatus, heavy_vdj_counts$cellType, heavy_vdj_counts$chain, heavy_vdj_counts$subType)
  
  sub_light_vdj_counts= cbind(light_vdj_counts$cell_id1,light_vdj_counts$total_reads, light_vdj_counts$unique_vdjcall_count,
                              light_vdj_counts$final_ReadStatus, light_vdj_counts$cellType, light_vdj_counts$chain, light_vdj_counts$subType)
  
  
  
  
  total_vdjdata = rbind( sub_heavy_vdj_counts ,  sub_light_vdj_counts )
  
  colnames(total_vdjdata) = c("cell_id1", "total_reads" , "unique_vdjcall_count" , "final_ReadStatus" , "cellType" , "chain" , "subType") 
  total_vdjdata = data.frame(total_vdjdata)
  total_vdjdata =  total_vdjdata[ total_vdjdata$final_ReadStatus=="doublet",]
  
  
  TRA= total_vdjdata[total_vdjdata$chain == "TRA",]
  TRB= total_vdjdata[total_vdjdata$chain == "TRB",]

  result <- get_overlap_by_type(TRA, TRB, "light_heavy")
  
  
  
  w1 <- ggplot(result$summary, aes(x = cellType, y = count, fill = overlap)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    geom_text(aes(label = count), 
              position = position_stack(vjust = 0.5), 
              color = "white", size = 4) +
    scale_fill_manual(values = c("Common" = "green", "Uncommon" = "brown")) +
    labs(title = "Overlap of cell_id1 in TRA and TRB by cellType",
         x = "Cell Type", y = "Number of Cells", fill = "Overlap Status") +
    theme_minimal(base_size = 14)
  
  
  write.csv(result$Tcells_common , paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_final_Tcell_doublet_common_light_heavy.csv"))
  write.csv(result$HRS_all , paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,'/',sample,"_final_HRS_doublet_either_light_or_heavy.csv"))
  
 

  pdf(paste0('/Users/saramoein/Documents/new_run_HL_May2025/',sub_folder,"/percentage/",sample,"_doublet_overlap_in_chains_per_celltype",".pdf"))
  print(w1)
  dev.off()
  
#}

