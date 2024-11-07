### code for finding doublets based on BCR data
library(data.table)
require(graphics)
require(pastecs)
library(dplyr)



sample="HL9"

setwd(paste0('/Users/saramoein/Documents/new_run_HL_August2024/',sample))
dir.create('./new_new_doublets_fixed')
dir.create('./new_new_doublets_fixed/combined_light_chains_automated')

### using the trust4 raw_cdr3.out file as the input
out_FR2_cdr3=read.table(Sys.glob("*cdr3.out"),sep='\t')
data_raw= out_FR2_cdr3

# extracting the cell_id and chain// the analysis on BCR is per chain: IGL/IGH/IGK
data_raw$cell_id1=substr(data_raw$V1,1,16)
data_raw$chain=substr(data_raw$V3,1,3)

kurtosis <- function (x) {
  n <- length(x)
  K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
  K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
  return(K)
}

#### extracting the Bcell, Tcell and HRS from GEX data

new_clusters= read.csv('/users/saramoein/downloads/2024-08-13_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)

ALLcells = new_clusters[new_clusters$Patient ==sample ,]


out_FR2_cdr3$cell_id1=substr(out_FR2_cdr3$V1,1,16)
complete_out_FR2_cdr3 = out_FR2_cdr3[which(out_FR2_cdr3$V10> 0 ),]


############################BCR_Bcells 

chain_names=c("IGL_IGK", "IGH")

for (j in 1:2){ ### THROUGH 3 CHAINS IGL_IGK/IGH
  
   if (j==1){
     chain= c("IGK","IGL")
   }else {
     chain="IGH"
   }
  
  ### FILTERING THE BCR DATA BASED ON BCELL, TCELL, HRS
  data_rawCDR3= data_raw [data_raw$chain %in% chain,]
  rawCDR3_ALLcells= data_rawCDR3[data_rawCDR3$cell_id1 %in% ALLcells$cell_id1,]
  
  
  ### VDJ GENES ARE COMBINED TO SHOW SIMILARITY OF VDJ IN CONTIGS// THIS IS ONE OF THE CONDITIONS... 
  ### the VDJ are based main group of vdj and the sub_part fater the star is not considered to define the VDJ name.
  rawCDR3_ALLcells$main_vcall_family= sub("(.*?)[\\.|*].*", "\\1", rawCDR3_ALLcells$V3)
  rawCDR3_ALLcells$main_dcall_family= sub("(.*?)[\\.|*].*", "\\1", rawCDR3_ALLcells$V4)
  rawCDR3_ALLcells$main_jcall_family= sub("(.*?)[\\.|*].*", "\\1", rawCDR3_ALLcells$V5)
  rawCDR3_ALLcells$main_vdjcall= paste0(rawCDR3_ALLcells$main_vcall_family," ", rawCDR3_ALLcells$main_dcall_family," ", rawCDR3_ALLcells$main_jcall_family)
  
  ## removing the contigs with incomplete CDR3 (removing the missing values)
  if (j==2){
    rawCDR3_ALLcells= rawCDR3_ALLcells[(rawCDR3_ALLcells$V3!="*" & rawCDR3_ALLcells$V4!="*" & rawCDR3_ALLcells$V5!="*"),]
  } else {
    rawCDR3_ALLcells= rawCDR3_ALLcells[(rawCDR3_ALLcells$V3!="*" & rawCDR3_ALLcells$V5!="*"),]
  }
  
  rawCDR3_ALLcells= rawCDR3_ALLcells[which(rawCDR3_ALLcells$V10 > 0), ]
  
  ### OBTAINING THE TOTAL CELLS ON EACH CHAIN
  cells=data.frame(unique(rawCDR3_ALLcells$cell_id1))
  colnames(cells)="id"
  
  #### AT THE END OF RUNNING, CELL ARE EDIVIDED TO 3 GROUPS: NOISE, DOUBLETS, SIGLETS
  singlet=vector()
  doublet=vector()
  low_confident_singlet=vector()
  
  
  #### FIRST CONDITION IS TO CHECK IF threE IS A DOMINANT CONTIG BASED ON THE THRESHOLD
  #### EACH CONTIG HAS A READ COUNT IN CDR3.OUT FILE, COLUMN V11
  #### IF A CONTIG IS BASED ON READS MORE THAN THE DEFINED THRESHOLD (thre) , then it is a singlet, since it has the 
  ### dominant number of reads
  
  
  #### one step before any step is to group_by contigs based on similar VDJ
  
  read_thre=1
  thre=0.8
  for (i in 1:length(cells$id)){ ## per cell we decide what is it. Each cell can be either singlet, or doublet or noise
    
    sub_raw_dat=rawCDR3_ALLcells[which(rawCDR3_ALLcells$cell_id1==cells$id[i]),] ### first all contigs for that cell are selected
    agg_tbl <- sub_raw_dat %>% group_by(main_vdjcall) %>% summarise(V11_group = sum(V11)) ## the VDJs of contigs per that cells are aggregated to know how many reads we have per each VDJ
    sum1= sum(agg_tbl$V11_group)
    agg_tbl$readRate=agg_tbl$V11_group/sum1 ## readRate shows what percentage of the total reads belong to each contig.
    dominant_singlet=agg_tbl[agg_tbl$readRate > 0.8,]
    ka= kurtosis(agg_tbl$V11_group)
    if (max(agg_tbl$V11_group) == read_thre){#### if the cell's maximum read count is less than the read_thre (here is 2) then that cell is noise
      
      low_confident_singlet=rbind(low_confident_singlet,cells$id[i])
    } else if (length(agg_tbl$V11_group)==1) {
      singlet=rbind(singlet,cells$id[i])
    
    } else if (nrow(dominant_singlet)>0){ ## if three is dominant contig, then the cell is singlet
      #print('singlet')
      singlet=rbind(singlet,cells$id[i])
    
    } else if (is.na(ka)){ 
      
      doublet=rbind(doublet,cells$id[i])
      
    } else if (ka <= 0){
      
      doublet=rbind(doublet,cells$id[i])
      
    } else if (ka > 0){
     
      singlet=rbind(singlet,cells$id[i])
      
    }  
    
    
    
  }
  ### The singlets should be only among Bcell and HRS. Thoise singlets that are not among HRS and Bcells are added to doublets.
  Bcell = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("Bcells"),]
  HRS = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("HRS"),]
  Tcell = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("Tcells"),]
  
  nonHRSnonBcell = new_clusters[new_clusters$Patient == sample & !(new_clusters$MainCelltype %in% c("HRS", "Bcells")),]
  
  
  all.singlet <- union(singlet,low_confident_singlet)
  
  T_doublet= all.singlet[!all.singlet %in% union(Bcell$cell_id1,HRS$cell_id1)]
  
  HRS_Bcell_singlet= all.singlet[all.singlet  %in% union(Bcell$cell_id1,HRS$cell_id1)]
  doublet= union(doublet,T_doublet)
  
  #### saving the noise/singlet/doublets to output files
   if (j==1){

     write.table(low_confident_singlet,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_low_confident_singlet_','IGL_IGK','_read_thre',read_thre,'.txt'),sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE)
     write.table(doublet,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_doublet_','IGL_IGK','_read_thre',read_thre,'.txt'),sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE)
     write.table(HRS_Bcell_singlet,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_singlet_','IGL_IGK','_read_thre',read_thre,'.txt'),sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE)
      
      } else {
      write.table(low_confident_singlet,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_low_confident_singlet_',chain,'_read_thre',read_thre,'.txt'),sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE)
      write.table(doublet,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_doublet_',chain,'_read_thre',read_thre,'.txt'),sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE)
      write.table(HRS_Bcell_singlet,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_singlet_',chain,'_read_thre',read_thre,'.txt'),sep='\t',row.names=FALSE, col.names = FALSE, quote=FALSE)
   }

  
  
  ###############################
  ############### we extract the maximum read contig per each cell and filter the original file
  
  
  rawCDR3_ALLcells_maximumReads = rawCDR3_ALLcells %>% group_by(cell_id1) %>% dplyr::slice(which.max(V11))  
 
  #rawCDR3_ALLcells_maximumReads = rawCDR3_ALLcells
  ###################### we add two columns to the cleaned file "cellType" for Tcell/Bcell/HRS and readStatus for signlet/doublet/noise  
  rawCDR3_ALLcells_maximumReads$ReadStatus=NA
  ind= match(low_confident_singlet,rawCDR3_ALLcells_maximumReads$cell_id1)
  rawCDR3_ALLcells_maximumReads$ReadStatus[ind]= "low_confident_singlet"
  
  ind= match(doublet,rawCDR3_ALLcells_maximumReads$cell_id1)
  rawCDR3_ALLcells_maximumReads$ReadStatus[ind]= "doublet"
  
  ind= match(HRS_Bcell_singlet,rawCDR3_ALLcells_maximumReads$cell_id1)
  rawCDR3_ALLcells_maximumReads$ReadStatus[ind]= "singlet"
  
  
  rawCDR3_ALLcells_maximumReads$cellType=NA
  Tcell = new_clusters[new_clusters$Patient ==sample & new_clusters$MainCelltype %in% c("Tcells"),]
  HRS = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("HRS"),]
  Bcell = new_clusters[new_clusters$Patient == sample & new_clusters$MainCelltype %in% c("Bcells"),]
  
  ind= match(Tcell$cell_id1,rawCDR3_ALLcells_maximumReads$cell_id1)
  rawCDR3_ALLcells_maximumReads$cellType[ind]= "Tcells"
  
  ind= match(Bcell$cell_id1,rawCDR3_ALLcells_maximumReads$cell_id1)
  rawCDR3_ALLcells_maximumReads$cellType[ind]= "Bcells"
  
  ind= match(HRS$cell_id1,rawCDR3_ALLcells_maximumReads$cell_id1)
  rawCDR3_ALLcells_maximumReads$cellType[ind]= "HRS"
  
  
  doublet_raw_file= out_FR2_cdr3[out_FR2_cdr3$cell_id1 %in% doublet,]
  
  #################### writing the filtered file after defining the readStatus
  if (j==2){
    write.csv(doublet_raw_file,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_BCR_data_cellType_doubletDetection_',chain,'_read_thre',read_thre,'_doubetRawFile.csv'))
    write.csv(rawCDR3_ALLcells_maximumReads,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_BCR_data_cellType_doubletDetection_',chain,'_read_thre',read_thre,'.csv'))
  } else {
      write.csv(rawCDR3_ALLcells_maximumReads,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_BCR_data_cellType_doubletDetection_','IGL_IGK','_read_thre',read_thre,'.csv'))
    write.csv(doublet_raw_file,paste0('./new_new_doublets_fixed/combined_light_chains_automated/',Sys.Date(),'_BCR_data_cellType_doubletDetection_','IGL_IGK','_read_thre',read_thre,'_doubetRawFile.csv'))
  }
  
  
  
}



