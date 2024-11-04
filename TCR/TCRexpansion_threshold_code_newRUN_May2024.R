################# CODE FOR DEFINING THE CLONE SIZE IN HRS TCELL CELLS.
################## IN THIS CODE IS BASED ON TWO DIFFERENT TCELL_HRS GROUPS: HRS_CD4T, HRS_CD8T 
################## A FISHER TEST BASED ON TCELL, HRST CELLS IN TWO DIFFERENT SIZES: SMALLER COLNE SIZE THAN A THRESHOD (FOR EXAMPLE THRESHOLD=3)
##### AND GREATER THAN THE THRESHOLD


require(dplyr)
library(ggplot2)

sample="HL20"
setwd(paste0('/Users/saramoein/Documents/new_run_HL_August2024/',sample,'_TCR'))
### reading the cell types from GEX

new_clusters=read.csv('/users/saramoein/downloads/2024-08-13_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv')
new_clusters$cell_id1=sub("^.*_([A-Z]+)-.*$", "\\1", new_clusters$Full.cell_id)


Tcell_HRS=new_clusters[new_clusters$Patient== sample & new_clusters$MainCelltype %in% c("HRS", "Tcells"),]

### after running TRUST4 clustering, we pull from HPC the outcome 
trust4_out_barcode_airr= read.table('out_clone_sara.tsv',header= FALSE, sep='\t')
colnames(trust4_out_barcode_airr)=c("consensus_id",	"index_within_consensus",	"V_gene",	"D_gene",	"J_gene",	"C_gene",	"CDR1",	"CDR2",	"CDR3",	"CDR3_score",	"read_fragment_count", "CDR3_germline_similarity", "full_length_assembly","contig_id","s" )


### from TRUST4 clustering we need to define the chain; V_gene is corresponded to V_gene
trust4_out_barcode_airr$chain=substr(trust4_out_barcode_airr$V_gene,1,3)
trust4_out_barcode_airr$cell_id1=substr(trust4_out_barcode_airr$contig_id,1,16)
trust4_out_barcode_airr= trust4_out_barcode_airr[trust4_out_barcode_airr$cell_id1 %in% Tcell_HRS$cell_id1,]



### seprating the locus per TRA, TRB



ind=match(trust4_out_barcode_airr$cell_id1,Tcell_HRS$cell_id1)
trust4_out_barcode_airr$malignantStatus=Tcell_HRS$MainCelltype[ind]



###### separating the data on TRA and TRB
TRA_new=trust4_out_barcode_airr[trust4_out_barcode_airr$chain=="TRA",] ### all the cells on TRA
TRA_new_Tcell=trust4_out_barcode_airr[trust4_out_barcode_airr$chain=="TRA" & trust4_out_barcode_airr$malignantStatus=="Tcells",] ###Tcells on TRA
chain="TRA"
### read contigs with 1 read (noise) and singlet
singlet_low_confidential= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/*_singlet_',chain,'_read_thre1.txt'))[1])
singlet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/*_singlet_',chain,'_read_thre1.txt'))[2])
doublet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/*_doublet_',chain,'_read_thre1.txt')))

TRA_new= TRA_new[TRA_new$cell_id1 %in% union(singlet_low_confidential$V1,singlet$V1),]
TRA_new_Tcell= TRA_new_Tcell[TRA_new_Tcell$cell_id1 %in% union(singlet_low_confidential$V1,singlet$V1),]
TRA_new= TRA_new %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
TRA_new_Tcell= TRA_new_Tcell %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 





TRB_new=trust4_out_barcode_airr[trust4_out_barcode_airr$chain=="TRB",]### all the cells on TRB
TRB_new_Tcell=trust4_out_barcode_airr[trust4_out_barcode_airr$chain=="TRB" & trust4_out_barcode_airr$malignantStatus=="Tcells",]###Tcells on TRB
TRB_new= TRB_new[TRB_new$D_gene!="",]
TRB_new_Tcell=TRB_new_Tcell[TRB_new_Tcell$D_gene!="",]
chain="TRB"
singlet_low_confidential= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/*_singlet_',chain,'_read_thre1.txt'))[1])
singlet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/*_singlet_',chain,'_read_thre1.txt'))[2])
doublet= read.table(Sys.glob(paste0('./new_new_doublets_fixed/combined_light_chains_automated/*_doublet_',chain,'_read_thre1.txt')))

TRB_new= TRB_new[TRB_new$cell_id1 %in% union(singlet_low_confidential$V1,singlet$V1),]
TRB_new_Tcell= TRB_new_Tcell[TRB_new_Tcell$cell_id1 %in% union(singlet_low_confidential$V1,singlet$V1),]
TRB_new= TRB_new %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count)) 
TRB_new_Tcell= TRB_new_Tcell %>% group_by(cell_id1) %>% dplyr::slice(which.max(read_fragment_count))



TRA_new$VJ=paste0(TRA_new$V_gene," ",TRA_new$J_gene)
TRA_new_Tcell$VJ=paste0(TRA_new_Tcell$V_gene," ",TRA_new_Tcell$J_gene)
TRB_new$VDJ=paste0(TRB_new$V_gene," ",TRB_new$D_gene," ",TRB_new$J_gene)
TRB_new_Tcell$VDJ=paste0(TRB_new_Tcell$V_gene," ",TRB_new_Tcell$D_gene," ",TRB_new_Tcell$J_gene)
trust4_out_barcode_airr$VJ=paste0(trust4_out_barcode_airr$V_gene," ",trust4_out_barcode_airr$J_gene)
trust4_out_barcode_airr$VDJ=paste0(trust4_out_barcode_airr$V_gene," ",trust4_out_barcode_airr$D_gene," ",trust4_out_barcode_airr$J_gene)

########
######### generating the clone size based on thresholds {3,4,5,6,7} and running fisher test per each threshold per chain per HRST cell type


total_CD4TwithTcell_TRA=matrix(NA,1,5)
total_Tcell_TRA=matrix(NA,1,5)
total_CD8TwithTcell_TRA=matrix(NA,1,5)
total_Tcell_TRA=matrix(NA,1,5)


total_CD8TwithTcell_TRB=matrix(NA,1,5)
total_Tcell_TRB=matrix(NA,1,5)
total_CD4TwithTcell_TRB=matrix(NA,1,5)
total_Tcell_TRB=matrix(NA,1,5)

###### clone size threshold {3:7}; that means what clones have the number of cells smaller than the threshold or greater than the threshold
##### we have two groups of TCR data: 1- those CD4T (or CD8T) segments with the same VDJ 2- Tcells with VDJs different from the CD4T VDJs
### so we have clonal (HRS_T cells) and non_clonal (Tcells). The fisher test table rows is based on the TOTAL cells in the clonal group and TOTAL cells in the non-clonal group
###### In the columns we divide based on the size of the clone. TOTAL cells in smaller than the thresholds, and TOTAL cells greater than the thresholds

thresh=c(3:7)

for (HRS_celltype in c("HRS-CD4T","HRS-CD8T")){
  
  for (i in c(1:5)){
    print(i)  
    ther=thresh[i]
    sample_HRS=Tcell_HRS[Tcell_HRS$SubtypeName == HRS_celltype,]
    TRB_new_HRS=TRB_new[TRB_new$cell_id1 %in% sample_HRS$cell_id1,]
    data= TRB_new_HRS
    
    ##### obtaining the VDJ of the cells belonging to HRS_cell type (CD4T or CD8T)
    ##### This step generate a data that cotains the clonal VDJS that has both HRS and Tcells
    data_withHRSandTcells=TRB_new[TRB_new$VDJ %in% data$VDJ,]
    
    #### generating the table based on cell count in clones per threshold
    ### grouping the clones to two groups: groups with clones that each one has cell number greater than threshold or 
    #### groups with clones that each one has cell number smaller than threshold 
    group_ALL_data_withHRSandTcells= data_withHRSandTcells %>% group_by(VDJ) %>% dplyr::count()
    group_ALL_data_withHRSandTcells_less_ther = group_ALL_data_withHRSandTcells[group_ALL_data_withHRSandTcells$n<ther,]
    group_ALL_data_withHRSandTcells_greater_ther = group_ALL_data_withHRSandTcells[group_ALL_data_withHRSandTcells$n>=ther,]
    
    
    ##### This step is to obtain any VDJ that does not have any HRS cell type in that (the Tcell group)
    #### After obtaining their VDJ, now we should extract any cells in TCR data (on TRB chain)  belonging to these VDJs
    data_Tcell_VDJ=setdiff(TRB_new_Tcell$VDJ,data_withHRSandTcells$VDJ)
    data_Tcell= TRB_new_Tcell[TRB_new_Tcell$VDJ %in% data_Tcell_VDJ,]
    
    group_data_Tcell= data_Tcell %>% group_by(VDJ) %>% dplyr::count()
    group_data_Tcell_less_ther=group_data_Tcell[group_data_Tcell$n<ther,]
    group_data_Tcell_greater_ther=group_data_Tcell[group_data_Tcell$n>=ther,]
    
    
    #########this step is to generate the matrix for running fisher test. It is a two by two matrix: 
    
    ##########################################################################################
    # number ofr cells from intersect the total  #  number of cells from intersect the total #
    # cells in clones with more                  #  cells in clones with less                #
    # than ther size and in HRSwith Tcells       #  than ther size and in HRSwith Tcell      #  
    ##########################################################################################
    #  number of cells from intersect the total  #   number of cell from intersect the total #
    #  cells in clones with more                 #   cells in clones with less               #
    #  than ther size and in Tcells              #   than ther size and in Tcells            #            
    ##########################################################################################
    
    total_cellCount_HRSwithTcell=sum(group_ALL_data_withHRSandTcells_less_ther$n)+sum(group_ALL_data_withHRSandTcells_greater_ther$n)
    total_Tcell= sum(group_data_Tcell_greater_ther$n) + sum(group_data_Tcell_less_ther$n)                                         
    mat_proportion=matrix(c(sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell , sum(group_ALL_data_withHRSandTcells_less_ther$n)/total_cellCount_HRSwithTcell,sum(group_data_Tcell_greater_ther$n)/total_Tcell , sum(group_data_Tcell_less_ther$n)/total_Tcell) , 2,2)
    mat_proportion=as.data.frame(mat_proportion)
    
    colnames(mat_proportion)<- c(paste0("Tcell_",HRS_celltype ," clone"), "Tcell_noHRS clone")
    rownames(mat_proportion)<- c(paste0("more than ", ther," cells"), paste0("less than ", ther," cells"))
    
    mat_count=matrix(c(sum(group_ALL_data_withHRSandTcells_greater_ther$n) , sum(group_ALL_data_withHRSandTcells_less_ther$n),sum(group_data_Tcell_greater_ther$n) , sum(group_data_Tcell_less_ther$n)) , 2,2)
    colnames(mat_count)<- c(paste0("Tcell_",HRS_celltype ,"clone"), "Tcell_noHRS clone")
    rownames(mat_count)<- c(paste0("more than ", ther," cells"), paste0("less than ", ther," cells"))
    fish_1= fisher.test(mat_count,alternative = "greater")
    write.table(mat_count,paste0('table_fisher_test_TRB_',HRS_celltype,'_ther',ther,'.csv'),row.names = TRUE, col.names = TRUE, quote=FALSE, sep=',')
    
    pdf(paste0('barplot_fisher_test_TRB_',HRS_celltype,'_ther',ther,'.pdf'))
    barplot(as.matrix(mat_proportion), legend= TRUE, main=paste0("p-value= ",fish_1$p.value))
    dev.off()
    
    if (HRS_celltype=="HRS-CD4T"){
      total_CD4TwithTcell_TRB[i]=sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell
      total_Tcell_TRB[i]=sum(group_data_Tcell_greater_ther$n)/total_Tcell ###noHRST : Tcell
      print(sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell)
      print(sum(group_data_Tcell_greater_ther$n)/total_Tcell)
    }
    else{
      total_CD8TwithTcell_TRB[i]=sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell
      total_Tcell_TRB[i]=sum(group_data_Tcell_greater_ther$n)/total_Tcell ###noHRST : Tcell
      print(sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell)
      print(sum(group_data_Tcell_greater_ther$n)/total_Tcell)
    }
    
    
    ########
    #########TRA###### SAME steps as above need to be repeated fro TRA. The only difference is that TRA does not have d_call gene.
    TRA_new_HRS=TRA_new[TRA_new$cell_id1 %in% sample_HRS$cell_id1,]
    data= TRA_new_HRS
    
    data_withHRSandTcells=TRA_new[TRA_new$VJ %in% data$VJ,]
    group_ALL_data_withHRSandTcells= data_withHRSandTcells %>% group_by(VJ) %>% dplyr::count()
    group_ALL_data_withHRSandTcells_less_ther=group_ALL_data_withHRSandTcells[group_ALL_data_withHRSandTcells$n<ther,]
    group_ALL_data_withHRSandTcells_greater_ther=group_ALL_data_withHRSandTcells[group_ALL_data_withHRSandTcells$n>=ther,]
    
    data_Tcell_VJ=setdiff(TRA_new_Tcell$VJ,data_withHRSandTcells$VJ)
    data_Tcell= TRA_new_Tcell[TRA_new_Tcell$VJ %in% data_Tcell_VJ,]
    
    group_data_Tcell= data_Tcell %>% group_by(VJ) %>% dplyr::count()
    group_data_Tcell_less_ther=group_data_Tcell[group_data_Tcell$n<ther,]
    group_data_Tcell_greater_ther=group_data_Tcell[group_data_Tcell$n>=ther,]
    
    ### normlized matrix
    total_cellCount_HRSwithTcell=sum(group_ALL_data_withHRSandTcells_less_ther$n)+sum(group_ALL_data_withHRSandTcells_greater_ther$n)
    total_Tcell= sum(group_data_Tcell_greater_ther$n) + sum(group_data_Tcell_less_ther$n)                                         
    mat_proportion=matrix(c(sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell , sum(group_ALL_data_withHRSandTcells_less_ther$n)/total_cellCount_HRSwithTcell,sum(group_data_Tcell_greater_ther$n)/total_Tcell , sum(group_data_Tcell_less_ther$n)/total_Tcell) , 2,2)
    mat_proportion=as.data.frame(mat_proportion)
    
    colnames(mat_proportion)<- c(paste0("Tcell_",HRS_celltype ,"clone"), "Tcell_noHRS clone")
    rownames(mat_proportion)<- c(paste0("more than ", ther," cells"), paste0("less than ", ther," cells"))
    
    #### original cell count matrix used by fisher test
    mat_count=matrix(c(sum(group_ALL_data_withHRSandTcells_greater_ther$n) , sum(group_ALL_data_withHRSandTcells_less_ther$n),sum(group_data_Tcell_greater_ther$n) , sum(group_data_Tcell_less_ther$n)) , 2,2)
    colnames(mat_count)<- c(paste0("Tcell_",HRS_celltype ," clone"), "Tcell_noHRS clone")
    rownames(mat_count)<- c(paste0("more than ", ther," cells"), paste0("less than ", ther," cells"))
    fish_1= fisher.test(mat_count, alternative = "greater")
    write.table(mat_count,paste0('table_fisher_test_TRA_',HRS_celltype,'_ther',ther,'.csv'),row.names = TRUE, col.names = TRUE, quote=FALSE, sep=',')
    
    pdf(paste0('barplot_fisher_test_TRA_',HRS_celltype,'_ther',ther,'.pdf'), height=10, width=15)
    barplot(as.matrix(mat_proportion), legend= TRUE, main=paste0("p-value= ",fish_1$p.value))
    dev.off()
    
    if (HRS_celltype=="HRS-CD4T"){
      total_CD4TwithTcell_TRA[i]=sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell
      total_Tcell_TRA[i]=sum(group_data_Tcell_greater_ther$n)/total_Tcell ###noHRST : Tcell
      print(sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell)
      print(sum(group_data_Tcell_greater_ther$n)/total_Tcell)
    }
    else{
      total_CD8TwithTcell_TRA[i]=sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell
      total_Tcell_TRA[i]=sum(group_data_Tcell_greater_ther$n)/total_Tcell ###noHRST : Tcell
      print(sum(group_ALL_data_withHRSandTcells_greater_ther$n)/total_cellCount_HRSwithTcell)
      print(sum(group_data_Tcell_greater_ther$n)/total_Tcell)
    }
    
    
  }
}


df_CD4TwithTcell_TRA=data.frame(c(3:7),t(total_CD4TwithTcell_TRA),t(total_Tcell_TRA))  
write.csv(df_CD4TwithTcell_TRA, 'threshold_clone_size_CD4T_TRA.csv')

df_CD4TwithTcell_TRB=data.frame(c(3:7),t(total_CD4TwithTcell_TRB),t(total_Tcell_TRB))  
write.csv(df_CD4TwithTcell_TRB, 'threshold_clone_size_CD4T_TRB.csv')

df_CD8TwithTcell_TRA=data.frame(c(3:7),t(total_CD8TwithTcell_TRA),t(total_Tcell_TRA))  
write.csv(df_CD8TwithTcell_TRA, 'threshold_clone_size_CD8T_TRA.csv')

df_CD8TwithTcell_TRB=data.frame(c(3:7),t(total_CD8TwithTcell_TRB),t(total_Tcell_TRB))  
write.csv(df_CD8TwithTcell_TRB, 'threshold_clone_size_CD8T_TRB.csv')



df=data.frame(c(3:7),t(total_CD4TwithTcell_TRA),t(total_Tcell_TRA))
colnames(df)<- c("x","CD4T_TCR","Tcell_TCR")
pdf('threshold_clone_size_CD4T_TRA.pdf')
ggplot()+
  geom_line(data=df,aes(y=CD4T_TCR,x= x,colour="CD4T_TCR"))+geom_point(data=df,aes(y=CD4T_TCR,x= x,colour="CD4T_TCR"),size=1 )+
  geom_line(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) + geom_point(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) +
  scale_color_manual(name = "groups", values = c("CD4T_TCR" = "red", "Tcell_TCR" = "green"))+
  labs(x = "threshold",y = "clone size (normalized)") + 
  ylim(0, 1)
dev.off()


df=data.frame(c(3:7),t(total_CD4TwithTcell_TRB),t(total_Tcell_TRB))
colnames(df)<- c("x","CD4T_TCR","Tcell_TCR")
pdf('threshold_clone_size_CD4T_TRB.pdf')
ggplot()+
  geom_line(data=df,aes(y=CD4T_TCR,x= x,colour="CD4T_TCR"))+geom_point(data=df,aes(y=CD4T_TCR,x= x,colour="CD4T_TCR"),size=1 )+
  geom_line(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) + geom_point(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) +
  scale_color_manual(name = "groups", values = c("CD4T_TCR" = "red", "Tcell_TCR" = "green"))+
  labs(x = "threshold",y = "clone size (normalized)") + 
  ylim(0, 1)
dev.off()


df=data.frame(c(3:7),t(total_CD8TwithTcell_TRA),t(total_Tcell_TRA))
colnames(df)<- c("x","CD8T_TCR","Tcell_TCR")
pdf('threshold_clone_size_CD8T_TRA.pdf')
ggplot()+
  geom_line(data=df,aes(y=CD8T_TCR,x= x,colour="CD8T_TCR"))+geom_point(data=df,aes(y=CD8T_TCR,x= x,colour="CD8T_TCR"),size=1 )+
  geom_line(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) + geom_point(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) +
  scale_color_manual(name = "groups", values = c("CD8T_TCR" = "red", "Tcell_TCR" = "green"))+
  labs(x = "threshold",y = "clone size (normalized)") + 
  ylim(0, 1)
dev.off()


df=data.frame(c(3:7),t(total_CD8TwithTcell_TRB),t(total_Tcell_TRB))
colnames(df)<- c("x","CD8T_TCR","Tcell_TCR")
pdf('threshold_clone_size_CD8T_TRB.pdf')
ggplot()+
  geom_line(data=df,aes(y=CD8T_TCR,x= x,colour="CD8T_TCR"))+geom_point(data=df,aes(y=CD8T_TCR,x= x,colour="CD8T_TCR"),size=1 )+
  geom_line(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) + geom_point(data=df,aes(y=Tcell_TCR,x= x,colour="Tcell_TCR")) +
  scale_color_manual(name = "groups", values = c("CD8T_TCR" = "red", "Tcell_TCR" = "green"))+
  labs(x = "threshold",y = "clone size (normalized)") + 
  ylim(0, 1)
dev.off()

