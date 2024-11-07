#!/bin/bash

#SBATCH --partition=scu-cpu   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=igblast
#SBATCH --time=48:00:00   # HH/MM/SS
#SBATCH --mem=512G   # memory requested, units available: K,M,G,T

#fastafile=./test/out_FR2_annot.fa

#cd /athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin


sample="HL8"
cd /athena/namlab/scratch/sam4032/HL8_s1s2/HL8_T4_Output_“2024_01_04”
sed '/^>/ s/ .*//' TRUST_all_BCR_R2_annot.fa > fix_my_sample.fa
grep -w -A 1 -f  contigs_HL8.txt fix_my_sample.fa --no-group-separator > sub_fix_my_sample.fa

/athena/namlab/scratch/sam4032/immcantation-master/scripts/ncbi-igblast-1.17.0/bin/igblastn -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGHV_clean.fasta -num_alignments_V 3 -germline_db_D /home/sam4032/share/igblast/database/imgt_human_IGHD_clean.fasta -num_alignments_D 3 -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGHJ_clean.fasta -num_alignments_J 3 -query sub_fix_my_sample.fa -organism human -outfmt 19 -out ${sample}_igblast_output2.txt



igblastn -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGHV_clean.fasta -num_alignments_V 3 -germline_db_D /home/sam4032/share/igblast/database/imgt_human_IGHD_clean.fasta -num_alignments_D 3 -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGHJ_clean.fasta -num_alignments_J 3 -query sub_fix_my_sample.fa -organism human -outfmt 7 -out ${sample}_igblast_output2_7.txt


