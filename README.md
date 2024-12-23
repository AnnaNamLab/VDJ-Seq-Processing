VDJ processing

VDJ (Variable, Diversity, and Joining) processing is a key mechanism in the generation of diverse B cell receptors (BCRs), which are essential for the adaptive immune response. The BCR is an antibody anchored on the surface of B cells, recognizing specific antigens through highly variable regions, generated during B cell development. The generation of a diverse B cell receptor (BCR) repertoire is fundamental to the immune system's ability to recognize and respond to a vast array of pathogens. This repository contains all codes related to generating BCR/TCR sequencing, and all downstream analysis.

Method)

Building the phylogenetic tree from BCR data has 5 main steps: generating BCR repertoire using TRUST4, removing the doublets, detecting the dominant clone, correcting the VDJ assignments errors in Hodgkin and Reed-Sternberg (HRS) cells, and building the phylogenetic tree

<img width="545" alt="BCR_steps" src="https://github.com/user-attachments/assets/2d5769f9-1764-46a4-9117-1ef7365244b8" />
1) 1)	Generating BCR Repertoire Using TRUST4

To reconstruct the BCR repertoire, we have applied TRUST4. TRUST4 is an open-source algorithm for reconstruction of immune receptor repertoire in αβ/γδ B cells and T cells from bulk and single-cell RNA-seq data [1]. TRUST4 can assemble full sequences based on VDJ assignments. This algorithm extracts the good TCR/BCR candidate reads from BAM or FATSQ files with fast speed and has a de-novo approach for assembling the candidate reads. It implements a greedy strategy by aligning the candidate reads to an existing contig. For contigs annotation, TRUST4 uses ImmunoGeneTics (IMGT) database [2]. TRUST4 has multiple output files including the extracted contigs and corresponding nucleotide weight, a fasta file for annotation of the consensus assembly, and files with the focus on the CDR3 and constructed full sequences. Also, TRUST4 provides an extra step for clustering the contigs based on similarity of VDJs and CDR3, and during this step they remove any partial or incomplete CDR3. This clustering file is the output we used for downstream. More information and the github for TRUST4 is available in: https://github.com/liulab-dfci/TRUST4

2)	Removing the Doublets

To remove the doublets based on BCR data, we have used the CDR3 raw file from TRUST4. This file contains all the contigs with number of reads that are supporting that. Each cell can be identified as singlet, or noise, or doublet. The histogram of read counts for contigs per each cell defines the status of each cell. For each cell: if the maximum number of supporting reads among all contigs is less than 2, then that cell is labels ‘noise’. If more than 55% of the reads are associated to one contig, then that cell is labeled as ‘singlet’. 
If more than 55% of the reads are not associated to one contig, then two conditions appear based on the difference of read counts in the top two contigs with highest number of reads: if the difference of read counts in top two contigs is more than a threshold (here we used 0.35), then the cell is labeled as ‘singlet’. But, if the difference of read counts in top two contigs is less than that threshold, then the cell is labeled as ‘doublet’. The singlets need to be intersected with HRS cells or B cells and we do this step per each chain (IGH/IGL/IGK). Then the clustering output from step 1 should be subsetted based on obtained singlets at this step. 
![image](https://github.com/user-attachments/assets/9315e448-a36e-44e3-af49-e4869fcbdf76)
