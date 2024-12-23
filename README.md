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


3)	Detecting the Dominant Clone

After removing the doublets from step 2, we identify the dominant clone per each chain. The dominant clone is the clone with the majority of HRS cells with same V/D/J assignment. We put no threshold for similarity of CDR3, since that causes to remove some of the HRS cells from the cancer clone, while we know that HRS cells are part of the clone. The dominant V/D/J on each chain will be recorded for the next step. 

4)	Correcting the V/D/J Assignments Errors

Form the output of the step 3 for detecting the dominant V/D/J, it happens that some of the HRS cells are carrying a different V/D/J comparing to dominant clone. One possible reason is about the errors during V/D/J gene identification after alignment of contigs to the sequences from IMGT database. To correct these errors, we extracted the contigs of HRS cells, and using the fasta file, which is one of the outcomes of TRUST4, started to extract new V/D/J assignments. From IgBlast, the top 3 V/D/J genes are extracted. Then for every HRS cell, if obtained IgBlast V/D/Js are the same of dominant V/D/J, then we correct the error from TRUST4. We, also explore the fasta file per each of the contigs, and if any of the top 3 V/D/Js are among the detected V/D/J for that contig, we correct that. That means per each contig any of the top 3 V/D/J obtained from IgBlast are neighbors of the dominant clone and that can still show that the cell belongs to the dominant clone. After correcting the V/D/Js per chain, we re-extracted all the contigs with similar V/D/J as the dominant clone and took it for generating the phylogenetic tree. 

5)	Building the Phylogenetic Tree

After extracting the HRS cells with similar V/D/J to the dominant clone and adding other HRS cells with not-similar V/D/J to the dominant clone, we use them to generate the phylogenetic tree. Phylogenetic tree is a diagram that shows the evolutionary relationships between CDR3/full-sequences of BCR. We used a maximum likelihood phylogenetic tree from Ape package in R [3]. The tree aims to put the sequences with highest similarity in the same clade based on genetic distances that are calculated using a substitution model. We used F81 substitution model for calculating the distances between nucleotides, since it is one of the wildly been used models for generating phylogenetic tree [4]. We generated the phylogenetic tree based on CDR3 and then annotated the cell states on each node on the tree. Also, we can identify the V/D/J assignments on the tree, or other BCR parameters to investigate their evolutionary process.



