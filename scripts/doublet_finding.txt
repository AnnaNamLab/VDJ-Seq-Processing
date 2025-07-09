To remove the doublets based on BCR data, we have used the CDR3 raw file from TRUST4. 
This file contains all the contigs with number of reads that are supporting that. Each cell can be identified as singlet, or not-assigned, or doublet. The histogram of read counts for contigs per each cell defines the status of each cell.
For each cell: if a cell has only one contig, with any number of reads, that cell will be called singlet. If not, then using excess kurtosis, normalized entropy, and confidence level of VDJ assignment, a cell will be called doublet, not-assigned, or singlet. 
If there are equal or less than 7 reads in cells with 2 or more contigs, then the cell can be called singlet if either of these conditions are met: either the normalized entropy is small (less than 0.8) or the CDR3 Levenshtein distance in contigs is less than 0.2.
If there are equal or less than 7 reads in cells with 2 or more contigs, but none of the conditions for CDR3 similarity or small normalized-entorpy are met, then the cell is called not-assigned. If the excess kurtosis is positive value (not +Inf), then the cell is called singlet. 
If the excess kurtosis is negative value (not -Inf), then the cell is called doublet. If the kurtosis in +Inf or -Inf, and the number of reads is more than 7, then if the normalized entropy is more than 0.8, the cell is called doublet.
There is an extra step to check the doublets based on the conidence levelof VDJ assignments. Do obtain the the confidence level of VDJ assignments, the length of CDR1 and CDR2 are interfered. 
If one of CDR1 or CDR2 were missing or both of them are incomplete, then that contig's VDJ assignment has medium-confident.
If there is complete CDR1 and CDR2, then the contig's VDJ assignment is high_confident. If both CDR1 and CDR2 are missing then the contig's VDJ assignment is low_confident.
For any doublet, we have checked the confidence level of VDJ assignment to mae sure the doublets are coming from different origin. To do this step, for any doublet, the confidence level of VDJ-assignments for two contigs (with highest number of reads) are checked. 
If both contigs are high-confident, then the cell-status stays doublet. If the two contigs are both low-confident, or one of them is low-confident and another is medium-confident, then we calculate the similarity of the contigs CDR3 based on Levenshtein distance . 
If the CDR3s were similar (Levenshtein distance <= 0.2), then the cell status will change to singlet. But if the CDR3s were not similar, then the cell status stays doublet.
The singlets need to be intersected with HRS cells or B cells and we do this step per each chain (IGH and IGL/IGK). Then the clustering output from step 1 should be subsetted based on obtained singlets at this step.


<div align="center">
  <img src="doublet_refinement.png" width="500" alt="doublet refinement" />
</div>
