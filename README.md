
# VDJ processing

VDJ (Variable, Diversity, and Joining) processing is a key mechanism in the generation of diverse B cell receptors (BCRs), which are essential for the adaptive immune response. The BCR is an antibody anchored on the surface of B cells, recognizing specific antigens through highly variable regions, generated during B cell development. The generation of a diverse B cell receptor (BCR) repertoire is fundamental to the immune system's ability to recognize and respond to a vast array of pathogens. This repository contains all codes related to generating BCR/TCR sequencing, and all downstream analysis.

## Requirements

1) TRUST4 (https://github.com/liulab-dfci/TRUST4)
2) igblast (https://ncbi.github.io/igblast/cook/How-to-set-up.html)
3) PhyloPart_v2.1 (https://sourceforge.net/projects/phylopart/)

## Method

Building the phylogenetic tree from BCR data has 5 main steps: generating BCR repertoire using TRUST4, removing the doublets, detecting the dominant clone, correcting the VDJ assignments errors in Hodgkin and Reed-Sternberg (HRS) cells, and building the phylogenetic tree

<div align="center">
  <img hwidth="545" alt="BCR_steps" src="https://github.com/user-attachments/assets/2d5769f9-1764-46a4-9117-1ef7365244b8" />
</div>
<p></p>
<p>1)  Generating BCR Repertoire Using TRUST4</p>
<p>2)	Removing the Doublets</p>
<p>3)	Detecting the Dominant Clone</p>
<p>4)	Correcting the V/D/J Assignments Errors</p>
<p>5)	Building the Phylogenetic Tree</p>


