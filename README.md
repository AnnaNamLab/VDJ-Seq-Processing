
# Single-Cell VDJ processing

The B cell receptor (BCR) and T cell receptor (TCR) are antigen-recognition receptors of the adaptive immune system. BCRs recognize antigens directly through membrane-bound immunoglobulins, whereas TCRs recognize peptide antigens presented by major histocompatibility complex (MHC) molecules. Both receptors achieve extraordinary diversity through V(D)J recombination during lymphocyte development, enabling precise and robust immune responses.

This repository contains code for:

### BCR processing steps:

1. Generating BCR Repertoire using TRUST4
1. Removing the Doublets using BCR sequences
1. Detecting the dominant clone
1. Refining the V/D/J assignments
1. Building the phylogenetic tree

The general scheme is illustrated below - More details about the steps are provided in the workflow folder.

<div align="left">
  <img
    width="800"
    height="800"
    alt="BCR workflow"
    src="workflow/BCR_flowchart.png"
  />
</div>

More details about the steps are provided in the workflow folder.


## Requirements

The following packages are required to be installed prior to running the pipeline:
1) TRUST4 (https://github.com/liulab-dfci/TRUST4) [1]
2) igblast (https://ncbi.github.io/igblast/cook/How-to-set-up.html) [2]
3) IgPhyML (https://igphyml.readthedocs.io/en/latest/install.html) [3, 4]
4) MAFFT (https://github.com/ GSLBiotech/mafft) 
5) R (v4.4.2)
    - <span style="color:red;">ggtree</span>
    - <span style="color:red;">ape</span>
    - <span style="color:red;">phangorn</span>
    - <span style="color:red;">DECIPHER</span>
    - <span style="color:red;">DECIPHER</span>
    
## Pipeline

### Generating BCR Repertoire using TRUST4

BCR repertoire is generated from BCR single cell fastq files using TRUST4 as follows:
```bash
chmod +x run_trust4_slurm.sh
sbatch run_trust4_slurm.sh <sample_name> <read1.fastq.gz> <read2.fastq.gz> <output_prefix>
```

TRUST4 generates multiple output files in `/path/to/trust4_output_dir/` including cdr3.out, which contains the CDR3 sequences.

Incomplete CDR3 sequences are then filtered out using the following command:

```bash
cd /path/to/trust4_output_dir/
python trust-cluster.py /path/to/trust4_output_dir/cdr3.out /path/to/trust4_output_dir/cluster_clone.tsv
```

The outcome file (cluster_clone.tsv) is used for more downstream analysis and generating the BCR phylogenetic tree.

### Removing the Doublets using BCR sequences

Doublets were identified based on the principle of allelic exclusion using the following command.
> `cell_metadata.csv`: four column table with cell barcode ("Full.cell_id") and cell type ("MainCelltype"), sub-cell type ("SubtypeName) and patient ID ("Patient")
Example `cell_metadata.csv`:
<details>
  <summary>See `cell_metadata.csv` example</summary>
  <table>
    <thead>
      <tr>
        <th>Full.cell_id</th>
        <th>MainCelltype</th>
        <th>SubtypeName</th>
        <th>Patient</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>AAACGGGCAAAGTGCG</td>
        <td>HRS</td>
        <td>State1</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AAACGGGCAGTCAGAG</td>
        <td>HRS</td>
        <td>State2</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AAACGGGGTATCAGTC</td>
        <td>HRS</td>
        <td>State1</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AAAGATGGTTGGTTTG</td>
        <td>HRS</td>
        <td>State3</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AAAGCAAAGTGGTAGC</td>
        <td>HRS</td>
        <td>State1</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AAATGCCTCTCTGTCG</td>
        <td>HRS</td>
        <td>State4</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AACGTTGGTAAGTTCC</td>
        <td>HRS</td>
        <td>State1</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AACTGGTCAAACGTGG</td>
        <td>HRS</td>
        <td>State2</td>
        <td>Sample1</td>
      </tr>
      <tr>
        <td>AACTGGTTCTGCTGTC</td>
        <td>HRS</td>
        <td>State1</td>
        <td>Sample1</td>
      </tr>
    </tbody>
  </table>
</details>

```bash
Rscript scripts/BCR/doublet_finding_BCR.R \
  --input /path/to/trust4_output_dir \
  --metadata /path/to/cell_metadata.csv \
  --sample Sample1 \
  --file /path/to/cdr3.out

```


### Detecting the dominant clone
The most expanded clone for each chain is identified, and the its V(D)J is recorded:
```bash
Rscript scripts/BCR/heatmap_before_VDJcorrection.R
```


### Refining the V/D/J assignments

For improved V(D)J assignment accuracy, we refine existing assignments from `TRUST4` by incorporating results from an additional `Igblast` method.

This step uses fasta file (`annot.fa`) from TRUST4 which contains multiple V(D)J assignments, and a list of contig IDs from cells of interest ("contigs_${sample}_${chain}.txt") from `/path/to/trust4_output_dir`.

<details>
  <summary>Example contig ID list: "contigs_Sample1_IGL.txt"</summary>
```
AAACGGGCAAAGTGCG_21556  
AAACGGGCAGTCAGAG_21980  
AAACGGGGTATCAGTC_21656  
AAAGATGGTTGGTTTG_22094  
AAAGCAAAGTGGTAGC_21648  
AAATGCCTCTCTGTCG_21695  
AACGTTGGTAAGTTCC_21525  
AACTGGTCAAACGTGG_22089  
AACTGGTTCTGCTGTC_21948  
```
</details>

First, we identify the V(D)J annotations using Igblast as follows:

```bash
chmod +x igblast_preprocess.sh
./igblast_preprocess.sh <sample> <clone_chain> <path to *annot.fa file as an output of TRUST4>
```

Here is an example: 
```bash
./igblast_preprocess.sh Sample1 IGL /path/to/trust4_output_dir/Sample1_T4_Output
```

Next, the V(D)J annotations are refined using the outputs of TRUST4 and Igblast with `refining_VDJ_BCR.R`

The above step results in `<TABLE_FILE_NAME>`. This table is cleaned to keep relevant columns and to re-evaluate the expanded clones after V(D)J refinement, using the following script `scripts/BCR/heatmap_after_VDJcorrection.R`

```bash
Rscript scripts/BCR/heatmap_after_VDJcorrection.R
```
The resulting table <span style="color:red">(<OUTPUT_FILENAME>)</span> contains only the complete CDR3 sequences and the corresponding refined V(D)J assignments.


### Building the phylogenetic tree

To generate the phylogenetic tree, the following steps are performed.

1. Data preparation (`IgPhyML_fasta_prepration.R`)
    - Input: <TABLE_FILE_NAME>
    - Output: "CDR3_preprocessed_MAFFT_aligned.fa"
    ```bash
    sbatch igphyml_run.sh Sample1
    ```

2. Tree construction using `IgPhyML`
    - `igphyml_run.sh` example:
      ```
      sample=Sample1
      path/to/igphyml \
      -i /path/to/CDR3_preprocessed_MAFFT_aligned.fa \
      -m GY \
      --run_id singleCDR3_${sample}
      ```

Finally, the resulting trees can be visualized and annotated using `scripts/BCR/igphyml_downStream.R`. This code also allows to assess the timing based on branch length and the CDR3 divergence based on the Levenshtein distances of the aligned CDR3 sequences.


## TCR Analysis

For TCR: 
```bash
Rscript scripts/BCR/doublet_finding_TCR.R \
  --input /path/to/trust4_output_dir \
  --metadata /path/to/cell_metadata.csv \
  --sample HL1 \
  --file /path/to/HL1_cdr3.out

```



## References
----------

1. Song L, Cohen D, Ouyang Z, Cao Y, Hu X, Liu XS. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nat Methods. 2021;18(6):627-630. doi:10.1038/s41592-021-01142-2.

2. Ye J, Ma N, Madden TL, Ostell JM. IgBLAST: an immunoglobulin variable domain sequence analysis tool. Nucleic Acids Res. 2013;41(W34–40). doi:10.1093/nar/gkt382.

3. Hoehn KB, Lunter G, Pybus OG. A Phylogenetic Codon Substitution Model for Antibody Lineages. Genetics. 2017;206(1):417–427. https://doi.org/10.1534/genetics.116.196303

4. Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SH. Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination. bioRxiv. 2019. https://doi.org/10.1101/558825
5. Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002 Jul 15;30(14):3059-66. doi: 10.1093/nar/gkf436. PMID: 12136088; PMCID: PMC135756.
EOF



