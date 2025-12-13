#! /bin/bash -l

#SBATCH --partition=scu-cpu   # cluster-specific
#SBATCH --nodes=2
#SBATCH --ntasks=14
#SBATCH --job-name=tree_MAFFT
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sam4032@med.cornell.edu

sample=$1
echo $1

cd /athena/namlab/scratch/sam4032/CDR3_tree/singleCDR3/$sample

/athena/namlab/scratch/sam4032/CDR3_tree/igphyml/src/igphyml \
 -i /athena/namlab/scratch/sam4032/CDR3_tree/singleCDR3/$sample/sampled_unique_sequences_${sample}.fasta \
 -m GY \
 --partitions combined.part.txt \
 --run_id singleCDR3_$sample 


