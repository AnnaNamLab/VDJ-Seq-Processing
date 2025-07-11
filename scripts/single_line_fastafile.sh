
### this command help to convert the fasta file to a merge the header and sequence of the fasta file in one line.
awk '/^>/&&NR>1{print "";}{printf "%s",/^>/ ? $0" " : $0}' TRUST_all_BCR_R2_annot.fa > singleline_TRUST_all_BCR_R2_annot.fa
