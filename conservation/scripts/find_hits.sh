#!/bin/bash

# retrieve taxids of all strains under "Bacillus subtilis" species rank (taxid:1423) from NCBI database using NCBI tool
get_species_taxids.sh -t 1423 > 1423.taxids

# Taxonkit tool (https://bioinf.shenwei.me/taxonkit/) was used to obtain ranks and names of all strains
taxonkit list --show-rank --show-name --indent "    " --ids 1423 > 1423_strains.taxids

# obtain all hits with E-value < 1E-3 from NCBI nr database (downloaded on 7/8/2020) restricted to Bacillus subtilis strains taxids <1423 (1423.taxids)
blastp -db $HOME/blast/db/nr/nr -query ../data/query.fa -taxidlist ../meta/1423.taxids -out ../data/hits.tab -evalue 1e-3 -max_target_seqs 100000 -outfmt "7 qstart qend stitle sacc staxids evalue pident bitscore sstart send sseq" -num_threads 4
