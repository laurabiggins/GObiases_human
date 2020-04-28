#!/bin/bash
index=$1;
#echo $index;
python ./generate_genelist_from_random_positions_human.py --output gene_lists_human/closest_genes_$index.txt;
#echo "closest_genes_$1.txt"
#echo "closest_genes_{$1}.txt"
