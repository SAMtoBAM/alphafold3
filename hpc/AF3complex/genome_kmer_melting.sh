#!/bin/bash
# genome_melting.sh
# This script runs PAQman/samtools/bedtools on FASTA files transferred via HTCondor

# Exit on error
set -euo pipefail

# Parameters (could also pass via HTCondor submit if you want)
kmersize=8 

#loop over the genome

ls | grep "\.fa$\|.fasta$" | while read genome
do
    name=$(basename "$genome" | sed 's/\.[^.]*$//')

    meryl k=${kmersize} count ${genome} output ${name}_k8.meryl
    meryl print ${name}_k8.meryl > ${name}_k8.and_count.txt
    cat ${name}_k8.and_count.txt | shuf | cut -f1 | awk 'BEGIN{srand()} !/^$/ {if(rand() <= 0.5) print $0}' > ${name}_k8.random_half.txt 
    cat ${name}_k8.random_half.txt | awk '{print ">"$0"\n"$0}' > ${name}_k8.random_half.fa

done
