#!/usr/bin/env bash

# In the terminal, execute `bash run-iterative.sh` to run all alignments.

# ALIGNMENT 0
# python3 MUSCLE/iterative_alignment.py -f data/example.fasta -s 0 -d 10 -e 0 -o MUSCLE/output/example.fasta
# echo -e "----------------------------------------------\n"

# ALIGNMENT 1
python3 MUSCLE/iterative_alignment.py -f data/sequences.fasta -s 0 -d 430 -e 0 -o MUSCLE/output/aligned
echo -e "----------------------------------------------\n"