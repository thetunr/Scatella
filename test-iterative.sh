#!/usr/bin/env bash

# In the terminal, execute `bash test-iterative.sh` to run all the test cases.

# TEST 0
python3 MUSCLE/iterative_alignment.py -f data/example.fasta -s 0 -d 10 -e 0 -o MUSCLE/output/example.fasta
echo -e "----------------------------------------------\n"

# TEST 1
# python3 MUSCLE/iterative_alignment.py -f data/sequences.fasta -s 0 -d 10 -e 0 -o MUSCLE/output/aligned.fasta
# echo -e "----------------------------------------------\n"