#!/bin/bash

# run script to test detect to modBAM functions

# stop execution if any command fails
set -e

# make a dummy fasta file
echo ">dummyI" > dummy.fa
echo "AGCTAGCTATCGTTTCTGTGAG" >> dummy.fa
echo ">dummyII" >> dummy.fa
echo "AGCTAGCTAGTCTCTAACGACCAA" >> dummy.fa
echo ">dummyIII" >> dummy.fa
echo "CCACACCACACCCACACACCCACACATCAAATCCACACCACACCACACCC" >> dummy.fa
echo "TGGGAGCCACCATAACGGCCTTATTG" >> dummy.fa

# make an index
samtools faidx dummy.fa

# test modBAM tools
python -m unittest test_modBAM_tools

# get rid of dummy files
rm dummy.fa
rm dummy.fa.fai