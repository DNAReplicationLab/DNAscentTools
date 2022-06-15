#!/bin/bash

# run script to convert detect to modBAM

# stop execution if any command fails
set -e

rootDir=/path/to/detect/dir/
detectFile=$rootDir/sample.detect
modBAMPrefix=$detectFile.mod

# set the script information
infoStr="#a sample comment"

# perform conversion
cat $detectFile |\
    sed "1i$infoStr" |\
    python convert_detect_to_modBAM.py --op $modBAMPrefix.bam

# index and sort files
samtools sort -o $modBAMPrefix.sorted.bam $modBAMPrefix.bam
samtools index $modBAMPrefix.sorted.bam