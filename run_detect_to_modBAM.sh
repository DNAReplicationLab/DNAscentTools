#!/bin/bash

# run script to convert detect to modBAM

# stop execution if any command fails
set -e

rootDir=/path/to/detect/dir/
detectFile=$rootDir/sample.detect
modBAMPrefix=$detectFile.mod

# set the script information
infoStr="#a sample comment"

# perform detect to modBAM conversion.
# NOTE: user can replace specified tag T, which means
# a generic thymidine modification. default tag
# is the CheBI code of BrdU and is used when no tag
# is specified.
sed "1i$infoStr" $detectFile |\
    python convert_detect_to_modBAM.py --op $modBAMPrefix.bam --tag T

# index and sort files
samtools sort -o $modBAMPrefix.sorted.bam $modBAMPrefix.bam
samtools index $modBAMPrefix.sorted.bam