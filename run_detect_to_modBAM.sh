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

# if using the two-analogue version of DNAscent, use code_1+code_2
# where code_1 and code_2 correspond to the second and third column
# of the detect file respectively. For e.g.:
# python convert_detect_to_modBAM.py --op $modBAMPrefix.bam --tag g+e

# index and sort files
samtools sort -o $modBAMPrefix.sorted.bam $modBAMPrefix.bam
samtools index $modBAMPrefix.sorted.bam