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
# There are variations of this command depending on which version of DNAscent produced the .detect file; please
# read below for more information.
# In every .detect file, there follows a table containing per-base information on the read after the header line
# reporting the read and its alignment information (header lines start with ">").
# In V2, the columns are position, prob BrdU, and a six-mer sequence.
# In V3, the columns are position, prob EdU, prob BrdU.
# In V4, the columns are position, prob EdU, prob BrdU, and a 9-mer sequence.
# Position refers to position on the reference genome where we make a call on a T (or on an A in the case of a reverse
# read), and prob EdU and prob BrdU are the probabilities of the base being an EdU or a BrdU, respectively.
# Unfortunately DNAscent detect offsets the position column by different amounts depending on the version:
# In V2, on a forward read, the offset is zero, but on a reversed read the offset is -5.
# In V3, the offset is zero on both forward and reversed reads.
# In V4, the offset is -4 on a forward read and +4 on a reversed read.
# (By 'offset' we mean if the probability call is at base i, then DNascent reports i + offset in the position column).
# This tool was first written for V2, so adjustments are made for V3 and V4 in the command below.
# We do not use the six-mer or nine-mer sequence information in the detect -> modBAM conversion (as of this date).
# (There is an optional construct_reference_from_detect.py which uses the 6-mer information but this is
#  not part of a standard workflow. It is not advisable to run this script as it attempts to create
#  a reference from the 6-mer in a v2 .detect file. It should only be used when a reference is missing)

# For v2:
# =======
# NOTE: user can replace specified tag T, which means
# a generic thymidine modification. default tag
# is the CheBI code of BrdU and is used when no tag
# is specified.
sed "1i$infoStr" $detectFile |\
    python convert_detect_to_modBAM.py --op $modBAMPrefix.bam --tag T

# For v3:
# =======
sed "1i$infoStr" $detectFile |\
  awk \
    -v IFS="\t" -v OFS="\t" \
      '{
      if ($1 ~ /^#/) {
          print $0
      } else if ($1 ~ /^>/) {
          print $0
      } else {
          print $1, $2, $3, "NNNNNN"
      }
  }' | python convert_detect_to_modBAM.py  --op $modBAMPrefix.bam --tag E+T --undo-shift-rev-strand-pos

# If you are only interested in one analogue, you can just use --tag T (or an appropriate numeric code) for the analogue
# and remove $2 or $3 from the print statement above ($2 is the second column that records probEdU,
# $3 is the third column that records probBrdU).

# E and T are not good tags for EdU and BrdU, respectively, because T means "ambiguous thymidine modification" and E
# is not a valid tag. The best tags for EdU and BrdU would be their CheBI codes.
# Unfortunately, there is no CheBI code for EdU (to our best knowledge as of this date),
# but there is one for BrdU (CHEBI:472552).
# Let's say for the purposes of argument that there is an EdU code, say 000000.
# Then, you can invoke the script using --tag 000000+472552 as the argument and this would be the best way to tag the
# EdU and BrdU modifications. As such an option is not available, we leave it to the user to decide what to do.

# For v4:
# =======
sed "1i$infoStr" $detectFile |\
  awk \
    -v IFS="\t" -v OFS="\t" -v flag=0 \
      '{
      if ($1 ~ /^#/) {
          print $0
      } else if ($1 ~ /^>/) {
          print $0
          if ($0 ~ /rev$/) {
              flag = 5
          } else {
              flag = 0
          }
      } else {
          if (flag == 0) {
              print $1 + 4, $2, $3, substr($4, 5, 5) "N"
          } else if (flag == 5) {
              print $1 - 4, $2, $3, "N" substr($4, 1, 5)
          }
      }
  }' | python convert_detect_to_modBAM.py  --op $modBAMPrefix.bam --tag E+T --undo-shift-rev-strand-pos

# The same comments we wrote for v3 regarding the tag information and one analogue mode apply here as well.

# index and sort files
samtools sort -o $modBAMPrefix.sorted.bam $modBAMPrefix.bam
samtools index $modBAMPrefix.sorted.bam