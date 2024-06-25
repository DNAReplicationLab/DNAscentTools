# DNAscentTools

Set of scripts to process files produced by
[DNAscent](https://github.com/MBoemo/DNAscent).

## Convert DNAscent detect to modBAM

Base modification probabilities can be added to BAM files using
[optional tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
Scripts here use a detect file and the corresponding reference
genome as input and produce a BAM file, with probabilities of
thymidine modification along each read specified using tags.

Please look at run_detect_to_modBAM.sh for an example on how to
perform a detect to modBAM conversion for different versions of DNAscent detect.
