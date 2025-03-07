import argparse
from modBAM_tools import convert_dnascent_detect_to_modBAM_file
from modBAM_tools_additional import convert_bed_to_detect_stream

if __name__ == "__main__":

    desc = """    Convert BED file to modBAM format.

    Sample usage:
        python <programName.py> --bed <bed_filename> --op sample.bam --fasta <filename> --tag T --base T 

    NOTE: The index associated with the fasta file called <filename>.fai must be available.
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--op', type=str, required=True,
                        help='output modBAM file')
    parser.add_argument('--bed', required=True, help='bed file to convert')
    parser.add_argument('--fasta', required=True, help='fasta file')
    parser.add_argument('--tag', type=str, required=True, help='ChEBI code or one letter code of base modification.')
    parser.add_argument('--base', type=str, required=True, help='base that is modified')
    parser.add_argument('--sam', required=False,
                        action='store_true',
                        help='(optional) make a sam file instead of a bam file',
                        default=False)

    args = parser.parse_args()

    # program info
    pgInfo = {
        'ID': 'convert_bed_to_modBAM',
        'PN': 'convert_bed_to_modBAM',
        'VN': 'v0.3'
    }

    # make modBAM/modSAM file
    convert_dnascent_detect_to_modBAM_file(
        iter(convert_bed_to_detect_stream(args.bed)),
        args.op, args.tag,
        not args.sam,
        pg_info=pgInfo,
        fasta=args.fasta,
        base=args.base,
        shift_rev_strand_pos=0)
    