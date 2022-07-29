import argparse
import sys
from modBAM_tools import convert_dnascent_detect_to_modBAM_file, \
    convert_detect_into_detect_stream

if __name__ == "__main__":

    desc = """    Convert detect file to modBAM format.

    Input: Pipe in detect file contents.
    
    Output is to specified modBAM file.

    Sample usage:
        cat sample.detect | python <programName.py> --op sample.bam

    NOTE: * The initial comment block in the detect file must
            have a line '#Genome <filename>' with the reference
            genome in fasta format. Both the fasta file and its
            associated index called <filename>.fai must be
            available.
          * The default T modification is BrdU. Set a different
            modification using the --tag option (see below).
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--op', type=str, required=True,
                        help='output modBAM file')
    parser.add_argument('--tag', type=str, required=False,
                        help=('(default: 472552, which is BrdU) ChEBI code or one letter '
                              'code of base modification. If using a two-analogue '
                              'version of DNAscent, specify two codes using the format '
                              'code_1+code_2. e.g. 62903+472252. code_1 and _2 correspond '
                              'to the second and third columns in the detect file respectively'),
                        default='472552')
    parser.add_argument('--sam', required=False,
                        action='store_true',
                        help='(optional) make a sam file instead of a bam file',
                        default=False)

    args = parser.parse_args()

    # ensure data is piped in

    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # program info
    pgInfo = {
        'ID': 'convert_detect_to_modBAM',
        'PN': 'convert_detect_to_modBAM',
        'VN': 'v0.2'
    }

    # make modBAM/modSAM file
    convert_dnascent_detect_to_modBAM_file(
        convert_detect_into_detect_stream(sys.stdin),
        args.op, args.tag,
        not args.sam,
        pg_info=pgInfo)
