import argparse
import sys


def convert_to_fasta_format(detect_lines, max_contig_len, max_line_len):
    r""" Convert detect file to fasta format.

    Args:
        detect_lines (iter of strings): lines from a detect file
        max_contig_len (int): maximum possible length of contig, please use a few
          more bases than needed
        max_line_len (int): maximum length of fasta line.

    Returns:
        a string = the contents of the fasta file.

    Examples:
        >>> input1 = ['#blah\n', '>b contig1 1 3 fwd\n', '1\t0.1\tNNNNNN\n']
        >>> convert_to_fasta_format(input1, 3, 2)
        '>contig1\nNT\nN'
        >>> input2 = ['#blah\n', '>b contig2 1 3 rev\n', '1\t0.1\tNNNNNN\n', '3\t0.1\tNNNNNN\n']
        >>> convert_to_fasta_format(input2, 9, 5)
        '>contig2\nNNNNN\nNANA'
        >>> convert_to_fasta_format(input1 + input2 + input1 + input2, 12, 5)
        '>contig1\nNTNNN\nNNNNN\nNN\n>contig2\nNNNNN\nNANAN\nNN'
        >>> input2[1] = input2[1].replace("contig2", "contig1")
        >>> convert_to_fasta_format(input1 + input2, 9, 5)
        '>contig1\nNTNNN\nNANA'
    """

    contig_seq = dict()
    contig = ""
    shift = 0
    nucleotide = "T"

    # go through each record and populate contigs
    for line in filter(lambda x: not x.startswith("#"), detect_lines):

        if line.startswith(">"):

            _, contig, _, _, strand = line.rstrip()[1:].split(" ")

            if contig not in contig_seq:
                contig_seq[contig] = ["N"] * max_contig_len

            if strand == "rev":
                shift = 5
                nucleotide = "A"
            elif strand == "fwd":
                shift = 0
                nucleotide = "T"
            else:
                raise ValueError("Incorrect orientations")

        else:
            pos, _ = line.split("\t", maxsplit=1)
            contig_seq[contig][int(pos) + shift] = nucleotide

    return "\n".join(">" + k + "\n" + "\n".join(
        "".join(v[i: i + max_line_len]) for i in range(0, len(v), max_line_len)
        ) for k, v in contig_seq.items())


if __name__ == "__main__":

    desc = """    Recreate a fasta reference from a detect file. Only use if reference genome mentioned in
    the detect file is unavailable. Contig sequences have A or T or N. As this is a quick-and-dirty way to get a
    reference genome, all contigs have the same length, and are padded with Ns as need be.

    Input: Pipe in detect file contents.

    Output is to stdout and is in fasta format.

    Sample usage:
        cat sample.detect | python <programName.py> --maxContigLen 100000 > sample.fa
    """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--maxContigLen', type=int, required=True,
                        help='All contigs are of this length. Set to a few bases more than needed.')
    parser.add_argument('--maxLineLen', type=int, required=False,
                        help='(default: 80) max length of fasta line',
                        default=80)

    args = parser.parse_args()

    # ensure data is piped in
    if sys.stdin.isatty():
        parser.print_help(sys.stdout)
        raise NotImplementedError("Do not run interactively. Pipe in inputs.")

    # make fasta file
    print(convert_to_fasta_format(
        sys.stdin, args.maxContigLen, args.maxLineLen))
