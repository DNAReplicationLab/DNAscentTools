import numpy as np
import itertools
import pysam
try:
    from modBAM_tools_additional import get_raw_data_from_modBAM, ModBamRecordProcessor
except ImportError:
    from .modBAM_tools_additional import get_raw_data_from_modBAM, ModBamRecordProcessor


def get_gaps_in_base_pos(pos, seq, base):
    """ Count occurrences of base in seq in open interval (pos_i, pos_{i+1})

    Args:
        pos (list of int): Positions of interest, in ascending order, no duplicates, and entries >= -1 & <= len(seq)
        seq (str): Sequence of interest
        base (str): one-character string

    Returns:
        list of ints with length 1 shy of length of pos.

    Raises:
        ValueError: If positions are not in ascending order, or if positions are inappropriate, or if positions do not match the base.

    Examples:
        >>> get_gaps_in_base_pos([1, 5], 'ATATATA', 'T')
        [1]
        >>> get_gaps_in_base_pos([1, 5], 'ATTTATA', 'T')
        [2]
        >>> get_gaps_in_base_pos([-1, 1, 5, 9, 14, 17], 'GGXXXGXXXGXXGXGXG', 'G')
        [1, 0, 0, 1, 1]
        >>> get_gaps_in_base_pos([1, 5, 9, 14], 'XTXXXGXXXGXXGXGXG', 'G')
        Traceback (most recent call last):
            ...
        ValueError: Positions of interest must be G
        >>> get_gaps_in_base_pos([1, 5, 9, 18], 'XGXXXGXXXGXXGXGXG', 'G')
        Traceback (most recent call last):
            ...
        ValueError: Inappropriate entries
        >>> get_gaps_in_base_pos([-2, 5, 9, 14], 'XGXXXGXXXGXXGXGXG', 'G')
        Traceback (most recent call last):
            ...
        ValueError: Inappropriate entries
        >>> get_gaps_in_base_pos([5, 1, 9, 14], 'XGXXXGXXXGXXGXGXG', 'G')
        Traceback (most recent call last):
            ...
        ValueError: Positions must be in ascending order
    """
    # raise error if first and last entries are inappropriate
    if pos[0] < -1 or pos[-1] > len(seq):
        raise ValueError("Inappropriate entries")

    # check pos is in ascending order
    if any(pos[k] - pos[k - 1] <= 0 for k in range(1, len(pos))):
        raise ValueError("Positions must be in ascending order")

    # raise error if any entry of seq has a different base at pos
    if not all(seq[m] == base for m in filter(lambda x: -1 < x < len(seq), pos)):
        raise ValueError(f'Positions of interest must be {base}')

    return [seq[k[0] + 1: k[1]].count(base) for k in
            zip(pos[:-1], pos[1:])]


def convert_data_per_T_to_modBAM_fmt(seq_data):
    """ Convert data per thymidine to modBAM format

    Using the SAM specification at
    https://samtools.github.io/hts-specs/SAMtags.pdf

    Args:
        seq_data (list of floats): BrdU probabilities

    Returns:
        list of probs scaled to 0 to 255

    Examples:
        >>> convert_data_per_T_to_modBAM_fmt([0.0, 0.5, 1.0])
        [0, 128, 255]
        >>> convert_data_per_T_to_modBAM_fmt([0.1, 0.2, 0.3])
        [25, 51, 76]
        >>> convert_data_per_T_to_modBAM_fmt([0.999, 0.001])
        [255, 0]
        >>> convert_data_per_T_to_modBAM_fmt([0.123, 0.456, 0.789])
        [31, 116, 201]
        >>> convert_data_per_T_to_modBAM_fmt([1.1, -0.1])
        [255, 0]
    """

    # convert entries to integers b/w 0 and 255, both included
    return [int(k) for k in
            np.clip(np.floor([r * 256 for r in seq_data]), 0, 255)]


def convert_detect_into_detect_stream(detect_obj, switch_2_and_3=False):
    """ Creates iterator that goes one detect record at a time. First item
    is a dictionary with keys comments and refFasta.

    Code adapted frm https://dnascent.readthedocs.io/en/latest/cookbook.html

    Args:
        detect_obj (iterable): each line of detect file
        switch_2_and_3 (bool): (default False) switch columns 2 and 3

    Yields:
        dict with keys readID, refContig, refStart, refEnd, strand,
          posOnRef, probBrdU, probEdU, sixMerOnRef. Meanings are as in the detect file.
          NOTE: First-ever item has keys comments (list of str) and
            refFasta (str).

    """
    col2, col3 = (1, 2) if not switch_2_and_3 else (2, 1)

    current_entry = {"comments": []}

    # split by '>' character, or per detection in other words
    for key, sub_iter in itertools.groupby(detect_obj,
                                           lambda x: x.startswith('>')):

        if key:

            # if header line, initialize record
            split_line = next(sub_iter).rstrip().split()

            # store header information
            current_entry = {
                'readID': split_line[0][1:],
                'refContig': split_line[1],
                'refStart': int(split_line[2]),
                'refEnd': int(split_line[3]),
                'strand': split_line[4],
                'posOnRef': [],
                'probBrdU': [],
                'probEdU': [],
                'sixMerOnRef': []
            }

        else:

            is_three_line_entry = False
            is_four_line_entry = False
            is_comment = False

            for line in sub_iter:

                # in dnascent detect files, a few comment lines
                # appear before the header. save these.
                # of particular interest is the line that starts
                # "#Genome" and contains the reference fasta file.
                if line.startswith('#'):

                    # store comments
                    is_comment = True
                    current_entry["comments"].append(line.rstrip())

                    # record fasta reference file
                    if line.startswith('#Genome'):
                        current_entry["refFasta"] = line.replace("#Genome ", "").rstrip()

                    continue

                # split data lines by whitespace
                split_line = line.rstrip().split()

                # store data
                # skip malformed lines if they exist
                if len(split_line) == 3:
                    current_entry['posOnRef'].append(int(split_line[0]))
                    current_entry['probBrdU'].append(float(split_line[col2]))
                    current_entry['sixMerOnRef'].append(split_line[col3])
                    is_three_line_entry = True
                elif len(split_line) == 4:
                    current_entry['posOnRef'].append(int(split_line[0]))
                    current_entry['probEdU'].append(float(split_line[col2]))
                    current_entry['probBrdU'].append(float(split_line[col3]))
                    current_entry['sixMerOnRef'].append(split_line[3])
                    is_four_line_entry = True
                else:
                    continue

            # return current record
            if is_three_line_entry or is_four_line_entry or is_comment:
                yield current_entry
            else:
                raise ValueError("Malformed detect file")


def convert_dnascent_detect_to_modBAM_file(detect_stream, filename,
                                           code="472552", bam_file=True, shift_rev_strand_pos=5, pg_info=None,
                                           fasta="", base='T'):
    """ Convert dnascent detect data to modBAM format

    NOTE: About dnascent reverse offset issue and the param shift_rev_strand_pos.
        In the DNAScent docs 
        (https://dnascent.readthedocs.io/en/latest/detect.html), an
        example of a reversed read is given.
        
        >c6785e1f-10d2-49cb-8ca3-e8d48979001b chrXIII 74003 81176 rev
        74010   0.012874        TCTCTA
        74011   0.012428        CTCTAA
        74014   0.016811        TAACGA
        74017   0.013372        CGACCA
        74018   0.013836        GACCAA

        So, although the leftmost column gives a reference coordinate, say i,
        the probability is for the base at i+5. Hence, we need to
        increment all pos on a reverse strand by 5, accomplished via the
        parameter shift_rev_strand_pos. If DNAscent output changes in future
        versions, set the parameter suitably.
    
    Args:
        detect_stream (iter of dicts): see below.
          * first entry must have a key comments = a list of comments
            from a dnascent detect file, and a key refFasta = location of
            fasta file. comments can be an empty list.
          * every other entry must have the keys readID (str), refStart (int), 
            refEnd (int), refContig (str), posOnRef (list of ints), 
            probBrdU (list of floats) which correspond to data per detect from
            the detect file.
        filename (str): name of output file
        code (str): default 472552, modification code (ChEBI code or 1 letter). If using a two-analogue
            version of DNAscent, specify two codes using the format code_1+code_2. e.g. 62903+472252.
            code_1 and _2 correspond to the second and third columns in the detect file respectively.
        bam_file (bool): default T, T/F = bam/sam file
        shift_rev_strand_pos (int): default 5. See note above, in fn docstr.
        pg_info (dict): (default None = no info). Info abt program, must
          have key 'ID'.
        fasta (str): (default ""). If supplied, this reference genome overrides
          reference genome in detect file
        base (str): (default 'T'). Base that is modified.

    Returns:
        None
    """

    # get comments and fasta file from detect
    # start making header to modBAM file
    init_info = next(detect_stream)
    header = {'CO': [
        'comment: ' + k for k in init_info['comments']
    ]
    }

    if fasta:
        ref_fasta_file = pysam.FastaFile(filename=fasta)

        # remove existing reference genome in comments if need be and replace with new reference genome
        for comment_line in header['CO']:
            if comment_line.startswith('comment: #Genome'):
                header['CO'].remove(comment_line)

        header['CO'].append(f'comment: #Genome {fasta}')

    elif 'refFasta' in init_info:
        ref_fasta_file = pysam.FastaFile(filename=init_info['refFasta'])
    else:
        raise ValueError("Fasta reference genome not found!")

    # load names of contigs, their order, and lengths from fasta file
    ref_name_length_map = dict(zip(ref_fasta_file.references,
                                   ref_fasta_file.lengths))
    ref_name_id_map = dict(zip(ref_fasta_file.references,
                               itertools.count()))

    # add contig and length information to the header
    header['SQ'] = [
        {'LN': v, 'SN': k} for k, v in ref_name_length_map.items()
    ]

    # add program information to the header
    if pg_info:
        header['PG'] = [pg_info]

    # iterate through each detect entry and add to modBAM/modSAM file
    with pysam.AlignmentFile(filename, "wb" if bam_file else "w",
                             header=header) as fb:

        for oneDetect in detect_stream:

            # set up reverse strand apparatus
            fwd_strand = oneDetect['strand'] == 'fwd'

            def fn_rev_if_needed(x):
                return x if fwd_strand else reversed(x)

            # reconstruct detect header but with underscores
            detect_header = "_".join(
                str(oneDetect[k]) for k in
                ['readID', 'refContig', 'refStart', 'refEnd', 'strand']
            )

            # get sequence from reference and convert to upper case
            seq = ref_fasta_file.fetch(oneDetect['refContig'],
                                       oneDetect['refStart'],
                                       oneDetect['refEnd']).upper()

            # check that length of retrieved sequence and DNAscent header match
            if not len(seq) == oneDetect['refEnd'] - oneDetect['refStart']:
                raise ValueError("Problems in retrieved sequence length")

            # create modBAM entry
            seg = pysam.AlignedSegment()
            seg.query_name = oneDetect['readID']
            seg.query_sequence = seq
            seg.flag = 0 if fwd_strand else 16
            seg.reference_id = ref_name_id_map[oneDetect['refContig']]
            seg.reference_start = oneDetect['refStart']
            seg.mapping_quality = 255
            seg.cigar = [(0, len(seq))]
            seg.template_length = len(seq)
            seg.query_qualities = pysam.qualitystring_to_array("")

            # scale probabilities from 0 to 255
            prob_frm0_to255 = convert_data_per_T_to_modBAM_fmt(
                oneDetect['probEdU'] + oneDetect['probBrdU'])

            # get number of skipped Ts b/w each T for which data is available
            # NOTE: the fn get_forward_sequence() returns seq if the reversed
            #   flag is not set, otherwise it returns the reversed complement.
            if fwd_strand:
                rel_pos = [k - oneDetect['refStart'] for k in
                           oneDetect['posOnRef']]
            else:
                rel_pos = [oneDetect['refEnd'] - 1 - (k + shift_rev_strand_pos)
                           for k in reversed(oneDetect['posOnRef'])]

            thym_coords = list(itertools.chain([-1], rel_pos))

            gaps_in_t = get_gaps_in_base_pos(thym_coords,
                                             seg.get_forward_sequence(), base)

            str_gaps_in_t = ",".join(str(int(k)) for k in gaps_in_t)

            # function to convert UUID to int
            def uuid_first7_char_to_int(x):
                return int(x[0:7], 16)

            # set modification positions and values
            # ChEBI code for BrdU is 472552
            # check if multiple analogues are used
            if '+' not in code:
                if len(oneDetect['probEdU']) == 0:
                    seg.set_tag("MM", f"{base}+{code}?,{str_gaps_in_t};", "Z")
                else:
                    raise ValueError('Please specify two codes if using a two-analogue version of DNAscent')
            else:
                code_1, code_2 = fn_rev_if_needed(code.split('+'))
                seg.set_tag("MM", f"{base}+{code_1}?,{str_gaps_in_t};" + f"{base}+{code_2}?,{str_gaps_in_t};", "Z")

            seg.set_tag("ML", list(fn_rev_if_needed(prob_frm0_to255)))
            seg.set_tag("XR", uuid_first7_char_to_int(oneDetect['readID']), "i")
            seg.set_tag("XA", f"{detect_header}", "Z")

            # store modBAM entry
            fb.write(seg)

    ref_fasta_file.close()


def get_mod_counts_per_interval(mod_bam_file, intervals, base,
                                mod_code, low_thres, high_thres):
    """ Get modification counts given reference-coordinate intervals

    Args:
        mod_bam_file (str): path to modBAM file
        intervals (iter): each entry has 3 elements = contig,start,end
        base (str): 'A', 'G', 'C' or 'T'
        mod_code (str): modification code
        low_thres (float): probability below which base marked as unmodified
        high_thres (float): probability above which base marked as modified

    Returns:
        Tuple of 4 iterators unmodified counts in rev and forward strand,
            modified counts in reverse and forward strands
    """

    if not low_thres == high_thres:
        raise NotImplementedError("Low and high thresholds must be equal")

    # set up mod bam parser
    mod_bam_parser = ModBamRecordProcessor(low_thres, mod_code, allow_non_na_mode=True, base=base)

    # set up output counts
    unmod_counts_fwd_list = []
    unmod_counts_rev_list = []
    mod_counts_fwd_list = []
    mod_counts_rev_list = []

    # iterate through the intervals
    for contig, start, end in intervals:

        unmod_counts_fwd = 0
        unmod_counts_rev = 0
        mod_counts_fwd = 0
        mod_counts_rev = 0

        # find the records
        alignments = pysam.view(mod_bam_file, f"{contig}:{start}-{end}")

        # iterate through records
        if len(alignments) > 0:
            for x in alignments.splitlines():
                if x:
                    # get modification counts and assign them suitably
                    mod_bam_parser.process_modbam_line(x)
                    count_mod, count_unmod, _ = mod_bam_parser.count_bases(mask_ref_interval=(start, end))
                    if mod_bam_parser.is_rev:
                        unmod_counts_rev += count_unmod
                        mod_counts_rev += count_mod
                    else:
                        unmod_counts_fwd += count_unmod
                        mod_counts_fwd += count_mod

        # store counts
        unmod_counts_fwd_list.append(unmod_counts_fwd)
        unmod_counts_rev_list.append(unmod_counts_rev)
        mod_counts_fwd_list.append(mod_counts_fwd)
        mod_counts_rev_list.append(mod_counts_rev)

    return (iter(unmod_counts_rev_list), iter(unmod_counts_fwd_list), iter(mod_counts_rev_list),
            iter(mod_counts_fwd_list))


def get_read_data_from_modBAM(mod_bam_file: str, read_id: str,
                              contig: str, start: int, end: int):
    """ Get modification probabilities corresponding to coords on a read

    Args:
        mod_bam_file: path to modBAM file
        read_id: requested read ID
        contig: contig on reference genome
        start: start pos on ref genome, 0-based
        end: end pos on ref genome, 0-based

    Returns:
        Iterator w each entry = (ref pos, probability)
    """

    return map(lambda x: (x[1], x[2]), get_raw_data_from_modBAM(mod_bam_file, contig, start, end,
                                                                base='T', code='T', read_id=read_id,
                                                                allow_multiple_entries_input_read_id=False))
