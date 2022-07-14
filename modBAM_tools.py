import numpy as np
import itertools
import pysam
from modbampy import ModBam


def get_gaps_in_base_pos(pos, seq, base):
    """ Count occurrences of base in seq in open interval (pos_i, pos_{i+1})

    Args:
        pos (list of int): Positions of interest, entries >= 0 & < len(seq)
        seq (str): Sequence of interest
        base (str): one-character string

    Returns:
        list of ints with length 1 shy of length of pos.
    """
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
    """

    # convert entries to integers b/w 0 and 255, both included
    return [int(k) for k in
            np.clip(np.floor([r * 256 for r in seq_data]), 0, 255)]


def convert_detect_into_detect_stream(detect_obj):
    """ Creates iterator that goes one detect record at a time. First item
    is a dictionary with keys comments and refFasta.

    Code adapted frm https://dnascent.readthedocs.io/en/latest/cookbook.html

    Args:
        detect_obj (iterable): each line of detect file

    Yields:
        dict with keys readID, refContig, refStart, refEnd, strand,
          posOnRef, probBrdU, sixMerOnRef. Meanings are as in the detect file.
          NOTE: First-ever item has keys comments (list of str) and
            refFasta (str).

    """

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
                'sixMerOnRef': []
            }

        else:

            for line in sub_iter:

                # in dnascent detect files, a few comment lines
                # appear before the header. save these.
                # of particular interest is the line that starts
                # "#Genome" and contains the reference fasta file.
                if line.startswith('#'):

                    # store comments
                    current_entry["comments"].append(line.rstrip())

                    # record fasta reference file
                    if line.startswith('#Genome'):
                        current_entry["refFasta"] = line.replace("#Genome ", "").rstrip()

                    continue

                # split data lines by whitespace
                split_line = line.rstrip().split()

                # skip malformed lines if they exist
                if not len(split_line) == 3:
                    continue

                # store data
                current_entry['posOnRef'].append(int(split_line[0]))
                current_entry['probBrdU'].append(float(split_line[1]))
                current_entry['sixMerOnRef'].append(split_line[2])

            # return current record
            yield current_entry


def convert_dnascent_detect_to_modBAM_file(detect_stream, filename,
                                           code="472552", bam_file=True, shift_rev_strand_pos=5, pg_info=None):
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
        code (str): default 472552, modification code (ChEBI code or 1 letter)
        bam_file (bool): default T, T/F = bam/sam file
        shift_rev_strand_pos (int): default 5. See note above, in fn docstr.
        pg_info (dict): (default None = no info). Info abt program, must
          have key 'ID'.

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

    ref_fasta_file = pysam.FastaFile(filename=init_info['refFasta'])

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

            # get sequence from reference
            seq = ref_fasta_file.fetch(oneDetect['refContig'],
                                       oneDetect['refStart'],
                                       oneDetect['refEnd'])

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
                oneDetect['probBrdU'])

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
                                             seg.get_forward_sequence(), 'T')

            str_gaps_in_t = ",".join(str(int(k)) for k in gaps_in_t)

            # function to convert UUID to int
            def uuid_first7_char_to_int(x):
                return int(x[0:7], 16)

            # set modification positions and values
            # ChEBI code for BrdU is 472552
            seg.set_tag("MM", f"T+{code}?,{str_gaps_in_t};", "Z")
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

    base_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    base_num = base_to_num[base]

    with ModBam(mod_bam_file) as bam:
        # get position and counts of bases in each interval
        t1 = list(map(lambda k:
                      bam.pileup(
                          k[0], k[1], k[2],
                          low_threshold=low_thres,
                          high_threshold=high_thres,
                          mod_base=mod_code),
                      intervals))

        # prepare to sum up the counts in each interval
        def f(index): return map(lambda k:
                                 int(sum(entry[index] for entry in k[1])), t1)

        # return values
        return f(base_num), f(base_num + 4), f(10), f(11)


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

    # function to convert UUID to int
    def uuid_first_7_char_to_int(x): return int(x[0:7], 16)

    # find relevant records and return data
    with ModBam(mod_bam_file) as bam:
        site_data = filter(
            lambda x: x[2] == read_id and start <= x[0] < end,
            ((
                r.rpos, r.qual,
                r.query_name
            ) for k in
                bam.reads(
                    contig, start, end, tag_name='XR',
                    tag_value=uuid_first_7_char_to_int(read_id)
                )
                for r in k.mod_sites))

        # return data
        # NOTE: in modBAM, a probability level x means the
        #   probability lies b/w x/256 and (x+1)/256 and x
        #   is a number b/w 0 and 255. So we return the midpoint
        #   of the interval given an x.
        return map(lambda x: (x[0], (x[1] + x[1] + 1) / (256 * 2)), site_data)
