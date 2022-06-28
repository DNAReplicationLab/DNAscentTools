import numpy as np
import itertools
import pysam
from modbampy import ModBam

def get_gaps_in_base_pos(posn, seq, base):
    """ Count occurences of base in seq in open interval (posn_i, posn_{i+1})

    Args:
        posn (list of int): Positions of interest, entries >= 0 & < len(seq)
        seq (str): Sequence of interest
        base (str): one-character string

    Returns:
        list of ints with length 1 shy of length of posn.
    """
    return [ seq[k[0] + 1 : k[1]].count(base) for k in 
        zip(posn[:-1], posn[1:]) ]

def convert_data_per_T_to_modBAM_fmt(seqData):
    """ Convert data per thymidine to modBAM format

    Using the SAM specification at
    https://samtools.github.io/hts-specs/SAMtags.pdf

    Args:
        seqData (list of floats): BrdU probabilities

    Returns:
        list of probs scaled to 0 to 255
    """

    # convert entries to integers b/w 0 and 255, both included
    return [ int(k) for k in 
                np.clip(np.floor([l * 256 for l in seqData]), 0, 255)]

def convert_detect_into_detect_stream(detectObj):
    """ Creates iterator that goes one detect record at a time. First item
    is a dictionary with keys comments and refFasta.

    Code adapted frm https://dnascent.readthedocs.io/en/latest/cookbook.html

    Args:
        detectObj (iterable): each line of detect file

    Yields:
        dict with keys readID, refContig, refStart, refEnd, strand,
          posOnRef, probBrdU, sixMerOnRef. Meanings are as in the detect file.
          NOTE: First ever item has keys comments (list of str) and
            refFasta (str).

    """

    currentEntry = {"comments": []}

    # split by '>' character, or per detection in other words
    for key, subiter in itertools.groupby(detectObj, 
        lambda x: x.startswith('>')):

        if key:

            # if header line, initialize record
            splitLine = next(subiter).rstrip().split()

            # store header information
            currentEntry = { 
                'readID': splitLine[0][1:],
                'refContig': splitLine[1],
                'refStart': int(splitLine[2]),
                'refEnd': int(splitLine[3]),
                'strand': splitLine[4],
                'posOnRef': [], 
                'probBrdU': [], 
                'sixMerOnRef': []
            }

        else:

            for line in subiter:

                # in dnascent detect files, a few comment lines
                # appear before the header. save these.
                # of particular interest is the line that starts
                # "#Genome" and contains the reference fasta file.
                if line.startswith('#'):

                    # store comments
                    currentEntry["comments"].append(line.rstrip())

                    # record fasta reference file
                    if line.startswith('#Genome'):
                        currentEntry["refFasta"] = line.replace("#Genome ","").rstrip()

                    continue

                # split data lines by whitespace
                splitLine = line.rstrip().split()

                # skip malformed lines if they exist
                if not len(splitLine) == 3:
                    continue

                # store data
                currentEntry['posOnRef'].append(int(splitLine[0]))
                currentEntry['probBrdU'].append(float(splitLine[1]))
                currentEntry['sixMerOnRef'].append(splitLine[2])

            # return current record
            yield currentEntry

def convert_dnascent_detect_to_modBAM_file(detectStream, filename,
    code = "472552", bamFile = True, shftRevStrndPos = 5, pgInfo = dict()):
    """ Convert dnascent detect data to modBAM format

    NOTE: About dnascent reverse offset issue and the param shftRevStrndPos.
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
        increment all posns on a reverse strand by 5, accomplished via the
        parameter shftRevStrndPos. If DNAscent output changes in future
        versions, set the parameter suitably.
    
    Args:
        detectStream (list of dicts): see below.
          * first entry in list must have a key comments = a list of comments
            from a dnascent detect file, and a key refFasta = location of
            fasta file. comments can be an empty list.
          * every other entry must have the keys readID (str), refStart (int), 
            refEnd (int), refContig (str), posOnRef (list of ints), 
            probBrdU (list of floats) which correspond to data per detect from
            the detect file.
        filename (str): name of output file
        code (str): default 472552, modification code (ChEBI code or 1 letter)
        bamFile (bool): default T, T/F = bam/sam file
        shftRevStrndPos (int): default 5. See note above, in fn docstr.
        pgInfo (dict): (default empty dict = no info). Info abt program, must
          have key 'ID'.

    Returns:
        None
    """

    # get comments and fasta file from detectfile
    # start making header to modBAM file
    initInfo = next(detectStream)
    header = { 'CO': [ 
                'comment: ' + k for k in initInfo['comments']
                ] 
            }

    refFastaFile = pysam.FastaFile(filename = initInfo['refFasta'])

    # load names of contigs, their order, and lengths from fasta file
    refNameLengthMap = dict(zip(refFastaFile.references, 
        refFastaFile.lengths))
    refNameIDMap = dict(zip(refFastaFile.references, 
        itertools.count()))

    # add contig and length information to the header
    header['SQ'] = [
        {'LN': v, 'SN': k} for k, v in refNameLengthMap.items()
    ]

    # add program information to the header
    if pgInfo:
        header['PG'] = [ pgInfo ]

    # iterate through each detect entry and add to modBAM/modSAM file
    with pysam.AlignmentFile(filename, "wb" if bamFile else "w", 
        header = header) as fb:

        for oneDetect in detectStream:

            # set up reverse strand apparatus
            fwdStrand = oneDetect['strand'] == 'fwd'
            fnRevIfNeeded = lambda x: x if fwdStrand else reversed(x)

            # reconstruct detect header but with underscores
            detectHeader = "_".join(
                str(oneDetect[k]) for k in 
                    ['readID','refContig','refStart','refEnd','strand']
                )

            # get sequence from reference
            seq = refFastaFile.fetch(oneDetect['refContig'],
                oneDetect['refStart'], 
                oneDetect['refEnd'])

            # create modBAM entry
            seg = pysam.AlignedSegment()
            seg.query_name = oneDetect['readID']
            seg.query_sequence = seq
            seg.flag = 0 if fwdStrand else 16
            seg.reference_id = refNameIDMap[oneDetect['refContig']]
            seg.reference_start = oneDetect['refStart']
            seg.mapping_quality = 255
            seg.cigar = [(0, len(seq))]
            seg.template_length = len(seq)
            seg.query_qualities = pysam.qualitystring_to_array("")

            # scale probabilities from 0 to 255
            probFrm0To255 = convert_data_per_T_to_modBAM_fmt(
                oneDetect['probBrdU'])

            # get number of skipped Ts b/w each T for which data is avlble
            # NOTE: the fn get_forward_sequence() returns seq if the reversed
            #   flag is not set, otherwise it returns the reversed complement.
            if fwdStrand:
                relPos = [k - oneDetect['refStart'] for k in 
                            oneDetect['posOnRef']]
            else:
                relPos = [oneDetect['refEnd'] - (k + shftRevStrndPos) 
                            for k in reversed(oneDetect['posOnRef'])]

            Tcoords = list(itertools.chain([-1], relPos))

            gapsInT = get_gaps_in_base_pos(Tcoords, 
                        seg.get_forward_sequence(), 'T')

            strGapsInT = ",".join(str(int(k)) for k in gapsInT)

            # function to convert UUID to int
            uuidFirst7CharToInt = lambda x: int(x[0:7], 16)

            # set modification positions and values
            # ChEBI code for BrdU is 472552
            seg.set_tag("MM", f"T+{code}?,{strGapsInT};" , "Z")
            seg.set_tag("ML", list(fnRevIfNeeded(probFrm0To255)))
            seg.set_tag("XR", uuidFirst7CharToInt(oneDetect['readID']), "i")
            seg.set_tag("XA", f"{detectHeader}" , "Z")

            # store modBAM entry
            fb.write(seg)

    refFastaFile.close()

def get_mod_counts_per_interval(modBamFile, intervals, base, 
    modCode, lowThres, highThres):
    """ Get modification counts given reference-coordinate intervals

    Args:
        modBamFile (str): path to modBAM file
        intervals (iter): each entry has 3 elements = contig,start,end
        base (str): 'A', 'G', 'C' or 'T'
        modCode (str): modification code
        lowThres (float): probability below which base marked as unmodified
        highThres (float): probability above which base marked as modified

    Returns:
        Tuple of 4 iterators: unmodCountInRevStrand, unmodCountInFwdStrand, 
            modCountInRevStrand, modCountInFwdStrand
    """

    baseToNum = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    baseNum = baseToNum[base]

    with ModBam(modBamFile) as bam:
        # get position and counts of bases in each interval
        t1 = list(map(lambda k: 
                         bam.pileup(
                             k[0], k[1], k[2],
                             low_threshold = lowThres, 
                             high_threshold = highThres,
                             mod_base = modCode), 
                     intervals))

        # prepare to sum up the counts in each interval
        f = lambda indx: map(lambda k: 
                                int(sum(entry[indx] for entry in k[1])), t1)
        
        # return values
        return f(baseNum), f(baseNum + 4), f(10), f(11)

def get_read_data_from_modBAM(modBamFile, readID, contig, start, end):
    """ Get modification probabilities corresponding to coords on a read

    Args:
        modBamFile (str): path to modBAM file
        readID (str): requested read ID
        contig (str): contig on reference genome
        start (num): start posn on ref genome, 0-based
        end (num): end posn on ref genome, 0-based

    Returns:
        Iterator w each entry = (ref posn, probability)
    """

    # function to convert UUID to int
    uuidFirst7CharToInt = lambda x: int(x[0:7], 16)

    # find relevant records and return data
    with ModBam(modBamFile) as bam:
        siteData = filter(
            lambda x: x[2] == readID,
            ((   
                l.rpos, l.qual,
                l.query_name
            ) for k in 
                bam.reads(
                    contig, start, end, tag_name = 'XR', 
                    tag_value = uuidFirst7CharToInt(readID)
                    ) 
            for l in k.mod_sites))

        # return data
        # NOTE: in modBAM, a probability level i means the 
        #   probability lies b/w i/256 and (i+1)/256 and i
        #   is a number b/w 0 and 255. So we return the midpoint
        #   of the interval given an i.
        return map(lambda x: (x[0], (x[1] + x[1] + 1)/(256*2)), siteData)