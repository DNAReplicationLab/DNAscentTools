import pandas as pd
import numpy as np
import pysam
from collections.abc import Callable, Iterable
from collections import namedtuple
from functools import reduce
from itertools import count, accumulate, pairwise

ModBase = namedtuple('ModBase', (
    'read_id', 'ref_pos', 'fwd_seq_pos', 'ref_strand', 'mod_strand',
    'can_base', 'mod_base', 'mod_qual'), defaults=('', -1, -1, 'unstranded', 'unstranded', '', '', -1))


# above tuple based on the class ModInfo of the now un-maintained modbampy package


class HeaderLines:
    """ Class to store and manipulate header lines of BAM/modBAM file """

    def __init__(self, program_name: str, program_version: str, program_call: str = ""):
        """ Initialize HeaderLines object

        Args:
            program_name: name of program
            program_version: version of program
            program_call: command line call of program, default is empty string

        Returns:
            HeaderLines object
        """
        self.header_lines = []  # header lines
        self.is_input_stopped = False  # whether input is stopped
        self.last_pg_line_position = -1  # position of the last PG line in header_lines (see note below)
        self.pg_data = []  # data from PG lines
        self.program_name = program_name  # name of current program that is manipulating the modBAM file
        self.program_version = program_version  # version of current program that is manipulating the modBAM file
        self.program_call = program_call  # command line call of current program that is manipulating the modBAM file

        """
        Note:
        I'm going to explain some parts of the BAM header @PG lines using an example.
        For more information, see https://samtools.github.io/hts-specs/SAMv1.pdf or a later version.
        
        Let's say we have a bam file that has so far been operated on by four programs in the order
        program1, program2, program3, and program1 again.
        BAM files represent this information as:
        - program1.1 <- program3
        - program3 <- program2
        - program2 <- program1
        - program1 <- N/A
        Please note that as each line contains information about itself and the preceding step, the line order
        can be changed without changing the information.
        
        The header lines of the bam file that contain the information above look like this
        (the big spaces are tabs; the order of the lines and the tags within the lines is not important):

        @PG     ID:program1     PN:program1     VN:1.0
        @PG     ID:program2     PN:program2     VN:1.0     PP:program1
        @PG     ID:program3     PN:program3     VN:1.0     PP:program2
        @PG     ID:program1.1   PN:program1     VN:1.0     PP:program3

        PN means program name and VN means program version. Programs can have any name and have any version.
        In the example above, we've chosen program names same as the ids and program versions of 1.0
        for all programs for simplicity.
        
        So, 
        - we've to locate the last run program, which is program1.1. This may not be the last line.
        - prepare an @PG line with the details of the current program using the ID of the last program as PP.
        - insert the new @PG line. Now, this line can be inserted anywhere in the header, but to maintain
        readability, I've chosen to insert it after the last @PG line (remember that this may not necessarily
        be the latest program to be run on the bam file). This position is called last_pg_line_position in our class. 
        """

    def add_line(self, header_line: str) -> None:
        """ Add header line to list of header lines """

        if self.is_input_stopped:
            raise ValueError("Cannot add header lines!")

        self.header_lines.append(header_line)

        # put id and pp data from PG lines in pg_data if available.
        if header_line.startswith("@PG"):
            self.last_pg_line_position = len(self.header_lines) - 1
            pg_id = None
            pg_pp = "*"
            for tag in header_line.split("\t"):
                if tag.startswith("ID:"):
                    pg_id = tag[3:]
                elif tag.startswith("PP:"):
                    pg_pp = tag[3:]
            if pg_id:
                self.pg_data.append((pg_id, pg_pp))

    def print_lines(self) -> None:
        """ Print header lines to stdout and stop input """
        for header_line in self.header_lines:
            print(header_line)
        self.is_input_stopped = True

    def stop_input(self) -> None:
        """ Stop input """
        self.is_input_stopped = True

    def has_print_action_completed(self) -> bool:
        """ Return whether print action has been completed """
        return self.is_input_stopped

    def reorder_pg_data(self) -> None:
        """ Reorder pg_data so that it starts with the first program and ends with the last program """

        # traverse the list to find the first program
        first_program = next((pg_id for pg_id, pg_pp in self.pg_data if pg_pp == "*"), None)

        # reorder pg_data from the first to the last program
        if first_program:
            new_pg_data = [(first_program, "*")]
            while True:
                for pg_id, pg_pp in self.pg_data:
                    if pg_pp == first_program:
                        new_pg_data.append((pg_id, pg_pp))
                        first_program = pg_id
                        break
                else:
                    break
            self.pg_data = new_pg_data

    def add_current_program_to_header(self) -> None:
        """ Add current program details to header lines """

        # arrange pg data so that it starts with the first program and ends with the last program
        self.reorder_pg_data()

        # prepare new @PG line
        header_line_entries = ["@PG", f"PN:{self.program_name}", f"VN:{self.program_version}"]

        # count number of occurrences of the current program
        n_program_calls = self.get_number_of_instances_of_current_program()

        # if program has been called before, append .n_program_calls to the program name
        if n_program_calls > 0:
            header_line_entries.insert(1, f"ID:{self.program_name}.{n_program_calls}")
        else:
            header_line_entries.insert(1, f"ID:{self.program_name}")

        # add get_modBAM_with_nascent_reads to header lines after last PG line
        if not self.pg_data == []:
            header_line_entries.insert(3, f"PP:{self.pg_data[-1][0]}")

        # if command call is available, append it to header line entries
        if self.program_call:
            header_line_entries.append(f"CL:{self.program_call}")

        self.header_lines.insert(self.last_pg_line_position + 1,
                                 "\t".join(header_line_entries))

    def get_number_of_instances_of_current_program(self) -> int:
        """ Return number of times the current program has been called so far"""
        return len([pg_id for pg_id, pg_pp in self.pg_data if self.program_name in pg_id])


class ModBamRecordProcessor:
    """Process data from data in one mod bam line
    To really understand what probability_modbam_format and thymidine_gaps mean, refer to
    https://samtools.github.io/hts-specs/SAMtags.pdf
    To understand what other BAM fields mean, refer to https://samtools.github.io/hts-specs/SAMv1.pdf

    Attributes:
        threshold (float): threshold above (below) which thymidine is regarded as modified (unmodified)
        code (str): code of thymidine modification
        use_xa_tag (bool): use contents of XA tag to set read_id instead of the usual first entry
        allow_non_na_mode (bool): can interpret missing bases as "unmodified" rather than "missing" depending on tag
        base (str): unmodified base to be considered
        force_missing (bool): force missing bases to be interpreted as "missing" rather than
                "unmodified", irrespective of dot/blank/question mark tag used in mod bam
        read_id (str): read id of the record
        is_alt_read_id (bool): if True, read id field stores something else e.g. from the XA tag.
        flag (int): flag of the record
        contig (str): contig of the record
        start (int): start position of the record, will be stored as 0-based
        end (int): end position of the record, will be stored as 0-based
        qual (int): quality of the record
        cigar (str): cigar string of the record
        seq (str): sequence of the record
        is_rev (bool): whether the record is reverse
        is_unmapped (bool): whether the record is unmapped
        fwd_seq (str): forward sequence of the record
        ref_to_query_tbl (list): list of tuples of reference and query positions calculated from cigar string
        raw_probability_modbam_format (list): entries are 0 to 255
        probability_modbam_format (list): entries are 0 or 1 meaning unmodified or modified
        thymidine_gaps (list): entries are ints, show number of skipped thymidines.
        fwd_seq_thymidine_coordinates (list): coordinates of thymidines in forward sequence
        fwd_seq_reference_coordinates (list): coordinates of thymidines in forward sequence mapped to reference
    """
    threshold: float
    code: str
    use_xa_tag: bool
    allow_non_na_mode: bool
    base: str
    force_missing: bool

    read_id: str
    is_alt_read_id: bool
    flag: int
    contig: str
    start: int
    end: int
    qual: int
    cigar: str
    seq: str
    is_rev: bool
    is_unmapped: bool

    fwd_seq: str
    ref_to_query_tbl: list[tuple[int, int]]

    raw_probability_modbam_format: list[int]
    probability_modbam_format: list[int]
    thymidine_gaps: list[int]
    fwd_seq_thymidine_coordinates: list[int]
    fwd_seq_reference_coordinates: list[int]

    def __init__(self, threshold, code, use_xa_tag=False, allow_non_na_mode=False, base="T", force_missing=False):
        """Initializes the instance

        Args:
            threshold (float): defines threshold of this instance
            code (str): defines code of this instance
            use_xa_tag (bool): (default False) defines use_xa_tag of this instance
            allow_non_na_mode (bool): (default False) defines allow_non_na_mode of this instance
            base (str): (default "T") defines base of this instance
            force_missing (bool): (default False) defines force_missing of this instance
        """

        self.threshold = threshold
        self.code = code
        self.use_xa_tag = use_xa_tag
        self.allow_non_na_mode = allow_non_na_mode
        self.base = base
        self.force_missing = force_missing

        self.delete_data()

    def delete_data(self) -> None:
        """ Remove stored data (but not parameters).

        Returns:
            None
        """
        self.read_id = ""
        self.is_alt_read_id = False
        self.flag = 0
        self.contig = ""
        self.start = 0
        self.end = 0
        self.qual = 0
        self.cigar = ""
        self.seq = ""
        self.is_rev = False
        self.is_unmapped = False

        self.fwd_seq = ""
        self.ref_to_query_tbl = []

        self.raw_probability_modbam_format = []
        self.probability_modbam_format = []
        self.thymidine_gaps = []
        self.fwd_seq_thymidine_coordinates = []
        self.fwd_seq_reference_coordinates = []

    def process_modbam_line(self, x) -> None:
        """Process data from one mod bam line"""

        # delete previously stored data
        self.delete_data()

        # variables to receive data from mod bam line
        ml_part = ""
        mm_part = ""

        for cnt, k in zip(count(), x.strip().split("\t")):

            if cnt >= 11 and k.startswith(("MM:Z:", "Mm:Z:")):
                mm_part = k
            elif cnt >= 11 and k.startswith(("Ml:B:C,", "ML:B:C,")):
                ml_part = k
            elif cnt == 0:
                self.read_id = k
            elif (not (cnt == 0 or self.is_unmapped)) and k.startswith("XA:Z:") and self.use_xa_tag:
                self.read_id = k[5:]
                self.is_alt_read_id = True
            elif 1 <= cnt <= 9:
                if cnt == 1:
                    self.flag = int(k)
                    self.is_unmapped = (self.flag & 4 == 4)
                    self.is_rev = (self.flag & 16 == 16)
                elif cnt == 2:
                    self.contig = k
                    if self.contig == "*":
                        self.is_unmapped = True
                elif cnt == 3:
                    self.start = int(k) - 1
                    if self.start == -1:
                        self.is_unmapped = True
                elif cnt == 4:
                    self.qual = int(k)
                elif cnt == 5:
                    self.cigar = k
                    if self.cigar == "*":
                        self.is_unmapped = True
                    if not self.is_unmapped:
                        self.ref_to_query_tbl = cigar_to_ref_to_query_tbl(self.cigar, self.start)
                        self.end = self.ref_to_query_tbl[-1][0]
                        assert self.end >= self.start
                elif cnt == 9:
                    self.seq = k
                    if self.seq == "*":
                        raise ValueError("We cannot deal with * sequences!")

        # get the forward sequence
        if self.is_rev:
            self.fwd_seq = reverse_complement(self.seq)
        else:
            self.fwd_seq = self.seq

        mod_data = list(filter(lambda y: y["base"] == self.base and y["mod_code"] == self.code,
                               parse_modBAM_modification_information(f"{mm_part}\t{ml_part}")))

        if len(mod_data) > 1:
            raise ValueError("It appears that there are many pieces corresponding to the same modification in a read!")

        elif len(mod_data) == 1:

            self.raw_probability_modbam_format = mod_data[0]["prob"]
            self.probability_modbam_format = list(map(lambda y: 1 if y >= self.threshold * 256 else 0,
                                                      self.raw_probability_modbam_format))
            self.thymidine_gaps = mod_data[0]["pos"]

            # insert zeroes for missing bases if non na mode is allowed
            if self.allow_non_na_mode and mod_data[0]['mode'] == '.' and not self.force_missing:
                self.insert_zeroes_according_to_gaps()

            # convert gap coordinates to normal coordinates
            self.fwd_seq_thymidine_coordinates = convert_gap_coordinates_to_normal_coordinates(self.fwd_seq,
                                                                                               self.thymidine_gaps,
                                                                                               self.base)
            # check that the positions all have the same base as the base of interest
            assert all(self.fwd_seq[k] == self.base for k in self.fwd_seq_thymidine_coordinates)

            # get reference coordinates if not unmapped
            if not self.is_unmapped:
                self.fwd_seq_reference_coordinates = convert_forward_seq_coordinates_to_ref_coordinates(
                    self.fwd_seq_thymidine_coordinates, self.ref_to_query_tbl,
                    self.is_rev)

    def has_data(self) -> bool:
        """Check if data is available"""

        return (len(self.read_id) > 0 and
                len(self.probability_modbam_format) > 0 and
                len(self.thymidine_gaps) > 0)

    def has_read(self) -> bool:
        """Check if a read has been read by the processor"""

        return len(self.read_id) > 0

    def count_bases(self, mask_ref_interval: tuple[int, int] = (-1, -1),
                    mask_fwd_seq_interval: tuple[int, int] = (-1, -1)) -> tuple[int, int, int]:
        """Count number of modified, unmodified, and bases where modification information is missing.

        Args:
            mask_ref_interval: (default (-1,-1)) only return counts for reference positions such that a_0 <= a < a_1,
                               where a_0 and a_1 are the two elements of mask_ref_interval. Default unused.
                               If you set this and the read is unmapped, you'll get an error.
            mask_fwd_seq_interval: (default (-1,-1)) only return counts for forward-sequence positions such that
                                   a_0 <= a < a_1, where a_0 and a_1 are the two elements of mask_fwd_seq_interval.
                                   Default unused. This works even if the read is unmapped.

        Returns:
            tuple of three ints in the order above.
        """

        if not self.has_data():
            if not self.has_read():
                raise ValueError("No data available!")
            else:
                return 0, 0, 0

        if mask_ref_interval != (-1, -1):

            if mask_fwd_seq_interval != (-1, -1):
                raise ValueError("Cannot mask both reference and forward-sequence intervals!")
            if self.is_unmapped:
                raise ValueError("Cannot mask reference interval for unmapped reads!")
            if mask_ref_interval[0] > mask_ref_interval[1]:
                raise ValueError("Bad mask_ref_interval!")

            indices = [k for k, a in enumerate(self.fwd_seq_reference_coordinates)
                       if mask_ref_interval[0] <= a < mask_ref_interval[1]]

        elif mask_fwd_seq_interval != (-1, -1):

            if mask_fwd_seq_interval[0] > mask_fwd_seq_interval[1]:
                raise ValueError("Bad mask_fwd_seq_interval!")

            indices = [k for k, a in enumerate(self.fwd_seq_thymidine_coordinates)
                       if mask_fwd_seq_interval[0] <= a < mask_fwd_seq_interval[1]]

        else:

            indices = range(len(self.probability_modbam_format))

        count_mod_plus_unmod = len(indices)
        count_mod = sum(self.probability_modbam_format[k] for k in indices)
        count_missing = sum(self.thymidine_gaps[k] for k in indices)

        return count_mod, count_mod_plus_unmod - count_mod, count_missing

    def insert_zeroes_according_to_gaps(self) -> None:
        """ Interpret missing thymidines as unmodified and adjust our data.

        Returns:
            None
        """
        if not self.has_data():
            raise ValueError("No data available!")

        # first, we account for implicit skips
        num_end_insertions = (len([k for k in self.fwd_seq if k == self.base]) - sum(self.thymidine_gaps)
                              - len(self.thymidine_gaps))

        # find positions where we need to insert zeroes and insert them.
        thymidine_zero_insertions = reduce(lambda x, y: x + [y[0]] * y[1],
                                           filter(lambda x: x[1] > 0,
                                                  enumerate(self.thymidine_gaps + [num_end_insertions])),
                                           [])
        self.probability_modbam_format = list(np.insert(self.probability_modbam_format, thymidine_zero_insertions, 0))
        self.raw_probability_modbam_format = list(np.insert(self.raw_probability_modbam_format,
                                                            thymidine_zero_insertions, 0))

        # after adjustment for unmodified bases, all thymidine gaps are set to zero
        self.thymidine_gaps = [0 for _ in range(len(self.probability_modbam_format))]

    def detect_header(self, force_unmapped: bool = False) -> str:
        """ Returns a header as would have been prepared by DNAscent detect

        Args:
            force_unmapped: (default False) if True, force unmapped header i.e. use query parameters and not
                            reference parameters, so contig is set to "unmapped" and strand is set to "unstranded",
                            and coordinates are set to 0 and read length - 1.

        Returns:
            detect header
        """
        if self.is_unmapped or force_unmapped:
            if self.is_alt_read_id:
                read_id, contig, start, end, orn = process_detect_index(self.read_id)
            else:
                read_id = self.read_id
            return f">{read_id} unmapped 0 {len(self.seq)} unstranded"
        elif self.use_xa_tag and self.is_alt_read_id:
            return f">{self.read_id.replace('_', ' ')}"
        else:
            return f">{self.read_id} {self.contig} {self.start} {self.end} " \
                   f"{'rev' if self.is_rev else 'fwd'}"

    def mod_data_to_table(self, move_parallel_top_ref_strand: bool = False) -> Iterable[ModBase]:
        """ Return modification data in a tabular format

        Args:
            move_parallel_top_ref_strand: (default False) if True, move parallel to top strand of reference genome.
                                          Default is to move parallel to the forward sequence, which is the same
                                          direction as the reference genome for fwd mapping and opposite for rev
                                          mapping. Option is ignored if the read is unmapped.

        Returns:
            iterator with dictionaries with keys read_id, fwd_seq_pos, ref_pos, mod_qual.
            We've chosen similar headers as `modkit extract` so that we can compare the two outputs.

        """

        if not self.has_data():
            return iter([])

        rev = move_parallel_top_ref_strand and (not self.is_unmapped) and self.is_rev
        fwd_seq_coords = reversed(self.fwd_seq_thymidine_coordinates) if rev else self.fwd_seq_thymidine_coordinates
        ref_coords = reversed(self.fwd_seq_reference_coordinates) if rev else self.fwd_seq_reference_coordinates
        probs = convert_probabilities_from_modBAM_to_normal(
            reversed(self.raw_probability_modbam_format) if rev else self.raw_probability_modbam_format)

        return (ModBase(read_id=self.read_id, fwd_seq_pos=b, ref_pos=c, mod_qual=d, can_base=self.base,
                        mod_base=self.code,
                        ref_strand='unmapped' if self.is_unmapped else ('-' if self.is_rev else '+'),
                        mod_strand='+')
                for b, c, d in zip(fwd_seq_coords, ref_coords, probs))

    def calculate_autocorrelations(self, chunk_size: int = 10000, window_size: int = 300,
                                   chunk_overlap_fraction: float = 0.0) -> tuple[list[float]]:
        """ Calculate autocorrelations in the thresholded and windowed modification probability per given chunk size.

        Args:
            chunk_size: (default 10000) split read into (non-overlapping by default) chunks of this many thymidines and
                        calculate correlations separately for each chunk. The last chunk will be smaller if the read
                        length is not a multiple of chunk_size (and may be rejected if it is too small).
                        Also see chunk_overlap_fraction.
            window_size: (default 300) window data in this window size (of thymidines, non-overlapping windows)
                         before calculating autocorrelations.
            chunk_overlap_fraction: (default 0.0) When we split read into chunks, stipulate that chunks should overlap
                                    by this fraction of chunk_size. This is a float >= 0 and < 1. Do not choose values
                                    close to 1 as this means a high level of chunk overlap.

        Returns:
            tuple of N lists of floats where N = read length // chunk_size.
            Each float list contains autocorrelations for the corresponding chunk at different lag distances.

        """

        # function to chunk and window data
        def window_non_overlapping(x, l_win, l_min):
            if len(x) < max(l_win, l_min):
                return
            v = np.lib.stride_tricks.sliding_window_view(x, l_win)[::l_win, :]
            return v.mean(axis=-1)

        # set stride length for chunking
        stride = int(chunk_size * (1 - chunk_overlap_fraction))

        # check that stride is between 1 and chunk_size
        if not 1 <= stride <= chunk_size:
            raise ValueError("Bad chunk_overlap_fraction!")

        # divide self.probability_modbam_format into chunks of size chunk_size and window each chunk
        chunked_data = [window_non_overlapping(np.array(self.probability_modbam_format[i:i + chunk_size], dtype=float),
                                               window_size, chunk_size / 2)
                        for i in range(0, len(self.probability_modbam_format), stride)]

        # normalize each chunk by its mean and the sd
        norm_chunked_data = [(chunk - np.mean(chunk)) / (np.std(chunk) + 1e-10) for chunk in chunked_data
                             if chunk is not None]

        # calculate autocorrelations for each chunk
        return tuple((np.correlate(chunk, chunk, mode='full')[(len(chunk) - 1):(2 * len(chunk) - 1)] *
                      1 / np.linspace(len(chunk), 1, num=len(chunk))).tolist()
                     for chunk in norm_chunked_data)


def parse_modBAM_modification_information(input_string: str) -> list[dict]:
    r""" Parse modification information from modBAM string

    Args:
        input_string: modBAM string containing modification information e.g. MM:Z:T+T?,1,2,3;\tML:B:C,1,2,3

    Returns:
        list of dictionaries containing the keys base, mod_strand, mod_code, mode, pos, prob

    Examples:
        >>> parse_modBAM_modification_information("MM:Z:T+T?,1,2,3;\tML:B:C,1,2,3")
        [{'base': 'T', 'mod_strand': '+', 'mod_code': 'T', 'mode': '?', 'pos': [1, 2, 3], 'prob': [1, 2, 3]}]
        >>> parse_modBAM_modification_information("MM:Z:G-h,1,2;\tML:B:C,100,120")
        [{'base': 'G', 'mod_strand': '-', 'mod_code': 'h', 'mode': '.', 'pos': [1, 2], 'prob': [100, 120]}]
    """

    def int_within_num(x: str, num: int) -> int:
        """Converts a string to an int if is a number between 0 and num (both inclusive), else raise Exception"""
        try:
            n = int(x)
            assert (0 <= n <= num)
            return n
        except (ValueError, AssertionError):
            raise ValueError("Bad probability!")

    # Initialize the output list
    result = []

    # get the MM and ML parts which correspond to modification position and probability
    temp_1, ml_part = input_string.split("\t")
    temp_2 = temp_1.split(";")
    mm_parts = temp_2[:-1]

    # ensure that ML and MM strings appear where they should and remove after check
    assert (mm_parts[0].startswith(("MM:Z:", "Mm:Z:")))
    assert (ml_part.startswith(("Ml:B:C,", "ML:B:C,")))
    mm_parts[0] = mm_parts[0][5:]
    ml_part = ml_part[7:]

    # make a list for probability data
    prob = list(map(lambda x: int_within_num(x, 255), ml_part.split(",")))

    # counter to slice probability data
    cnt_prob = 0

    for mm_part in mm_parts:

        # split the MM part into its components: first bit is about the modification code, the rest are positions
        l = mm_part.split(",")
        assert (len(l) > 1)

        # set modification code. This could be a single letter, a number, or a string of letters
        is_ending_with_mode = l[0].endswith(("?", "."))
        mod_code_str = l[0][2:-1] if is_ending_with_mode else l[0][2:]
        assert (len(mod_code_str) > 0)

        if mod_code_str.isdigit():
            mod_code_list = [int(mod_code_str)]
        elif mod_code_str.isalpha():
            mod_code_list = list(mod_code_str)
            # check that the modification codes are unique
            assert (len(mod_code_list) == len(set(mod_code_list)))
        else:
            raise ValueError("Bad modification code!")

        # store modification information per modification code
        large_num = 1_000_000_000
        for mod_code in mod_code_list:
            base = l[0][0]
            mod_strand = l[0][1]
            assert (base in ["A", "G", "C", "T", "N"])
            assert (mod_strand in ["+", "-"])
            result.append({"base": base,
                           "mod_strand": mod_strand,
                           "mod_code": mod_code,
                           "mode": l[0][-1] if is_ending_with_mode else ".",
                           "pos": list(map(lambda x: int_within_num(x, large_num), l[1:])),  # large_num is arbitrary
                           "prob": prob[cnt_prob: cnt_prob + len(l) - 1]
                           })
        cnt_prob += len(l) - 1

    return result


def reverse_complement(x: str) -> str:
    """ Get reverse complement of a DNA sequence using samtools reversed read apparatus

    Args:
        x: DNA sequence

    Returns:
        Reverse complement of DNA sequence

    Examples:
        >>> reverse_complement("AGTATGGCT")
        'AGCCATACT'
    """

    seg = pysam.AlignedSegment()
    seg.flag = 16
    seg.query_sequence = x
    return seg.get_forward_sequence()


def convert_forward_seq_coordinates_to_ref_coordinates(fwd_seq_coordinates: list[int],
                                                       ref_to_query_tbl: list[tuple[int, int]],
                                                       is_rev: bool = False):
    """ Convert forward sequence coordinates to reference coordinates

    Args:
        fwd_seq_coordinates: a list of coordinates along the forward sequence
            (forward sequence is the same as the query sequence for a forward read and is the reverse complement
            of the query sequence for a reverse read). Some coordinates can be missing.
        ref_to_query_tbl: a list where each entry is a tuple of two integers: reference and query coordinates.
            This can be calculated using the CIGAR string. Not all reference and query coordinates need to be
            present as we can use linear interpolation.
            e.g. [(100, 0), (110, 10), (120, 10), (120, 20), (130, 30)] means the first ten bases of the query map
            to the first ten bases of the reference starting at position 100. The next ten bases of the reference
            do not map to the query, and the next ten bases of the query do not map to the reference, and the
            last ten bases of the query map to the last ten bases of the reference starting at position 120.
        is_rev: (default False) if True, forward sequence is the reverse complement of query sequence,
                else it is the same

    Returns:
        list of reference coordinates corresponding to the forward sequence coordinates, with -1 for coordinates
        on the forward sequence that do not map to the reference sequence.

    Examples:
        >>> fwd_seq_T_coordinates = [0, 9, 10, 15, 20, 29, 30]
        >>> ref_to_query = [(100, 0), (110, 10), (120, 10), (120, 20), (130, 30)]
        >>> print(convert_forward_seq_coordinates_to_ref_coordinates(fwd_seq_T_coordinates, ref_to_query))
        [100, 109, -1, -1, 120, 129, -1]
        >>> print(convert_forward_seq_coordinates_to_ref_coordinates([30, 9, 10, 29, 20, 15, 0], ref_to_query)) # shuf
        [-1, 109, -1, 129, 120, -1, 100]
        >>> print(convert_forward_seq_coordinates_to_ref_coordinates([0, 9, 10, 15, 20, 29, 30], ref_to_query, True))
        [129, 120, -1, -1, 109, 100, -1]
        >>> print(convert_forward_seq_coordinates_to_ref_coordinates([40, -1], ref_to_query))
        [-1, -1]
        >>> print(convert_forward_seq_coordinates_to_ref_coordinates([40, 104, 204, 404], [(200, 0), (400, 200)]))
        [240, 304, -1, -1]
        >>> print(convert_forward_seq_coordinates_to_ref_coordinates([0, 20, 40, 50, 80],
        ...  [(200, 0), (200, 33), (429, 262)]))
        [-1, -1, 207, 217, 247]
    """
    # Sort the ref_to_query_tbl based on query coordinates for easy interpolation
    ref_to_query_tbl = sorted(ref_to_query_tbl, key=lambda x: x[1])

    # get the starting reference coordinate
    ref_start = ref_to_query_tbl[0][0]

    # get the query size
    query_size = ref_to_query_tbl[-1][1]

    # perform some simple checks first
    assert ref_to_query_tbl[0][1] == 0
    assert ref_start >= 0
    assert query_size > 0

    # Examine the N-1 intervals between the N query coordinates in ref_to_query_tbl and only retain
    # those where reference and query change by equal amounts.
    def retain_matched_intervals(x):
        a = x[0]
        b = x[1]
        if a[0] < b[0] and a[1] < b[1] and a[0] - b[0] == a[1] - b[1]:
            return a[1], b[1], a[0]
        elif not ((a[0] == b[0] and a[1] < b[1]) or (a[1] == b[1] and a[0] < b[0])):
            raise ValueError("Bad ref_to_query_tbl!")
        return None

    intervals = [z for z in map(retain_matched_intervals, pairwise(ref_to_query_tbl)) if z is not None]

    # sort forward sequence coordinates remembering the original order
    fwd_seq_coordinates_sorted = sorted(enumerate([query_size - 1 - k if is_rev else k for k in fwd_seq_coordinates]),
                                        key=lambda x: x[1])

    # set up an output list
    output = [(-1, -1) for _ in range(len(fwd_seq_coordinates))]

    # set up two pointers, one for the intervals and one for the forward sequence coordinates
    i = 0
    j = 0

    # iterate over the intervals and obtain desired mapping
    while j < len(fwd_seq_coordinates_sorted):
        if i >= len(intervals):
            output[j] = (fwd_seq_coordinates_sorted[j][0], -1)
            j += 1
        elif intervals[i][0] > fwd_seq_coordinates_sorted[j][1]:
            output[j] = (fwd_seq_coordinates_sorted[j][0], -1)
            j += 1
        elif intervals[i][0] <= fwd_seq_coordinates_sorted[j][1] < intervals[i][1]:
            output[j] = (fwd_seq_coordinates_sorted[j][0],
                         intervals[i][2] + fwd_seq_coordinates_sorted[j][1] - intervals[i][0])
            j += 1
        else:
            i += 1

    return [m[1] for m in sorted(output, key=lambda x: x[0])]


def get_mod_mean_per_read_per_interval(mod_bam_file: str, contig: str, start: int, end: int, thres: float):
    """ Get mean modification probability per read confined to given interval

    Args:
        mod_bam_file: path to modBAM file
        contig: contig on reference genome
        start: start position on ref genome, 0-based
        end: end position on ref genome, 0-based
        thres: value above (below) which T is called as BrdU (T)

    Returns:
        Iterator w each entry = (read id, mean BrdU, num bases)

    Examples:
        >>> output = list(get_mod_mean_per_read_per_interval("sample.bam", "dummyI", 9, 13, 0.01))
        >>> len(output)
        1
        >>> output[0]
        ('5d10eb9a-aae1-4db8-8ec6-7ebb34d32575', 1.0, 2)
        >>> output = list(get_mod_mean_per_read_per_interval("sample.bam", "dummyII", 6, 30, 0.5))
        >>> len(output)
        1
        >>> output[0]
        ('fffffff1-10d2-49cb-8ca3-e8d48979001b', 0.2, 5)
    """

    df = pd.DataFrame(
        map(lambda x: (x[0], 1) if x[2] >= thres else (x[0], 0),
            get_raw_data_from_modBAM(mod_bam_file, contig, start, end)),
        columns=['read_id', 'brdU']
    )

    return ((k[0], k[1][0], k[1][1])
            for k in df.groupby('read_id').agg(lambda x: (x.mean(), len(x))).itertuples(name=None))


def get_raw_data_from_modBAM(mod_bam_file: str, contig: str, start: int, end: int,
                             base: str = 'T', code: str = 'T',
                             read_id: str = '',
                             allow_multiple_entries_input_read_id: bool = False) -> Iterable[tuple[str, int, float]]:
    """ Gets modification probabilities from coords in interval on all reads 

    Args:
        mod_bam_file: path to modBAM file
        contig: contig on reference genome
        start: start position on ref genome, 0-based
        end: end position on ref genome, 0-based
        base: base to be considered for modification, default is "T"
        code: code of modification, default is "T"
        read_id: (default "") read id to be retained, if not provided, all reads are considered
        allow_multiple_entries_input_read_id: (default False) if True, allow multiple entries for the given input
                                              read id e.g. multiple alignments for a read id that map to the same
                                              region on the reference genome.

    Returns:
        Iterator w each entry = (read id, ref position, probability)

    Examples:
        >>> output = list(get_raw_data_from_modBAM("sample.bam", "dummyI", 9, 13))
        >>> len(output)
        2
        >>> output[0]
        ('5d10eb9a-aae1-4db8-8ec6-7ebb34d32575', 9, 0.017578125)
        >>> output[1]
        ('5d10eb9a-aae1-4db8-8ec6-7ebb34d32575', 12, 0.029296875)
    """

    # find relevant records and return data
    if read_id:
        alignments = pysam.view("-e", f"qname==\"{read_id}\"", mod_bam_file, f"{contig}:{start}-{end}")
    else:
        alignments = pysam.view(mod_bam_file, f"{contig}:{start}-{end}")

    # threshold below is unused in this script
    mod_bam_parser = ModBamRecordProcessor(0.01, code, allow_non_na_mode=True, base=base)

    if len(alignments) > 0:

        split_lines = alignments.splitlines()
        if sum([1 for k in split_lines if k]) > 1 and read_id and not allow_multiple_entries_input_read_id:
            raise ValueError(f"More than one record found for read_id {read_id}!")

        for x in split_lines:
            if x:
                mod_bam_parser.process_modbam_line(x)
                if mod_bam_parser.has_data():
                    for k in filter(lambda y: start <= y.ref_pos < end,
                                    mod_bam_parser.mod_data_to_table(move_parallel_top_ref_strand=True)):
                        yield k.read_id, k.ref_pos, k.mod_qual
    else:
        return iter([])


def convert_gap_coordinates_to_normal_coordinates(seq: str, gaps: list[int], base: str = 'T') -> list[int]:
    """ Convert coordinates from gap notation to standard notation

    Args:
        seq: reference sequence
        gaps: each entry a number >= 0 indicating number of bases to skip
        base: base to be considered for modification, default is "T"

    Returns:
        list of coordinates in standard notation

    Examples:
        >>> convert_gap_coordinates_to_normal_coordinates("AGTATGGCT", [0, 0, 0])
        [2, 4, 8]
        >>> convert_gap_coordinates_to_normal_coordinates("AGTATGGCT", [0, 1])
        [2, 8]
        >>> convert_gap_coordinates_to_normal_coordinates("AGTATGGCTA", [0, 1], "A")
        [0, 9]
        >>> convert_gap_coordinates_to_normal_coordinates("AGTATGGCT", [0, 1], "G")
        [1, 6]
    """

    candidate_coordinates = [k for k in range(len(seq)) if seq[k] == base]
    bases_to_retain = [m - 1 for m in accumulate(k + 1 for k in gaps)]

    return list(map(lambda x: candidate_coordinates[x], bases_to_retain))


def convert_probabilities_from_modBAM_to_normal(prob_0_to_255: list[int]) -> list[float]:
    """ Convert modification probabilities from modBAM notation to normal notation.
    In modBAM notation, an integer x means the modification probability lies between
    x/256 and (x+1)/256 (more precisely, x/256 <= prob < (x+1)/256). So, when we convert
    the probability to a normal notation, we report it as the midpoint of the interval.

    Args:
        prob_0_to_255: each entry a number from 0 to 255

    Returns:
        list of probabilities in normal notation

    Examples:
        >>> convert_probabilities_from_modBAM_to_normal([128, 129])
        [0.501953125, 0.505859375]
    """
    return list(map(lambda x: (x + x + 1) / (256 * 2), prob_0_to_255))


def modBAM_record_to_detect(modbam_line: str, fasta_file_name: str = "", reverse_shift: int = 5, code: str = 'T') \
        -> str:
    r""" Convert raw data from one read of modbam to detect format

    Args:
        modbam_line: one line of data from modBAM file. must contain MM, ML and XA tags
        fasta_file_name: (default "") fasta file containing reference genome. If not provided, will operate
                         in unmapped mode, where we report modification data along
                         the fwd-seq read coordinate and not the reference coordinate and mention unmapped and
                         unstranded for the contig and the strand respectively in the header.
        reverse_shift: (default 5) on reversed reads, detect format reports shifted coordinates. undo if necessary.
        code: (default 'T') thymidine modification code.

    Returns:
        detect record for given modBAM line

    Examples:
        >>> sample_first_bam_fields = "loremipsum\t0\tdummyI\t2\t0\t12M\t*\t0\t12\tGCTAGCTATCGT\t*"
        >>> sample_line = f"{sample_first_bam_fields}\tMM:Z:T+T?,1;\tML:B:C,10\tXA:Z:loremipsum_dummyI_1_13_fwd"
        >>> modBAM_record_to_detect(sample_line, "sample.fa")
        '>loremipsum dummyI 1 13 fwd\n7\t0.041015\tTATCGT'
        >>> modBAM_record_to_detect(sample_line)
        '>loremipsum unmapped 0 12 unstranded\n6\t0.041015\tTATCGT'
        >>> sample_second_bam_fields = "loremipsum2\t16\tdummyI\t3\t0\t19M\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*"
        >>> sample_line = f"{sample_second_bam_fields}\tMM:Z:T+g?,1;\tML:B:C,50\tXA:Z:loremipsum2_dummyI_2_21_rev"
        >>> modBAM_record_to_detect(sample_line, "sample.fa", code='g')
        '>loremipsum2 dummyI 2 21 rev\n3\t0.197265\tTAGCTA'
        >>> modBAM_record_to_detect(sample_line, "sample.fa", reverse_shift=0, code='g')
        '>loremipsum2 dummyI 2 21 rev\n8\t0.197265\tATCGTT'
        >>> modBAM_record_to_detect(sample_line, reverse_shift=0, code='g')
        '>loremipsum2 unmapped 0 19 unstranded\n12\t0.197265\tTAGCTA'
        >>> modBAM_record_to_detect(sample_line, code='g')
        '>loremipsum2 unmapped 0 19 unstranded\n7\t0.197265\tAACGAT'
    """

    # process modBAM line (NOTE: although we've specified a threshold, the thresholded probabilities are not used here)
    mod_bam_record = ModBamRecordProcessor(0.5, code, allow_non_na_mode=True, use_xa_tag=True)
    mod_bam_record.process_modbam_line(modbam_line)

    # set unmapped mode if fasta file is not provided or the record is unmapped
    is_unmapped_mode = (not fasta_file_name) or mod_bam_record.is_unmapped

    # if fasta file is not provided and the record is unmapped, return empty string
    if fasta_file_name and is_unmapped_mode:
        return ""

    # get normal formats for thymidine probabilities and positions
    probability_normal_format = convert_probabilities_from_modBAM_to_normal(
        mod_bam_record.raw_probability_modbam_format)

    if is_unmapped_mode:
        thymidine_pos_normal_format = mod_bam_record.fwd_seq_thymidine_coordinates
        start = 0
        end = len(mod_bam_record.fwd_seq)
    else:
        thymidine_pos_normal_format = mod_bam_record.fwd_seq_reference_coordinates
        start = mod_bam_record.start
        end = mod_bam_record.end

    if mod_bam_record.is_rev and not is_unmapped_mode:
        probability_normal_format = list(reversed(probability_normal_format))
        thymidine_pos_normal_format = list(reversed(thymidine_pos_normal_format))

    # get reference sequence
    # NOTE: We may need to get 6mers beyond the end, so we get 5 more bases than necessary
    if not is_unmapped_mode:
        fasta_file = pysam.FastaFile(fasta_file_name)
        ref_seq = fasta_file.fetch(mod_bam_record.contig, start,
                                   min(end + 5, fasta_file.get_reference_length(mod_bam_record.contig))).upper()
    else:
        ref_seq = mod_bam_record.fwd_seq + "-----"

    # function that gets sixmers
    def get_six_mer(pos0, seq):
        return seq[pos0:pos0 + 6]

    # make and return detect data
    detect_body = "\n".join("\t".join(
        (
            str(k[0]) if (not mod_bam_record.is_rev) else str(k[0] - reverse_shift),
            str(k[1])[0: 8],
            get_six_mer(k[0] - start if (not mod_bam_record.is_rev)
                        else k[0] - start - reverse_shift, ref_seq)
        )
    ) for k in filter(lambda x: x[0] != -1, zip(thymidine_pos_normal_format, probability_normal_format)))

    return mod_bam_record.detect_header(is_unmapped_mode) + ("\n" + detect_body if detect_body else "")


def modBAM_record_windowed_average(modbam_line: str, window_size: int = 0, threshold: float = 0.5,
                                   operation_per_win: Callable[[list[int]], any] = lambda x: sum(x) / max(1, len(x)),
                                   code: str = 'T', use_xa_tag: bool = False, rev: bool = False, base: str = 'T',
                                   force_missing: bool = False) \
        -> tuple[str, list[float]]:
    r""" Window thresholded modification data and report one value per window per read id.
    Caution: Windows are over thymidines where the data is available,
             i.e. bases marked as missing in modbam are not included in window size or averaging calculations.
            Let's say there are six thymidines, with data missing at the second thymidine,
            and modification probabilities at the first, third,..., sixth thymidines.
            If the user requests averages over windows of size 2, the function reports the average of the first & third
            and fourth & fifth thymidines, omitting the sixth as it occupies a window of size 1.
            If missing data is treated as NA, which this function does not do,
            then the output will be NA, average of third & fourth, average of fifth & sixth,
            as data is missing at the second thymidine.

    Args:
        modbam_line: one line of data from modBAM file. must contain MM, ML tags and only one type of modification.
        window_size: (default 0) report one average per window of given size. default value of 0 reports average
            over entire read.
        threshold: (default 0.5) threshold above (below) which base is regarded as modified (unmodified)
        operation_per_win: (defaults to get mean) function that operates on list of integers and produces some output
        code: (default T) code of modification
        use_xa_tag: (default F) if T, report detect header format but with underscores (readid_contig_start_end_orn)
            instead of read ids
        rev: (default F) reverse the analogue-incorporation data before windowing
        base: (default T) base to be considered for modification
        force_missing: (default F) if True, missing bases are treated as missing and not as unmodified,
            irrespective of blank/dot/question mark tags used in mod bam

    Returns:
        tuple of read id and list of values per window

    Examples:
        >>> sample_first_bam_fields = "loremipsum\t0\t*\t0\t0\t*\t*\t0\t10\tTGTGTGTGTG\t*"
        >>> sample_first_bam_fields_aligned = "loremipsum\t0\taabb\t4\t0\t10M\t*\t0\t10\tTGTGTGTGTG\t*"
        >>> sample_second_bam_fields = "loremipsum\t0\t*\t0\t0\t*\t*\t0\t10\tTGTTTTTTTTTTGTGTGTTTTTTGT\t*"
        >>> sample_first_bam_fields_meth = "loremipsum\t0\t*\t0\t0\t*\t*\t0\t10\tCGCGCGCGCG\t*"
        >>> sample_first_bam_fields_meth_aligned = "loremipsum\t0\tccdd\t5\t0\t10M\t*\t0\t10\tCGCGCGCGCG\t*"
        >>> sample_line = f"{sample_first_bam_fields}\tMM:Z:T+T?,1,0,0,0;\tML:B:C,10,11,12,13\tblahblah"
        >>> modBAM_record_windowed_average(sample_line, 1, 0.01)
        ('loremipsum', [1.0, 1.0, 1.0, 1.0])
        >>> sample_line = f"{sample_first_bam_fields}\tMM:Z:T+T,1,0,0,0;\tML:B:C,10,11,12,13\tblahblah"
        >>> modBAM_record_windowed_average(sample_line, 1, 0.01)
        ('loremipsum', [0.0, 1.0, 1.0, 1.0, 1.0])
        >>> sample_line = f"{sample_first_bam_fields}\tMM:Z:T+T.,1,0,0,0;\tML:B:C,10,11,12,13\tblahblah"
        >>> modBAM_record_windowed_average(sample_line, 1, 0.01)
        ('loremipsum', [0.0, 1.0, 1.0, 1.0, 1.0])
        >>> sample_line = f"{sample_first_bam_fields}\tMM:Z:T+T.,1,0,0,0;\tML:B:C,10,11,12,13\tblahblah"
        >>> modBAM_record_windowed_average(sample_line, 1, 0.01, force_missing=True)
        ('loremipsum', [1.0, 1.0, 1.0, 1.0])
        >>> sample_line = f"{sample_second_bam_fields}\tMM:Z:T+T.,0,10,0,6;\tML:B:C,10,11,12,13\tblahblah"
        >>> modBAM_record_windowed_average(sample_line, 4, 0.01)
        ('loremipsum', [0.25, 0.0, 0.25, 0.25, 0.25])
        >>> sample_line = f"{sample_first_bam_fields_aligned}\tMm:Z:T+g?,1,0,0,0;\tMl:B:C,10,11,12,13\tXA:Z:loremipsum_aabb_4_14_fwd"
        >>> modBAM_record_windowed_average(sample_line, 1, 0.01, use_xa_tag=True, code='g')
        ('loremipsum_aabb_4_14_fwd', [1.0, 1.0, 1.0, 1.0])
        >>> sample_line_meth = f"{sample_first_bam_fields_meth_aligned}\tMm:Z:C+m?,1,0,0,0;\tMl:B:C,10,11,12,13\tXA:Z:loremipsum_ccdd_5_15_fwd"
        >>> modBAM_record_windowed_average(sample_line_meth, 1, 0.01, use_xa_tag=True, code='m', base='C')
        ('loremipsum_ccdd_5_15_fwd', [1.0, 1.0, 1.0, 1.0])
        >>> modBAM_record_windowed_average(sample_line, 2, 11/256, code='g')
        ('loremipsum', [0.5, 1.0])
        >>> modBAM_record_windowed_average(sample_line, 3, 11/256, lambda x: sum(x), code='g')
        ('loremipsum', [2])
        >>> modBAM_record_windowed_average(sample_line, 3, 11/256, lambda x: sum(x), code='g', rev=True)
        ('loremipsum', [3])
        >>> modBAM_record_windowed_average(sample_line, 0, 11/256, code='g')
        ('loremipsum', [0.75])
        >>> modBAM_record_windowed_average(sample_line, 0, 11/256, lambda x: len(x), code='g')
        ('loremipsum', [4])
    """

    window_size = int(window_size)

    modbam_record = ModBamRecordProcessor(threshold, code, use_xa_tag, True, base, force_missing)
    modbam_record.process_modbam_line(modbam_line)

    if modbam_record.has_data():

        read_id = modbam_record.read_id
        probability_modbam_format = modbam_record.probability_modbam_format

        if window_size <= 0:
            return read_id, [operation_per_win(probability_modbam_format)]

        # reverse data if requested
        if rev:
            probability_modbam_format = list(reversed(probability_modbam_format))

        n_data = len(probability_modbam_format)
        return read_id, [operation_per_win(probability_modbam_format[k: k + window_size]) for k in
                         range(0, n_data, window_size) if k + window_size <= n_data]

    elif modbam_record.has_read():
        return modbam_record.read_id, []
    else:
        raise ValueError("Problem with modbam format. "
                         "If multiple modifications are present, then do not use.")


def modBAM_add_XR_tag(modbam_line: str) -> str:
    """ Add an XR tag to a .mod.bam data line.

    An XR tag is something we use in our workflow.
    It is a custom tag and is not a part of the BAM standard.
    We need it because we think it makes it easier to query a .mod.bam file using pre-existing tools.
    Let's say we want data from a read id in the region chrI:1000-2000.
    Previously, we have to query that region and then search through all the read ids for a match.
    With the XR tag, we can query that region and for an XR value simultaneously.
    As we use only the first 7 characters for the XR tag, it's possible that two reads will have the same tag value.
    But, the probability is very small, about 1 in 250 million (16^7).
    (the calculation above doesn't even take into account the probability that reads with the same tag will
    overlap on the genome, which will decrease the probability even further).

    Args:
        modbam_line: one line of data from modBAM file which does not contain an XR tag.

    Returns:
        modbam_line (with new lines removed) + "\t" + "XR:i:f(read_id)" where
          f(x) is a function that converts first 7 characters of read_id into an integer

    """

    read_id, non_read_id_part = modbam_line.split("\t", 1)

    if len(read_id) == 0 or "XR:i:" in non_read_id_part:
        raise ValueError("Problem with input!")

    # function to convert UUID to int
    def uuid_first7_char_to_int(x):
        return int(x[0:7], 16)

    return modbam_line.strip() + f"\tXR:i:{uuid_first7_char_to_int(read_id)}\n"


def modBAM_add_XA_tag(modbam_line: str) -> str:
    r""" Add an XA tag to a reference-anchored .mod.bam data line.

    An XA tag is something we use in our workflow. It looks like
    XA:Z:readid_contig_start_end_orn where orn=fwd/rev and the other meanings are clear.

    Args:
        modbam_line: one line of data from modBAM file which does not contain an XA tag.

    Returns:
        modbam_line (with new lines removed) + "\t" + "XA:z:string" where string is as explained above.

    Examples:
        >>> sample_line = "loremipsum\t0\tdummyI\t2\t0\t12M\t*\t0\t12\tGCTAGCTATCGT\t*\n"
        >>> modBAM_add_XA_tag(sample_line)
        'loremipsum\t0\tdummyI\t2\t0\t12M\t*\t0\t12\tGCTAGCTATCGT\t*\tXA:Z:loremipsum_dummyI_1_13_fwd\n'
        >>> sample_line = "loremipsum\t0\tdummyI\t2\t0\t12M\t*\t0\t12\tGCTAGCTATCGT\t*\tXA:Z:loremipsum_dummyI_1_13_fwd\n"
        >>> modBAM_add_XA_tag(sample_line)
        Traceback (most recent call last):
          ...
        ValueError: Problem with input!
        >>> sample_line = "loremipsum2\t16\tdummyII\t3\t0\t19M\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*\n"
        >>> modBAM_add_XA_tag(sample_line)
        'loremipsum2\t16\tdummyII\t3\t0\t19M\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*\tXA:Z:loremipsum2_dummyII_2_21_rev\n'
        >>> sample_line = "loremipsum2\t4\tdummyII\t3\t0\t19M\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*\n"
        >>> modBAM_add_XA_tag(sample_line)
        'loremipsum2\t4\tdummyII\t3\t0\t19M\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*\n'
        >>> sample_line = "loremipsum2\t16\t*\t0\t0\t*\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*\n"
        >>> modBAM_add_XA_tag(sample_line)
        'loremipsum2\t16\t*\t0\t0\t*\t*\t0\t19\tCTAGCTATCGTTTCTGTGA\t*\n'
    """

    read_id, flag, contig, start, quality, cigar, rest = modbam_line.split("\t", 6)

    if len(read_id) == 0 or "XA:Z:" in rest:
        raise ValueError("Problem with input!")

    flag = int(flag)
    orn = "rev" if flag & 16 == 16 else "fwd"
    start = int(start) - 1
    is_unmapped = (contig == "*") or (start == 0) or (cigar == "*") or (flag & 4 == 4)

    if not is_unmapped:
        end = cigar_to_ref_to_query_tbl(cigar, start)[-1][0]
        return modbam_line.strip() + f"\tXA:Z:{read_id}_{contig}_{start}_{end}_{orn}\n"
    else:
        return modbam_line.strip() + "\n"


def modBAM_add_XP_tag(modbam_line: str) -> str:
    """ Add an XP tag to a .mod.bam data line.

    XP tag is used for pulse data. It looks like
    XP:Z:plot_only_reads_with_pulse

    Args:
        modbam_line: one line of data from modBAM file which does not contain an XP tag.

    Returns:
        modbam_line (with new lines removed) + "\t" + "XP:Z:string" where string is as explained above.

    """

    read_id, rest = modbam_line.split("\t", 1)

    if len(read_id) == 0 or "XP:Z:" in rest:
        raise ValueError("Problem with input!")

    return modbam_line.strip() + f"\tXP:Z:plot_only_self\n"


def cigar_to_ref_to_query_tbl(cigar_str, ref_start):
    """Use cigar string to map query coord to ref coord

    Args:
        cigar_str (string): cigar string
        ref_start (int): starting position along reference

    Returns:
        list of two-element tuples (ref coord, query coord)

    Examples:
        >>> cigar_to_ref_to_query_tbl("100M",0)
        [(0, 0), (100, 100)]
        >>> cigar_to_ref_to_query_tbl("21=",0)
        [(0, 0), (21, 21)]
        >>> cigar_to_ref_to_query_tbl("42=5D7M12D",2)
        [(2, 0), (44, 42), (49, 42), (56, 49), (68, 49)]
    """

    latest_pair = (ref_start, 0)
    ref_q_map_tbl = [latest_pair]
    num = 0

    for k in cigar_str:
        if '0' <= k <= '9':
            num = num * 10 + int(k)
        else:
            latest_num = num
            num = 0
            if k == 'M' or k == '=' or k == 'X':
                mov_tup = (1, 1)
            elif k == 'D' or k == 'N':
                mov_tup = (1, 0)
            elif k == 'I' or k == 'S':
                mov_tup = (0, 1)
            elif k == 'P' or k == 'H':
                # mov_tup = (0,0)
                # not gonna add redundant entries
                continue
            else:
                raise ValueError('Invalid cigar string!')

            latest_pair = (latest_pair[0] + latest_num * mov_tup[0],
                           latest_pair[1] + latest_num * mov_tup[1])

            ref_q_map_tbl.append(latest_pair)

    return ref_q_map_tbl


def process_detect_index(x: str) -> tuple[str, str, int, int, str]:
    """ Process detect index

    Args:
        x: detect index

    Returns:
        read_id, contig, start, end, strand (=fwd or rev)

    Examples:
        >>> process_detect_index("mm_dummyI_1_13_fwd")
        ('mm', 'dummyI', 1, 13, 'fwd')
        >>> process_detect_index("mm_dummyI_dummyII_100_4000_rev_L_5000_70000")
        ('mm', 'dummyI_dummyII', 100, 4000, 'rev')
    """
    return process_fork_index(x + "_L_0_10", False)[0:5]


def process_fork_index(x: str, check_if_fork_within_read: bool = True) -> tuple[str, str, int, int, str, str, int, int]:
    """ Process fork index

    Args:
        x: fork index
        check_if_fork_within_read: (default True) if True, check if fork coordinates are within read coordinates

    Returns:
        read_id, contig, start, end, strand (=fwd or rev), fork_type (=L or R), fork_start, fork_end

    Examples:
        >>> process_fork_index("mm_dummyI_dummyII_100_400000_rev_L_5000_70000")
        ('mm', 'dummyI_dummyII', 100, 400000, 'rev', 'L', 5000, 70000)
        >>> process_fork_index("mm_dummyI_1_130000_fwd_R_90000_100000")
        ('mm', 'dummyI', 1, 130000, 'fwd', 'R', 90000, 100000)
    """
    # split by underscore
    b = x.split("_")

    # find the index where strand information is stored
    strand_index = [i for i, s in enumerate(b) if s in ["fwd", "rev"]][0]

    # assign start, end, strand, startFork, endFork, forkType
    start = int(b[strand_index - 2])
    end = int(b[strand_index - 1])
    strand = b[strand_index]
    startFork = int(b[strand_index + 2])
    endFork = int(b[strand_index + 3])
    forkType = b[strand_index + 1]

    # check that fork type is either L or R
    assert forkType in ["L", "R"]

    # check that start < end and startFork < endFork
    assert start < end
    assert startFork < endFork

    # ensure that startFork and endFork are within start and end
    assert (not check_if_fork_within_read) or (start <= startFork <= endFork <= end)

    # return
    return b[0], "_".join(b[1:strand_index - 2]), start, end, strand, forkType, startFork, endFork
