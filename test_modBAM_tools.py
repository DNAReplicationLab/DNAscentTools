import unittest
import os
import pysam
from modBAM_tools import get_gaps_in_base_pos, \
    convert_detect_into_detect_stream, \
    convert_data_per_T_to_modBAM_fmt, \
    convert_dnascent_detect_to_modBAM_file, \
    get_mod_counts_per_interval, \
    get_read_data_from_modBAM
from modBAM_tools_additional import get_raw_data_from_modBAM, \
    ModBamRecordProcessor, parse_modBAM_modification_information, \
    convert_bed_to_detect_stream


class TestConvertBedToDetectStream(unittest.TestCase):

    def setUp(self):
        """
        Set up the test environment by creating a temporary BED file.
        This method creates a BED file with the following content:
        # BED file content:
        # sample comment
        track something
        chr1    100    101    name1    0.5    +
        chr1    150    151    name1    0.7    +
        chr1    200    201    name2    0.9    -
        chr1    250    251    name2    0.8    -
        The file is saved as 'test.bed' in the current directory.
        """
        self.bed_file_content = (
            "# sample comment\n"
            "track something\n"
            "chr1\t100\t101\tname1\t0.5\t+\n"
            "chr1\t150\t151\tname1\t0.7\t+\n"
            "chr1\t200\t201\tname2\t0.9\t-\n"
            "chr1\t250\t251\tname2\t0.8\t-\n"
        )
        self.bed_file_path = "test.bed"
        with open(self.bed_file_path, "w") as bed_file:
            bed_file.write(self.bed_file_content)

    def tearDown(self):
        """ Remove the temporary BED file """
        os.remove(self.bed_file_path)

    def test_convert_bed_to_detect_stream(self):
        """ Test the convert_bed_to_detect_stream function """
        expected_output = [
            {"comments": ["converted from bed file test.bed"]},
            {
                "readID": "name1",
                "refContig": "chr1",
                "refStart": 100,
                "refEnd": 151,
                "strand": "fwd",
                "posOnRef": [100, 150],
                "probBrdU": [0.5, 0.7],
                'probEdU': [],
                "sixMerOnRef": ["NNNNNN", "NNNNNN"]
            },
            {
                "readID": "name2",
                "refContig": "chr1",
                "refStart": 200,
                "refEnd": 251,
                "strand": "rev",
                "posOnRef": [200, 250],
                "probBrdU": [0.9, 0.8],
                'probEdU': [],
                "sixMerOnRef": ["NNNNNN", "NNNNNN"]
            }
        ]

        result = convert_bed_to_detect_stream(self.bed_file_path)
        self.assertEqual(result, expected_output)


class TestDetectToModBAMSuite(unittest.TestCase):
    pgInfo = None
    fakeDetect = None
    fakeFaFile = None
    seqContig1 = None
    seqContig2 = None
    seqContig3 = None
    comments = None

    @classmethod
    def setUpClass(cls):
        """ Make fake detect and fasta data

        Returns:
            None
        """

        cls.comments = [
            "# dummy data adapted from dnascent docs",
            "#Genome sample.fa",
            "#Index dummy.index"
        ]

        cls.seqContig1 = "AGCTAGCTATCGTTTCTGTGAG"
        cls.seqContig2 = "GGGGGGGGGGTCTCTAACGACCAAGGGGGGGGGGGGGGGGGGGGGGGG"
        cls.seqContig3 = (
            "CCACACCACACCCACACACCCACACATCAAATCCACACCACACCACACCC"
            "TGGGAGCCACCATAACGGCCTTATTG"
        )

        cls.fakeDetect = (f"{cls.comments[0]}\n"
                          f"{cls.comments[1]}\n"
                          f"{cls.comments[2]}\n"
                          ">5d10eb9a-aae1-4db8-8ec6-7ebb34d32575 dummyI 9 17 fwd\n"
                          f"9\t0.017496\t{cls.seqContig1[9:15]}\n"
                          f"12\t0.029483\t{cls.seqContig1[12:18]}\n"
                          f"13\t0.039008\t{cls.seqContig1[13:19]}\n"
                          f"16\t0.026997\t{cls.seqContig1[16:22]}\n"
                          ">a4f36092-b4d5-47a9-813e-c22c3b477a0c dummyIII 23 71 fwd\n"
                          f"26\t0.866907\t{cls.seqContig3[26:32]}\n"
                          f"31\t0.947935\t{cls.seqContig3[31:37]}\n"
                          f"50\t0.014683\t{cls.seqContig3[50:56]}\n"
                          f"62\t0.186812\t{cls.seqContig3[62:68]}\n"
                          f"70\t0.934850\t{cls.seqContig3[70:76]}\n"
                          ">fffffff1-10d2-49cb-8ca3-e8d48979001b dummyII 3 36 rev\n"
                          f"10\t0.012874\t{cls.seqContig2[10:16]}\n"
                          f"11\t0.012428\t{cls.seqContig2[11:17]}\n"
                          f"14\t0.016811\t{cls.seqContig2[14:20]}\n"
                          f"17\t0.013372\t{cls.seqContig2[17:23]}\n"
                          f"18\t0.713836\t{cls.seqContig2[18:24]}\n"
                          )

        cls.query1 = cls.seqContig1[9: 17]
        cls.query2 = cls.seqContig2[3: 36]
        cls.query3 = cls.seqContig3[23: 71]

        cls.fakeFaFile = (">dummyI\n"
                          f"{cls.seqContig1}\n"
                          ">dummyII\n"
                          f"{cls.seqContig2}\n"
                          ">dummyIII\n"
                          f"{cls.seqContig3[0:50]}\n"
                          f"{cls.seqContig3[50:]}\n"
                          )

        # make fake fasta file and index it
        with open("sample.fa", "w") as dummyFa:
            dummyFa.write(cls.fakeFaFile)

        pysam.faidx("sample.fa")

        # program info
        cls.pgInfo = {
            'PN': 'convert_detect_to_modBAM',
            'ID': 'convert_detect_to_modBAM',
            'VN': 'dummyVersion'
        }

        # make modBAM file
        convert_dnascent_detect_to_modBAM_file(
            convert_detect_into_detect_stream(cls.fakeDetect.split("\n")),
            'sample.bam', 'T',
            True,
            pg_info=cls.pgInfo)

        # make modSAM file
        convert_dnascent_detect_to_modBAM_file(
            convert_detect_into_detect_stream(cls.fakeDetect.split("\n")),
            'sample.sam', 'T',
            False,
            pg_info=cls.pgInfo)

        # index file
        pysam.index("sample.bam")

    def test_get_gaps_in_base_pos(self):
        """ Test if finding gaps in T works """

        self.assertEqual([1], get_gaps_in_base_pos([1, 5], 'ATATATA', 'T'))
        self.assertEqual([2], get_gaps_in_base_pos([1, 5], 'ATTTATA', 'T'))
        self.assertEqual([1, 0, 0, 1, 1], get_gaps_in_base_pos([-1, 1, 5, 9, 14, 17],
                                                               'GGXXXGXXXGXXGXGXG', 'G'))

        # are erroneous inputs caught?
        self.assertRaises(ValueError, get_gaps_in_base_pos, [1, 5, 9, 14],
                          'XTXXXGXXXGXXGXGXG', 'G')
        self.assertRaises(ValueError, get_gaps_in_base_pos, [1, 5, 9, 18],
                          'XGXXXGXXXGXXGXGXG', 'G')
        self.assertRaises(ValueError, get_gaps_in_base_pos, [-2, 5, 9, 14],
                          'XGXXXGXXXGXXGXGXG', 'G')
        self.assertRaises(ValueError, get_gaps_in_base_pos, [5, 1, 9, 14],
                          'XGXXXGXXXGXXGXGXG', 'G')

    def test_convert_detect_into_detect_stream(self):
        """ Test if iterating over detect records works """

        expected_detect_stream = [
            {
                "comments": [
                    "# dummy data adapted from dnascent docs",
                    "#Genome sample.fa",
                    "#Index dummy.index"
                ],
                "refFasta": "sample.fa"
            },
            {
                "readID": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "refContig": "dummyI",
                "refStart": 9,
                "refEnd": 17,
                "strand": "fwd",
                "posOnRef": [9, 12, 13, 16],
                "probBrdU": [0.017496, 0.029483, 0.039008,
                             0.026997],
                "probEdU": [],
                "sixMerOnRef": ["TCGTTT", "TTTCTG", "TTCTGT",
                                "TGTGAG"]
            },
            {
                "readID": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
                "refContig": "dummyIII",
                "refStart": 23,
                "refEnd": 71,
                "strand": "fwd",
                "posOnRef": [26, 31, 50, 62, 70],
                "probBrdU": [0.866907, 0.947935, 0.014683,
                             0.186812, 0.934850],
                "probEdU": [],
                "sixMerOnRef": ["TCAAAT", "TCCACA", "TGGGAG",
                                "TAACGG", "TTATTG"]
            },
            {
                "readID": "fffffff1-10d2-49cb-8ca3-e8d48979001b",
                "refContig": "dummyII",
                "refStart": 3,
                "refEnd": 36,
                "strand": "rev",
                "posOnRef": [10, 11, 14, 17, 18],
                "probBrdU": [0.012874, 0.012428, 0.016811,
                             0.013372, 0.713836],
                "probEdU": [],
                "sixMerOnRef": ["TCTCTA", "CTCTAA", "TAACGA",
                                "CGACCA", "GACCAA"]
            }
        ]

        self.assertEqual(
            expected_detect_stream,
            list(
                convert_detect_into_detect_stream(
                    self.fakeDetect.split("\n")
                ))
        )

    def test_convert_data_per_T_to_modBAM_str(self):
        """ Test converting data to modBAM strings """

        expected_op = [128, 129, 0, 255]

        op = convert_data_per_T_to_modBAM_fmt(
            [0.5, 0.5 + 1 / 256, 0, 0.9999999]
        )

        self.assertEqual(expected_op, op)

    def test_modBAM_retrieval_1(self):
        """ Test that modBAM retrieval works """

        # get mod counts
        t = get_mod_counts_per_interval('sample.bam',
                                        [('dummyIII', 26, 30),
                                         ('dummyIII', 31, 35),
                                         ('dummyIII', 26, 35)],
                                        'T', 'T', 0.5, 0.5)

        # check we picked up only modified bases in fwd strand
        self.assertEqual([[0, 0, 0], [0, 0, 0], [0, 0, 0], [1, 1, 2]],
                         list(list(k) for k in t))

        # with high thresholds, we revert to unmodified
        t = get_mod_counts_per_interval('sample.bam',
                                        [('dummyIII', 26, 30),
                                         ('dummyIII', 31, 35),
                                         ('dummyIII', 26, 35)],
                                        'T', 'T', 0.99, 0.99)

        self.assertEqual([[0, 0, 0], [1, 1, 2], [0, 0, 0], [0, 0, 0]],
                         list(list(k) for k in t))

    def test_modBAM_retrieval_1a(self):
        """ Test that modBAM retrieval works """

        # test that individual sites look ok using the mod_data_to_table method
        alignments = pysam.view("-e", f"qname==\"fffffff1-10d2-49cb-8ca3-e8d48979001b\"",
                                "sample.bam", "dummyII:3-36")
        split_lines = alignments.splitlines()
        if sum([1 for k in split_lines if k]) > 1:
            raise ValueError("More than one alignment found")

        mod_bam_record = ModBamRecordProcessor(0.01, 'T', base='T', allow_non_na_mode=True)
        mod_bam_record.process_modbam_line(split_lines[0])

        site_data = [
            (
                p.ref_pos, p.fwd_seq_pos, p.mod_qual,
                p.read_id, p.ref_strand, p.mod_strand, p.can_base, p.mod_base
            ) for p in mod_bam_record.mod_data_to_table(move_parallel_top_ref_strand=True)]

        site_tuple = ("fffffff1-10d2-49cb-8ca3-e8d48979001b",
                      "-", "+", 'T', 'T')

        self.assertEqual(site_data,
                         [
                             (15, 20, (3 + 4) / (2 * 256), *site_tuple),
                             (16, 19, (3 + 4) / (2 * 256), *site_tuple),
                             (19, 16, (4 + 5) / (2 * 256), *site_tuple),
                             (22, 13, (3 + 4) / (2 * 256), *site_tuple),
                             (23, 12, (182 + 183) / (2 * 256), *site_tuple),
                         ])

    def test_modBAM_retrieval_1b(self):
        """ Test that modBAM retrieval works """

        self.assertEqual(list(get_raw_data_from_modBAM("sample.bam",
                                                       "dummyII", 3, 36, "T", "T",
                                                       "fffffff1-10d2-49cb-8ca3-e8d48979001b")),
                         [
                             ("fffffff1-10d2-49cb-8ca3-e8d48979001b", 15, (3 + 4) / (2 * 256)),
                             ("fffffff1-10d2-49cb-8ca3-e8d48979001b", 16, (3 + 4) / (2 * 256)),
                             ("fffffff1-10d2-49cb-8ca3-e8d48979001b", 19, (4 + 5) / (2 * 256)),
                             ("fffffff1-10d2-49cb-8ca3-e8d48979001b", 22, (3 + 4) / (2 * 256)),
                             ("fffffff1-10d2-49cb-8ca3-e8d48979001b", 23, (182 + 183) / (2 * 256)),
                         ])

    def test_modBAM_retrieval_2(self):
        """ Test retrieval of data from modbam file """

        # test getting read data
        for detectRecord in filter(lambda x: 'readID' in x,
                                   convert_detect_into_detect_stream(
                                       self.fakeDetect.split("\n")
                                   )):

            # expected positions and probabilities
            expected_op = zip(detectRecord["posOnRef"],
                              detectRecord["probBrdU"])

            # get read data and compare to expectation
            for k in zip(expected_op,
                         get_read_data_from_modBAM(
                             "sample.bam",
                             detectRecord["readID"],
                             detectRecord["refContig"],
                             detectRecord["refStart"],
                             detectRecord["refEnd"],
                         )
                         ):

                # if reverse strand, increment dnascent index by 5
                if detectRecord["strand"] == "rev":
                    self.assertEqual(k[0][0] + 5, k[1][0])
                else:
                    self.assertEqual(k[0][0], k[1][0])

                # in going to modBAM, resolution takes a hit, so
                # we can only ascertain equality within 1/256ths
                self.assertAlmostEqual(k[0][1], k[1][1], delta=1 / 256)

    def test_modBAM_retrieval_3(self):
        """ More tests for retrieval of data from modBAM file """

        for detectRecord in filter(lambda x: 'readID' in x,
                                   convert_detect_into_detect_stream(
                                       self.fakeDetect.split("\n")
                                   )):

            cnt = 0

            # get read data
            for _ in get_read_data_from_modBAM(
                    "sample.bam",
                    detectRecord["readID"],
                    detectRecord["refContig"],
                    detectRecord["refStart"],
                    detectRecord["refStart"] + 1,
            ):
                cnt += 1

            # at most one record should have gotten through as we are querying
            # start to start + 1
            self.assertTrue(cnt == 0 or cnt == 1)

    def test_headers_modBAM(self):
        """ Test that headers are fine in our dummy modBAM file """

        # function to extract that portion of a header corresponding to a tag
        def gets_header_given_tag(tag): return \
            [*map(
                lambda x: x[4:],
                filter(
                    lambda x: x.startswith(f"@{tag}\t"),
                    header_str.split("\n")
                )
            )]

        # get header
        header_str = pysam.view("-H", "sample.bam")

        # prepare contig and contig length portion of header and compare
        contigs = ["dummyI", "dummyII", "dummyIII"]
        contig_lens = [len(self.seqContig1), len(self.seqContig2),
                       len(self.seqContig3)]

        contig_str_list = [*map(
            lambda x: f"SN:{x[0]}\tLN:{x[1]}",
            zip(contigs, contig_lens)
        )]

        self.assertEqual(contig_str_list, gets_header_given_tag("SQ"))

        # prepare comment portion of header and compare
        self.assertEqual(
            [*map(lambda x: "comment: " + x, self.comments)],
            gets_header_given_tag("CO"))

        # compare program id part of header
        pg_info_header = "\t".join(
            f"{k[0]}:{k[1]}" for k in self.pgInfo.items()
        )
        self.assertTrue(any(map(
            lambda x: x == pg_info_header,
            gets_header_given_tag("PG"))))

    def test_alternate_reference_and_column_switch(self):
        """ Test using an alternate fasta reference for detect file """

        detect_file_contents = ("#Genome this_file_does_not_exist.fa\n"
                                ">5d10eb8a-aae1-5db8-9ec6-6ebb34d32575 fake 10 20 fwd\n"
                                "11\tAA\t0.220000\n"
                                "13\tBB\t0.420000\n"
                                "15\tCC\t0.620000\n"
                                "17\tDD\t0.820000\n")

        fasta_file_contents = ">fake\nCGCGCGCGCGCTGTGTGTGGGGGGG"

        expected_detect_stream = [
            {
                "comments": [
                    "#Genome this_file_does_not_exist.fa"
                ],
                "refFasta": "this_file_does_not_exist.fa"
            },
            {
                "readID": "5d10eb8a-aae1-5db8-9ec6-6ebb34d32575",
                "refContig": "fake",
                "refStart": 10,
                "refEnd": 20,
                "strand": "fwd",
                "posOnRef": [11, 13, 15, 17],
                "probBrdU": [0.22, 0.42, 0.62, 0.82],
                "probEdU": [],
                "sixMerOnRef": ["AA", "BB", "CC", "DD"]
            }
        ]

        self.assertEqual(
            expected_detect_stream,
            list(
                convert_detect_into_detect_stream(
                    detect_file_contents.split("\n"), switch_2_and_3=True
                ))
        )

        # make fake fasta file and index it
        with open("sample_2.fa", "w") as dummyFa:
            dummyFa.write(fasta_file_contents)

        pysam.faidx("sample_2.fa")

        # make modBAM file
        convert_dnascent_detect_to_modBAM_file(
            convert_detect_into_detect_stream(detect_file_contents.split("\n"), True),
            'sample_2.bam', 'T', True, fasta='sample_2.fa')

        # index file
        pysam.index("sample_2.bam")

    def test_parse_modBAM_modification_information(self):
        """ Test parsing of modBAM modification information """

        mod_information_1 = "MM:Z:T+m?,1,0,0;\tML:B:C,100,200,100"

        expected_result_1 = [
            {
                "base": "T",
                "mod_strand": "+",
                "mod_code": "m",
                "mode": "?",
                "pos": [1, 0, 0],
                "prob": [100, 200, 100]
            }
        ]

        mod_information_2 = "MM:Z:T+99001.,1,0,0;\tML:B:C,100,200,100"
        expected_result_2 = [{"base": "T", "mod_code": 99001, "mode": ".", "mod_strand": "+",
                              "pos": [1, 0, 0], "prob": [100, 200, 100]}]

        mod_information_3a = "MM:Z:T+de,1,0,0;C-h?,2,0;\tML:B:C,100,200,101,40,90,120,50,22"
        expected_result_3a = [{"base": "T", "mod_code": "d", "mode": ".", "mod_strand": "+",
                              "pos": [1, 0, 0], "prob": [100, 101, 90]},
                             {"base": "T", "mod_code": "e", "mode": ".", "mod_strand": "+",
                              "pos": [1, 0, 0], "prob": [200, 40, 120]},
                             {"base": "C", "mod_code": "h", "mode": "?", "mod_strand": "-",
                              "pos": [2, 0], "prob": [50, 22]}]

        mod_information_3b = "MM:Z:C-h?,2,0;T+de,1,0,0;\tML:B:C,100,200,101,40,90,120,50,22"
        expected_result_3b = [{"base": "C", "mod_code": "h", "mode": "?", "mod_strand": "-",
                               "pos": [2, 0], "prob": [100, 200]},
                              {"base": "T", "mod_code": "d", "mode": ".", "mod_strand": "+",
                               "pos": [1, 0, 0], "prob": [101, 90, 50]},
                              {"base": "T", "mod_code": "e", "mode": ".", "mod_strand": "+",
                               "pos": [1, 0, 0], "prob": [40, 120, 22]}]

        mod_information_4 = "MM:Z:T+d,1,0,0;\tML:B:C,100,200,100"
        expected_result_4 = [{"base": "T", "mod_code": "d", "mode": ".", "mod_strand": "+",
                              "pos": [1, 0, 0], "prob": [100, 200, 100]}]

        mod_information_5 = "MM:Z:T+c.,1,0,0;\tML:B:C,100,200,100"
        expected_result_5 = [{"base": "T", "mod_code": "c", "mode": ".", "mod_strand": "+",
                              "pos": [1, 0, 0], "prob": [100, 200, 100]}]

        self.assertEqual(
            expected_result_1,
            parse_modBAM_modification_information(mod_information_1)
        )
        self.assertEqual(
            expected_result_2,
            parse_modBAM_modification_information(mod_information_2)
        )
        self.assertEqual(
            expected_result_3a,
            parse_modBAM_modification_information(mod_information_3a)
        )
        self.assertEqual(
            expected_result_3b,
            parse_modBAM_modification_information(mod_information_3b)
        )
        self.assertEqual(
            expected_result_4,
            parse_modBAM_modification_information(mod_information_4)
        )
        self.assertEqual(
            expected_result_5,
            parse_modBAM_modification_information(mod_information_5)
        )

    @classmethod
    def tearDownClass(cls):
        """ Delete some temporary files """
        os.system("rm sample.fa.fai")
        os.system("rm sample_2.*")


class TestDetectToModBAMWithEdU(unittest.TestCase):
    pgInfo = None
    fakeDetect = None
    fakeFaFile = None
    seqContig1 = None
    seqContig2 = None
    seqContig3 = None
    comments = None

    @classmethod
    def setUpClass(cls):
        """ Make fake detect and fasta data with EdU

        Returns:
            None
        """

        cls.comments = [
            "# dummy data adapted from dnascent docs",
            "#Genome sample.fa",
            "#Index dummy.index"
        ]

        cls.seqContig1 = "AGCTAGCTATCGTTTCTGTGAG"
        cls.seqContig2 = "GGGGGGGGGGTCTCTAACGACCAAGGGGGGGGGGGGGGGGGGGGGGGG"
        cls.seqContig3 = (
            "CCACACCACACCCACACACCCACACATCAAATCCACACCACACCACACCC"
            "TGGGAGCCACCATAACGGCCTTATTG"
        )

        cls.fakeDetect = (f"{cls.comments[0]}\n"
                          f"{cls.comments[1]}\n"
                          f"{cls.comments[2]}\n"
                          ">5d10eb9a-aae1-4db8-8ec6-7ebb34d32575 dummyI 9 17 fwd\n"
                          f"9\t0.017496\t0.1\t{cls.seqContig1[9:15]}\n"
                          f"12\t0.029483\t0.2\t{cls.seqContig1[12:18]}\n"
                          f"13\t0.039008\t0.3\t{cls.seqContig1[13:19]}\n"
                          f"16\t0.026997\t0.4\t{cls.seqContig1[16:22]}\n"
                          ">a4f36092-b4d5-47a9-813e-c22c3b477a0c dummyIII 23 71 fwd\n"
                          f"26\t0.866907\t0.5\t{cls.seqContig3[26:32]}\n"
                          f"31\t0.947935\t0.6\t{cls.seqContig3[31:37]}\n"
                          f"50\t0.014683\t0.7\t{cls.seqContig3[50:56]}\n"
                          f"62\t0.186812\t0.8\t{cls.seqContig3[62:68]}\n"
                          f"70\t0.934850\t0.9\t{cls.seqContig3[70:76]}\n"
                          ">fffffff1-10d2-49cb-8ca3-e8d48979001b dummyII 3 36 rev\n"
                          f"10\t0.012874\t0.1\t{cls.seqContig2[10:16]}\n"
                          f"11\t0.012428\t0.2\t{cls.seqContig2[11:17]}\n"
                          f"14\t0.016811\t0.3\t{cls.seqContig2[14:20]}\n"
                          f"17\t0.013372\t0.4\t{cls.seqContig2[17:23]}\n"
                          f"18\t0.713836\t0.5\t{cls.seqContig2[18:24]}\n"
                          )

        cls.query1 = cls.seqContig1[9: 17]
        cls.query2 = cls.seqContig2[3: 36]
        cls.query3 = cls.seqContig3[23: 71]

        cls.fakeFaFile = (">dummyI\n"
                          f"{cls.seqContig1}\n"
                          ">dummyII\n"
                          f"{cls.seqContig2}\n"
                          ">dummyIII\n"
                          f"{cls.seqContig3[0:50]}\n"
                          f"{cls.seqContig3[50:]}\n"
                          )

        # make fake fasta file and index it
        with open("sample.fa", "w") as dummyFa:
            dummyFa.write(cls.fakeFaFile)

        pysam.faidx("sample.fa")

        # program info
        cls.pgInfo = {
            'PN': 'convert_detect_to_modBAM',
            'ID': 'convert_detect_to_modBAM',
            'VN': 'dummyVersion'
        }

        # make modBAM file
        convert_dnascent_detect_to_modBAM_file(
            convert_detect_into_detect_stream(cls.fakeDetect.split("\n")),
            'sample_with_edu.bam', '472552+472553',
            True,
            pg_info=cls.pgInfo)

        # index file
        pysam.index("sample_with_edu.bam")

    def test_convert_detect_into_detect_stream_with_edu(self):
        """ Test if iterating over detect records with EdU data works """
    
        expected_detect_stream = [
            {
                "comments": [
                    "# dummy data adapted from dnascent docs",
                    "#Genome sample.fa",
                    "#Index dummy.index"
                ],
                "refFasta": "sample.fa"
            },
            {
                "readID": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "refContig": "dummyI",
                "refStart": 9,
                "refEnd": 17,
                "strand": "fwd",
                "posOnRef": [9, 12, 13, 16],
                "probBrdU": [0.017496, 0.029483, 0.039008, 0.026997],
                "probEdU": [0.1, 0.2, 0.3, 0.4],
                "sixMerOnRef": ["TCGTTT", "TTTCTG", "TTCTGT", "TGTGAG"]
            },
            {
                "readID": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
                "refContig": "dummyIII",
                "refStart": 23,
                "refEnd": 71,
                "strand": "fwd",
                "posOnRef": [26, 31, 50, 62, 70],
                "probBrdU": [0.866907, 0.947935, 0.014683, 0.186812, 0.934850],
                "probEdU": [0.5, 0.6, 0.7, 0.8, 0.9],
                "sixMerOnRef": ["TCAAAT", "TCCACA", "TGGGAG", "TAACGG", "TTATTG"]
            },
            {
                "readID": "fffffff1-10d2-49cb-8ca3-e8d48979001b",
                "refContig": "dummyII",
                "refStart": 3,
                "refEnd": 36,
                "strand": "rev",
                "posOnRef": [10, 11, 14, 17, 18],
                "probBrdU": [0.012874, 0.012428, 0.016811, 0.013372, 0.713836],
                "probEdU": [0.1, 0.2, 0.3, 0.4, 0.5],
                "sixMerOnRef": ["TCTCTA", "CTCTAA", "TAACGA", "CGACCA", "GACCAA"]
            }
        ]
    
        self.assertEqual(
            expected_detect_stream,
            list(
                convert_detect_into_detect_stream(
                    self.fakeDetect.split("\n"), switch_2_and_3=True
                ))
        )

    def test_modBAM_retrieval_brdu(self):
        """ Test retrieval of data from modbam file """

        # get read data and compare to expectation
        positions = []
        probabilities = []
        for k in get_read_data_from_modBAM(
                "sample_with_edu.bam",
                "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "dummyI",
                9,
                17,
                "T",
                "472552"
        ):
            positions.append(k[0])
            probabilities.append(k[1])

        self.assertEqual(positions, [9, 12, 13, 16])

        for entries in zip(probabilities, [0.017496, 0.029483, 0.039008, 0.026997]):
            self.assertAlmostEqual(entries[0], entries[1], delta=1 / 256)

    def test_modBAM_retrieval_edu(self):
        """ Test retrieval of data from modbam file """

        # get read data and compare to expectation
        positions = []
        probabilities = []
        for k in get_read_data_from_modBAM(
                "sample_with_edu.bam",
                "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "dummyI",
                9,
                17,
                "T",
                "472553"
        ):
            positions.append(k[0])
            probabilities.append(k[1])

        self.assertEqual(positions, [9, 12, 13, 16])

        for entries in zip(probabilities, [0.1, 0.2, 0.3, 0.4]):
            self.assertAlmostEqual(entries[0], entries[1], delta=1 / 256)

    @classmethod
    def tearDownClass(cls):
        """ Delete some temporary files """
        os.system("rm sample_with_edu.bam.bai")
        os.system("rm sample_with_edu.bam")
        os.system("rm sample.fa.fai")


if __name__ == '__main__':
    unittest.main()
