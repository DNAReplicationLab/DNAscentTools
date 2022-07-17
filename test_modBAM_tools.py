import unittest
import os
import pysam
from modbampy import ModBam
from modBAM_tools import get_gaps_in_base_pos, \
    convert_detect_into_detect_stream, \
    convert_data_per_T_to_modBAM_fmt, \
    convert_dnascent_detect_to_modBAM_file, \
    get_mod_counts_per_interval, \
    get_read_data_from_modBAM


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
            "#Genome dummy.fa",
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
        with open("dummy.fa", "w") as dummyFa:
            dummyFa.write(cls.fakeFaFile)

        pysam.faidx("dummy.fa")

        # program info
        cls.pgInfo = {
            'PN': 'convert_detect_to_modBAM',
            'ID': 'convert_detect_to_modBAM',
            'VN': 'dummyVersion'
        }

        # make modBAM file
        convert_dnascent_detect_to_modBAM_file(
            convert_detect_into_detect_stream(cls.fakeDetect.split("\n")),
            'dummy.bam', 'T',
            True,
            pg_info=cls.pgInfo)

        # index file
        pysam.index("dummy.bam")

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
                    "#Genome dummy.fa",
                    "#Index dummy.index"
                ],
                "refFasta": "dummy.fa"
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

        expected_op = [128, 129]

        op = convert_data_per_T_to_modBAM_fmt(
            [0.5, 0.5 + 1 / 256]
        )

        self.assertEqual(expected_op, op)

    def test_modBAM_retrieval_1(self):
        """ Test that modBAM retrieval works """

        # get mod counts
        t = get_mod_counts_per_interval('dummy.bam',
                                        [('dummyIII', 26, 30),
                                         ('dummyIII', 31, 35),
                                         ('dummyIII', 26, 35)],
                                        'T', 'T', 0.5, 0.5)

        # check we picked up only modified bases in fwd strand
        self.assertEqual([[0, 0, 0], [0, 0, 0], [0, 0, 0], [1, 1, 2]],
                         list(list(k) for k in t))

        # with high thresholds, we revert to unmodified
        t = get_mod_counts_per_interval('dummy.bam',
                                        [('dummyIII', 26, 30),
                                         ('dummyIII', 31, 35),
                                         ('dummyIII', 26, 35)],
                                        'T', 'T', 0.99, 0.99)

        self.assertEqual([[0, 0, 0], [1, 1, 2], [0, 0, 0], [0, 0, 0]],
                         list(list(k) for k in t))

        # test that individual sites look ok
        # note: fffffff, the first 7 letters of read id, become
        # 16 ** 7 - 1 when converted to decimal

        with ModBam("dummy.bam") as bam:
            site_data = [
                (
                    p.rpos, p.qpos, p.qual,
                    p.query_name, p.strand, p.mstrand, p.cbase, p.mbase
                ) for k in
                bam.reads(
                    'dummyII', 3, 36, tag_name='XR',
                    tag_value=16 ** 7 - 1
                )
                for p in k.mod_sites]

            site_tuple = ("fffffff1-10d2-49cb-8ca3-e8d48979001b",
                          "-", 0, 'T', 'T')

            self.assertEqual(site_data,
                             [
                                 (15, 12, 3, *site_tuple),
                                 (16, 13, 3, *site_tuple),
                                 (19, 16, 4, *site_tuple),
                                 (22, 19, 3, *site_tuple),
                                 (23, 20, 182, *site_tuple),
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
                             "dummy.bam",
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
                    "dummy.bam",
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
        header_str = pysam.view("-H", "dummy.bam")

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

    @classmethod
    def tearDownClass(cls):
        """ Delete temporary files """
        os.system("rm dummy.bam")
        os.system("rm dummy.bam.bai")
        os.system("rm dummy.fa")
        os.system("rm dummy.fa.fai")


if __name__ == '__main__':
    unittest.main()
