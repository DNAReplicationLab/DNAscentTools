import unittest
import os
import pysam
import itertools
from modBAM_tools import get_gaps_in_base_pos,\
        convert_detect_into_detect_stream, \
        convert_data_per_T_to_modBAM_fmt, \
        convert_dnascent_detect_to_modBAM_file, \
        get_mod_counts_per_interval

class TestDetectToModBAMSuite(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Make fake detect and fasta data

        Args:
            None

        Returns:
            None
        """
        
        cls.fakeDetect = ("# dummy data adapted from dnascent docs\n"
            "#Genome dummy.fa\n"
            "#Index dummy.index\n"
            ">5d10eb9a-aae1-4db8-8ec6-7ebb34d32575 dummyI 9 17 fwd\n"
            "9\t0.017496\tTCGTTT\n"
            "12\t0.029483\tTTTCTG\n"
            "13\t0.039008\tTTCTGT\n"
            "14\t0.031474\tTCTGTG\n"
            "16\t0.026997\tTGTGAG\n"
            ">a4f36092-b4d5-47a9-813e-c22c3b477a0c dummyIII 23 71 fwd\n"
            "26\t0.866907\tTCAAAT\n"
            "31\t0.947935\tTCCACA\n"
            "50\t0.014683\tTGGGAG\n"
            "62\t0.186812\tTAACGG\n"
            "70\t0.934850\tTTATTG\n"
            ">c6785e1f-10d2-49cb-8ca3-e8d48979001b dummyII 3 36 rev\n"
            "10\t0.012874\tTCTCTA\n"
            "11\t0.012428\tCTCTAA\n"
            "14\t0.016811\tTAACGA\n"
            "17\t0.013372\tCGACCA\n"
            "18\t0.013836\tGACCAA\n"
            )

        cls.seq1 = "AGCTAGCTATCGTTTCTGTGAG"
        cls.seq2 = "AGCTAGCTAGTCTCTAACGACCAA"
        cls.seq3 = (
            "CCACACCACACCCACACACCCACACATCAAATCCACACCACACCACACCC"
            "TGGGAGCCACCATAACGGCCTTATTG"
            )

        cls.fakeFaFile = (">dummyI\n"
            f"{cls.seq1}\n"
            ">dummyII\n"
            f"{cls.seq2}\n"
            ">dummyIII\n"
            f"{cls.seq3[0:50]}\n"
            f"{cls.seq3[50:]}\n"
            )

        # make fake fasta file and index it
        with open("dummy.fa", "w") as dummyFa:
            dummyFa.write(cls.fakeFaFile)

        pysam.faidx("dummy.fa")

    def test_get_gaps_in_base_pos(self):
        """ Test if finding gaps in T works """

        self.assertEqual([1],get_gaps_in_base_pos([1,5], 'ATATATA', 'T'))
        self.assertEqual([2],get_gaps_in_base_pos([1,5], 'ATTTATA', 'T'))
        self.assertEqual([0,0,1],get_gaps_in_base_pos([1,5,9,14], 
            'XGXXXGXXXGXXGXGXG', 'G'))

    def test_convert_detect_into_detect_stream(self):
        """ Test if iterating over detect records works """

        expectedDetectStream = [
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
                "posOnRef": [9, 12, 13, 14, 16],
                "probBrdU": [0.017496, 0.029483, 0.039008, 
                    0.031474, 0.026997],
                "sixMerOnRef": ["TCGTTT", "TTTCTG", "TTCTGT", 
                    "TCTGTG", "TGTGAG"]
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
                "readID": "c6785e1f-10d2-49cb-8ca3-e8d48979001b",
                "refContig": "dummyII",
                "refStart": 3,
                "refEnd": 36,
                "strand": "rev",
                "posOnRef": [10, 11, 14, 17, 18],
                "probBrdU": [0.012874, 0.012428, 0.016811, 
                    0.013372, 0.013836],
                "sixMerOnRef": ["TCTCTA", "CTCTAA", "TAACGA", 
                    "CGACCA", "GACCAA"]
            }
        ]

        self.assertEqual(
            expectedDetectStream,
            list(
                convert_detect_into_detect_stream(
                    self.fakeDetect.split("\n")
                ))
            )

    def test_convert_data_per_T_to_modBAM_str(self):
        """ Test converting data to modBAM strings """

        expectedOp = [128,129]

        op = convert_data_per_T_to_modBAM_fmt(
                [0.5, 0.5 + 1/256]
            )

        self.assertEqual(expectedOp, op)

    def test_modBAM_making(self):
        """ Test that modBAM making works """

        # make modBAM file
        convert_dnascent_detect_to_modBAM_file(
            convert_detect_into_detect_stream(self.fakeDetect.split("\n")), 
            'dummy.bam', 'T', 
            True,
            pgInfo = {
                'ID': 'convert_detect_to_modBAM',
                'PN': 'convert_detect_to_modBAM',
                'VN': 'dummy'
            })

        # index file
        pysam.index("dummy.bam")

        # get mod counts
        t = get_mod_counts_per_interval('dummy.bam', 
                    [('dummyIII', 26, 30),
                     ('dummyIII', 31, 35),
                     ('dummyIII', 26, 35)], 
                    'T', 'T', 0.5, 0.5)

        # check we picked up only modified bases in fwd strand
        self.assertEqual([[0,0,0],[0,0,0],[0,0,0],[1,1,2]],
            list(list(k) for k in t) )

        # with high thresholds, we revert to unmodified
        t = get_mod_counts_per_interval('dummy.bam', 
                    [('dummyIII', 26, 30),
                     ('dummyIII', 31, 35),
                     ('dummyIII', 26, 35)], 
                    'T', 'T', 0.99, 0.99)

        self.assertEqual([[0,0,0],[1,1,2],[0,0,0],[0,0,0]],
            list(list(k) for k in t) )

    @classmethod
    def tearDownClass(cls):
        """ Delete temporary files """
        os.system("rm dummy.bam")
        os.system("rm dummy.bam.bai")
        os.system("rm dummy.fa")
        os.system("rm dummy.fa.fai")

if __name__ == '__main__':
    unittest.main()
