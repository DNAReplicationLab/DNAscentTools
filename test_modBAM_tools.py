import unittest
from modBAM_tools import get_gaps_in_base_pos,\
        convert_detect_into_detect_stream, \
            convert_data_per_T_to_modBAM_fmt

class TestDetectToModBAMSuite(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Make a fake detect file

        Args:
            None

        Returns:
            None
        """
        
        cls.fakeDetect = ("# dummy data adapted from dnascent docs\n"
            "#Genome dummy.fasta\n"
            "#Index dummy.index\n"
            ">5d10eb9a-aae1-4db8-8ec6-7ebb34d32575 chrXIII 72319 72327 fwd\n"
            "72319\t0.017496\tTCGTTT\n"
            "72322\t0.029483\tTTTCTG\n"
            "72323\t0.039008\tTTCTGT\n"
            "72324\t0.031474\tTCTGTG\n"
            "72326\t0.026997\tTGTGAG\n"
            ">a4f36092-b4d5-47a9-813e-c22c3b477a0c chrXVI 899273 899321 fwd\n"
            "899276\t0.866907\tTCAAAT\n"
            "899281\t0.947935\tTCCACA\n"
            "899300\t0.014683\tTGGGAG\n"
            "899312\t0.186812\tTAACGG\n"
            "899320\t0.934850\tTTATTG\n"
            ">c6785e1f-10d2-49cb-8ca3-e8d48979001b chrXIII 74003 81176 rev\n"
            "74010\t0.012874\tTCTCTA\n"
            "74011\t0.012428\tCTCTAA\n"
            "74014\t0.016811\tTAACGA\n"
            "74017\t0.013372\tCGACCA\n"
            "74018\t0.013836\tGACCAA\n"
            )

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
                    "#Genome dummy.fasta",
                    "#Index dummy.index"
                ],
                "refFasta": "dummy.fasta"
            },
            {
                "readID": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "refContig": "chrXIII",
                "refStart": 72319,
                "refEnd": 72327,
                "strand": "fwd",
                "posOnRef": [72319, 72322, 72323, 72324, 72326],
                "probBrdU": [0.017496, 0.029483, 0.039008, 
                    0.031474, 0.026997],
                "sixMerOnRef": ["TCGTTT", "TTTCTG", "TTCTGT", 
                    "TCTGTG", "TGTGAG"]
            },
            {
                "readID": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
                "refContig": "chrXVI",
                "refStart": 899273,
                "refEnd": 899321,
                "strand": "fwd",
                "posOnRef": [899276, 899281, 899300, 899312, 899320],
                "probBrdU": [0.866907, 0.947935, 0.014683, 
                    0.186812, 0.934850],
                "sixMerOnRef": ["TCAAAT", "TCCACA", "TGGGAG", 
                    "TAACGG", "TTATTG"]
            },
            {
                "readID": "c6785e1f-10d2-49cb-8ca3-e8d48979001b",
                "refContig": "chrXIII",
                "refStart": 74003,
                "refEnd": 81176,
                "strand": "rev",
                "posOnRef": [74010,74011,74014,74017,74018],
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

if __name__ == '__main__':
    unittest.main()