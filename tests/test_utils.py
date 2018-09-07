import os
import json
import unittest

from worker import utils

class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestUtils, cls).setUpClass()
        cls.cfg = prepare_cfg()
        
    def test_conn(self):
        self.assertEqual(utils.check_internet_conn(), True)

    def test_string_to_dict(self):
        self.assertDictEqual(
            utils.string_to_dict("(Man)3(GlcNAc)2"),
            {"(GlcNAc)": 2, "(Man)": 3}
        )
        self.assertEqual(
            utils.string_to_dict("(HexNAc)14 (NeuAc)4 (NeuGc)4 (Sulph)1"),
            {"(HexNAc)": 14, "(NeuAc)": 4, "(NeuGc)": 4, "(Sulph)": 1}
        )
        self.assertDictEqual(
            utils.string_to_dict("(Hex)10 (HexNAc)5 (NeuAc)5 + (Man)3(GlcNAc)2"),
            {"(Hex)": 10, "(HexNAc)": 5, "(GlcNAc)": 2, "(NeuAc)": 5, "(Man)": 3}
        )
        self.assertDictEqual(
            utils.string_to_dict("(Hex)3 (HexNAc)9 (Deoxyhexose)1 (NeuAc)3 (NeuGc)2 (Pent)3 (Phos)1 (Sulph)1 + (Man)3(GlcNAc)2"),
            {"(Man)": 3, "(Hex)": 3, "(HexNAc)": 9, "(GlcNAc)": 2, "(NeuAc)": 3, "(NeuGc)": 2, "(Pent)": 3, "(Sulph)": 1, "(Phos)": 1, "(Deoxyhexose)": 1}
        )
 
    def test_truncate_str(self):
        self.assertEqual(utils.truncate_str("(Hex)4 (HexNAc)2"), "H4N2")
        self.assertEqual(utils.truncate_str("(Hex)10 (HexNAc)5 (NeuAc)5 + (Man)3(GlcNAc)2"), "H13N7S5")
        self.assertEqual(
            utils.truncate_str("(Hex)2 (HexNAc)9 (Deoxyhexose)2 (NeuAc)2 (NeuGc)3 (KDN)3 (Pent)3 (Phos)1 (Sulph)1 (HexA)1 + (Man)3(GlcNAc)2"),
            "H5N11F2P3S2Sg3KDN3Sulph1Phos1HexA1"
        )
    
    def test_calc_theor_mono_mass(self):
        self.assertAlmostEqual(utils.calc_theor_mono_mass({"(GlcNAc)": 2, "(Man)": 3}, TestUtils.cfg), 911.335 , places=3)        
        self.assertAlmostEqual(utils.calc_theor_mono_mass({"(HexNAc)": 4, "(NeuAc)": 4, "(NeuGc)": 3, "(Sulph)": 1, "(Hex)": 5}, TestUtils.cfg), 3807.2087, places=3)
        # ProA labeled
        self.assertAlmostEqual(utils.calc_theor_mono_mass({"(GlcNAc)": 2, "(Man)": 3}, TestUtils.cfg, reducing_end="ProA"), 1130.507 , places=2)  
        self.assertAlmostEqual(
            utils.calc_theor_mono_mass(
                {"(HexNAc)": 4, "(NeuAc)": 4, "(NeuGc)": 3, "(Sulph)": 1, "(Hex)": 5}, TestUtils.cfg, reducing_end="ProA"),
                4026.381 , places=2
        )
        
    def test_calc_theor_avg_mass(self):
        self.assertAlmostEqual(utils.calc_theor_avg_mass({"(GlcNAc)": 2, "(Man)": 3}, TestUtils.cfg), 911.83987 , places=3)        
        self.assertAlmostEqual(utils.calc_theor_avg_mass({"(HexNAc)": 4, "(NeuAc)": 4, "(NeuGc)": 3, "(Sulph)": 1, "(Hex)": 5}, TestUtils.cfg), 3809.38237, places=3)
        # ProA labeled
        self.assertAlmostEqual(utils.calc_theor_avg_mass({"(GlcNAc)": 2, "(Man)": 3}, TestUtils.cfg, reducing_end="ProA"), 1131.165, places=2)  
        self.assertAlmostEqual(
            utils.calc_theor_avg_mass(
                {"(HexNAc)": 4, "(NeuAc)": 4, "(NeuGc)": 3, "(Sulph)": 1, "(Hex)": 5}, TestUtils.cfg, reducing_end="ProA"),
                4028.708, places=2
        )

    def test_validate_filename(self):
        self.assertEqual(utils.validate_filename("proper_file_name"), "proper_file_name")
        self.assertEqual(utils.validate_filename("pr0p3r_f1l3_n4m3"), "pr0p3r_f1l3_n4m3")
        with self.assertRaises(ValueError):
            utils.validate_filename("_@!2unsafe_file_name213-5112493@#$%^&")

def prepare_cfg():
    with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_cfg.json")), "r") as f:
        cfg = json.load(f)
    return cfg