import os
import unittest
import json
from bs4 import BeautifulSoup

from worker.worker import GlycomodWorker as GW
from worker.data_types import GlycomodComposition, SubmittedMass

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_CFG_PATH = os.path.join(THIS_DIR, 'test_cfg.json')

class TestGlycomodWorker(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestGlycomodWorker, cls).setUpClass()
        cls.cfg = prepare_cfg()
    
    def test_parsing_small(self):
        worker = GW(cfg=TestGlycomodWorker.cfg)
        # mocking html
        worker.soup = prepare_mock_html('minimal_test.html')
        worker._parse_gm_html()
        self.assertEqual(
            worker.parsed_data, 
            [
                [
                    'User mass: 911.30', 
                    'Adduct ([M+H]+): 1.00727', 
                    'Derivative mass (Free reducing end): 18.0105546', 
                    '892.317-0.034(Hex)3 (HexNAc)2  ', '1 structure'
                    ], 
                [
                    'User mass: 1057.33', 
                    'Adduct ([M+H]+): 1.00727', 
                    'Derivative mass (Free reducing end): 18.0105546', 
                    '1038.375-0.062(Hex)3 (HexNAc)2 (Deoxyhexose)1  ', 
                    '1 structure found.'
                    ]
            ]
        )

    def test_create_glycan_objects(self):
        worker = GW(cfg=TestGlycomodWorker.cfg)
        worker.parsed_data = prepare_mock_parsed_data()
        worker._create_glycan_objects()
        self.assertIsInstance(worker.compositions[0], SubmittedMass)
        self.assertIsInstance(worker.compositions[0].glycomod_structures[0], GlycomodComposition)


def prepare_mock_html(file_name):
    file_path = os.path.join(THIS_DIR, 'test_html', file_name)  
    with open(file_path, "r") as f:
        html_string = f.read()
    data = BeautifulSoup(html_string, 'html5lib')
    return data

def prepare_mock_parsed_data():
    return  [
                [
                    'User mass: 1057.33', 
                    'Adduct ([M+H]+): 1.00727', 
                    'Derivative mass (Free reducing end): 18.0105546', 
                    '1038.375-0.062(Hex)3 (HexNAc)2 (Deoxyhexose)1  ', 
                    '1 structure found.'
                    ]
            ]


def prepare_cfg():
    with open(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_cfg.json")), "r") as f:
        cfg = json.load(f)
    return cfg