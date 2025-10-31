import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd

from pipeline import ProteinPipeline
from config import PipelineConfig


class TestProteinPipeline(unittest.TestCase):
    def setUp(self):
        self.config = PipelineConfig.from_args(
            base="test_base",
            input_template="tests/fixtures/test_input.xlsx",
            output_excel="test_output.xlsx",
        )
        self.pipeline = ProteinPipeline(config=self.config)

    def test_initialization(self):
        self.assertEqual(self.pipeline.cfg, self.config)
        self.assertIsNone(self.pipeline.ids)
        self.assertFalse(self.pipeline.verbose)

    def test_load_uniprot_ids(self):
        with patch("builtins.open", MagicMock()):
            with patch("pandas.read_excel", return_value=pd.DataFrame({"UniProt ID": ["P12345", "Q67890"]})):
                ids = self.pipeline._load_uniprot_ids()
                self.assertEqual(ids, ["P12345", "Q67890"])

    @patch("pipeline.fetch_uniprot_json")
    @patch("pipeline.parse_uniprot_data")
    @patch("pipeline.find_related_uniprots")
    @patch("pipeline.get_alphafold_pdb_info")
    @patch("pipeline.download_pdb_structure")
    @patch("pipeline.fix_structure")
    @patch("pipeline.detect_gaps_and_discard_coils")
    @patch("pipeline.parse_helix_sheet_annotations")
    @patch("pipeline.get_pdb_details")
    @patch("pipeline.parse_small_molecules")
    @patch("pipeline.download_alphafold_pdb")
    @patch("pipeline.get_alphafold_details")
    def test_process_uniprot(self, mock_get_alphafold_details, mock_download_alphafold_pdb, mock_parse_small_molecules, mock_get_pdb_details, mock_parse_helix_sheet_annotations, mock_detect_gaps_and_discard_coils, mock_fix_structure, mock_download_pdb_structure, mock_get_alphafold_pdb_info, mock_find_related_uniprots, mock_parse_uniprot_data, mock_fetch_uniprot_json):
        mock_fetch_uniprot_json.return_value = {}
        mock_parse_uniprot_data.return_value = {
            "gene": "TEST",
            "protein_name": "Test Protein",
            "sequence": "ABC",
            "pdb_ids": ["1ABC"],
        }
        mock_find_related_uniprots.return_value = []
        mock_get_alphafold_pdb_info.return_value = "AF-P12345-F1"
        mock_download_pdb_structure.return_value = (
            "test_base/pdb_structures_exp/TEST_P12345/TEST_P12345_1ABC.pdb", "RCSB")
        mock_fix_structure.return_value = True
        mock_detect_gaps_and_discard_coils.return_value = []
        mock_parse_helix_sheet_annotations.return_value = {}
        mock_get_pdb_details.return_value = {}
        mock_parse_small_molecules.return_value = []
        mock_download_alphafold_pdb.return_value = "test_base/pdb_structures_pred/TEST_P12345/TEST_P12345_AF-P12345-F1.pdb"
        mock_get_alphafold_details.return_value = {}

        self.pipeline._process_uniprot("P12345")

        self.assertEqual(len(self.pipeline._uniprot_records), 1)
        self.assertEqual(self.pipeline._uniprot_records[0]["Gene Name"], self.pipeline._hyperlink(
            "https://www.genecards.org/cgi-bin/carddisp.pl?gene=TEST", "TEST"))
        self.assertEqual(len(self.pipeline._pdb_detail_records), 2)


if __name__ == "__main__":
    unittest.main()
