#!/usr/bin/env python3
"""Basic functionality test for PDBCreatorAO pipeline."""

import sys
from pathlib import Path

# Test 1: Import all core modules
print("Test 1: Importing core modules...")
try:
    import pipeline
    import uniprot
    import alignment
    import alphafold
    import oriented
    import structure_repair
    import structure_analysis
    import pdb_download
    import pdb_details
    import config
    print("✅ All modules imported successfully")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

# Test 2: Verify pipeline configuration
print("\nTest 2: Testing pipeline configuration...")
try:
    from config import PipelineConfig
    cfg = PipelineConfig.from_args(base="/tmp/test_pipeline")
    assert cfg.base_dir == Path("/tmp/test_pipeline")
    print("✅ Pipeline configuration works")
except Exception as e:
    print(f"❌ Configuration failed: {e}")
    sys.exit(1)

# Test 3: Test pipeline initialization
print("\nTest 3: Testing pipeline initialization...")
try:
    from pipeline import ProteinPipeline
    pipe = ProteinPipeline(config=cfg, ids=["P12345"], verbose=False)
    assert pipe.ids == ["P12345"]
    print("✅ Pipeline initialization works")
except Exception as e:
    print(f"❌ Pipeline initialization failed: {e}")
    sys.exit(1)

# Test 4: Test helper functions
print("\nTest 4: Testing helper functions...")
try:
    from pipeline import extract_pdb_from_cell
    from unittest.mock import Mock
    
    cell = Mock()
    cell.hyperlink = None
    cell.value = "1ABC"
    result = extract_pdb_from_cell(cell)
    assert result == "1ABC"
    print("✅ Helper functions work")
except Exception as e:
    print(f"❌ Helper functions failed: {e}")
    sys.exit(1)

# Test 5: Test argument parser
print("\nTest 5: Testing CLI argument parser...")
try:
    from pipeline import build_arg_parser
    parser = build_arg_parser()
    args = parser.parse_args(['--ids', 'P12345', 'Q67890', '--window', '10'])
    assert args.ids == ['P12345', 'Q67890']
    assert args.window == 10
    print("✅ CLI argument parser works")
except Exception as e:
    print(f"❌ Argument parser failed: {e}")
    sys.exit(1)

print("\n" + "="*60)
print("ALL BASIC FUNCTIONALITY TESTS PASSED ✅")
print("The code is functional and ready for use.")
print("="*60)
