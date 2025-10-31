from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest
from Bio.PDB import PDBParser

import oriented
import structure_repair
from structure_repair import (
    _measure_thickness,
    orient_with_memembed,
    orient_with_tmdet,
)


def _write_linear_membrane_model(pdb_path: Path, spacing: float = 4.0) -> None:
    residues = [-16 + spacing * i for i in range(9)]
    lines = []
    atom_id = 1
    for idx, x in enumerate(residues, start=1):
        lines.append(
            f"ATOM  {atom_id:5d}  N   LEU A{idx:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 50.00           N"
        )
        atom_id += 1
        lines.append(
            f"ATOM  {atom_id:5d}  CA  LEU A{idx:4d}    {x + 1.2:8.3f}{0.5:8.3f}{0.2:8.3f}  1.00 50.00           C"
        )
        atom_id += 1
        lines.append(
            f"ATOM  {atom_id:5d}  C   LEU A{idx:4d}    {x + 2.2:8.3f}{0.8:8.3f}{-0.3:8.3f}  1.00 50.00           C"
        )
        atom_id += 1
    lines.append("TER")
    lines.append("END")
    pdb_path.write_text("\n".join(lines) + "\n")


@pytest.fixture()
def linear_membrane(tmp_path: Path) -> Path:
    pdb = tmp_path / "linear.pdb"
    _write_linear_membrane_model(pdb)
    return pdb


def _load_coords(path: Path) -> np.ndarray:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    coords = [atom.get_coord() for atom in structure.get_atoms() if getattr(
        atom, "element", atom.get_name()[0]).upper() != "H"]
    return np.asarray(coords)


def test_memembed_orientation_outputs_file(linear_membrane: Path, tmp_path: Path):
    out_file = tmp_path / "memembed_oriented.pdb"
    meta = orient_with_memembed(linear_membrane, out_file)
    assert meta is not None
    assert out_file.is_file()
    coords = _load_coords(out_file)
    thickness = _measure_thickness(coords)
    assert pytest.approx(thickness, rel=1e-3) == meta["thickness"]
    assert 20.0 <= meta["thickness"] <= 45.0


def test_tmdet_enforces_thickness_window(linear_membrane: Path, tmp_path: Path):
    out_file = tmp_path / "tmdet_oriented.pdb"
    meta = orient_with_tmdet(linear_membrane, out_file)
    assert meta is not None
    assert out_file.is_file()
    assert 25.0 <= meta["thickness"] <= 40.0
    assert abs(meta["midplane_offset"]) < 1.0


def test_fallback_uses_memembed_before_tmdet(monkeypatch, tmp_path: Path):
    order: list[str] = []

    def _fail(name):
        def _handler(*args, **kwargs):
            order.append(name)
            raise RuntimeError(name)

        return _handler

    def _memembed(dest: Path):
        order.append("memembed")
        Path(dest).write_text("HEADER\nEND\n")
        return {"method": "Memembed"}

    monkeypatch.setattr(structure_repair, "orient_with_ppm", _fail("ppm"))
    monkeypatch.setattr(structure_repair, "orient_with_memembed",
                        lambda pdb, dest: _memembed(dest))
    monkeypatch.setattr(structure_repair, "orient_with_tmdet", _fail("tmdet"))
    monkeypatch.setattr(
        structure_repair, "orient_with_principal_axes", _fail("axes"))

    cfg = SimpleNamespace(ppm_path=None, prefer_native_sources=False)
    dummy = tmp_path / "dummy.pdb"
    dummy.write_text("HEADER\nEND\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    new_path, method = oriented._orient_with_fallback(
        "TEST", dummy, out_dir, cfg)

    assert method == "Memembed"
    assert Path(new_path).is_file()
    assert order == ["memembed"]


def test_ensure_repair_dependencies_python_guard(monkeypatch):
    monkeypatch.setattr(structure_repair, "PDBFixer", None, raising=False)
    monkeypatch.setattr(structure_repair, "PDBFile", None, raising=False)
    monkeypatch.setattr(
        structure_repair, "_REPAIR_AVAILABLE", None, raising=False)
    monkeypatch.setattr(
        structure_repair, "_REPAIR_WARNING_EMITTED", False, raising=False)
    monkeypatch.setattr(
        structure_repair, "_INSTALL_ATTEMPTED", False, raising=False)

    fake_sys = SimpleNamespace(version_info=(3, 13, 0), executable="python")
    monkeypatch.setattr(structure_repair, "sys", fake_sys, raising=False)
    monkeypatch.setattr(structure_repair.subprocess, "check_call",
                        lambda *a, **k: pytest.fail("should not install"))

    assert structure_repair._ensure_repair_dependencies(
        auto_install=True) is False
    assert structure_repair._REPAIR_AVAILABLE is False
