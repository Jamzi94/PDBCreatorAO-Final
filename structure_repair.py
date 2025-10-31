"""Structure repair module using PDBFixer, with additional membrane protein orientation tools."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import subprocess
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from datetime import datetime, timedelta
from bs4 import BeautifulSoup

try:
    from openmm.app import PDBFile  # type: ignore
    from pdbfixer import PDBFixer  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    PDBFile = None  # type: ignore
    PDBFixer = None  # type: ignore

try:  # pragma: no cover - optional dependency
    from Bio.PDB import PDBIO, PDBParser  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    PDBIO = None  # type: ignore
    PDBParser = None  # type: ignore


log = logging.getLogger(__name__)


def _create_http_session(*, retries: int = 3, backoff: float = 0.8) -> requests.Session:
    """Return a requests session configured with retries and friendly headers."""
    session = requests.Session()
    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=("GET", "POST"),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy,
                          pool_connections=10, pool_maxsize=20)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(
        {
            "User-Agent": "MembraneProteinPipeline/1.0 (+https://github.com/yourusername/yourrepo)",
            "Accept": "application/xml,application/json,text/html,application/xhtml+xml;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5",
            "Accept-Encoding": "gzip, deflate",
            "Connection": "keep-alive",
            "Upgrade-Insecure-Requests": "1",
        }
    )
    return session


session = _create_http_session()


def _download_file(url: str, out_path: Path, *, timeout: tuple[int, int] = (5, 30)) -> bool:
    """Download a file to disk using the shared HTTP session."""
    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(response.content)
        return True
    except Exception as exc:  # pragma: no cover - defensive log
        log.warning("Download failed for %s: %s", url, exc)
        return False


def _download_file_with_size_check(
    url: str,
    out_path: Path,
    *,
    timeout: tuple[int, int] = (5, 30),
    max_size_mb: int = 20,
) -> bool:
    """Download a file ensuring it does not exceed the configured size limit."""
    try:
        with session.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            content_length = response.headers.get("content-length")
            if content_length:
                size_mb = int(content_length) / (1024 * 1024)
                if size_mb > max_size_mb:
                    log.warning(
                        "File too large (%.2f MB > %d MB): %s",
                        size_mb,
                        max_size_mb,
                        url,
                    )
                    return False

            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "wb") as handle:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        handle.write(chunk)

        return out_path.is_file() and out_path.stat().st_size > 0
    except Exception as exc:  # pragma: no cover - defensive log
        log.warning("Download failed for %s: %s", url, exc)
        if out_path.is_file():
            out_path.unlink()
        return False


class OrientationCache:
    """Persist orientation metadata for reuse across runs."""

    def __init__(self, cache_dir: Path = Path(".orientation_cache")) -> None:
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(exist_ok=True)
        self.cache_file = cache_dir / "orientation_cache.json"
        self.cache_data = self._load_cache()

    def _load_cache(self) -> Dict[str, Dict[str, Any]]:
        if self.cache_file.is_file():
            try:
                with self.cache_file.open("r", encoding="utf-8") as handle:
                    return json.load(handle)
            except Exception:  # pragma: no cover - ignore corrupt cache
                return {}
        return {}

    def _save_cache(self) -> None:
        with self.cache_file.open("w", encoding="utf-8") as handle:
            json.dump(self.cache_data, handle, indent=2)

    def get_cached_result(self, pdb_id: str, service: str) -> Optional[Dict[str, Any]]:
        cache_key = f"{pdb_id}_{service}"
        entry = self.cache_data.get(cache_key)
        if not entry:
            return None
        try:
            timestamp = datetime.fromisoformat(entry["timestamp"])
        except Exception:  # pragma: no cover - invalid timestamp
            return None
        if timestamp < datetime.now() - timedelta(days=7):
            return None
        return entry.get("result")

    def set_cached_result(self, pdb_id: str, service: str, result: Dict[str, Any]) -> None:
        cache_key = f"{pdb_id}_{service}"
        self.cache_data[cache_key] = {
            "timestamp": datetime.now().isoformat(),
            "result": result,
        }
        self._save_cache()


_orientation_cache = OrientationCache()


def _ensure_repair_dependencies(auto_install: bool) -> bool:
    """Ensure PDBFixer/OpenMM are available, installing them when allowed."""
    global PDBFixer, PDBFile, _REPAIR_WARNING_EMITTED, _REPAIR_AVAILABLE, _INSTALL_ATTEMPTED

    if PDBFixer is not None and PDBFile is not None:
        _REPAIR_AVAILABLE = True
        return True

    if _REPAIR_AVAILABLE is False:
        return False

    major, minor = sys.version_info[:2]
    if (major, minor) >= (3, 13):
        if not _REPAIR_WARNING_EMITTED:
            log.warning(
                "PDBFixer/OpenMM wheels are unavailable for Python %s.%s; use the provided conda environment (python 3.11) or install manually.",
                major,
                minor,
            )
            _REPAIR_WARNING_EMITTED = True
        _REPAIR_AVAILABLE = False
        return False

    if auto_install and not _INSTALL_ATTEMPTED:
        _INSTALL_ATTEMPTED = True
        log.info(
            "Attempting automatic installation of PDBFixer/OpenMM for structure repair.")
        try:
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "openmm", "pdbfixer"])
        except Exception as exc:  # pragma: no cover - defensive log
            log.error("Automatic installation of repair tools failed: %s", exc)
            _REPAIR_AVAILABLE = False
            return False

    try:
        from openmm.app import PDBFile as _PDBFile  # type: ignore
        from pdbfixer import PDBFixer as _PDBFixer  # type: ignore
    except ImportError as exc:
        if not _REPAIR_WARNING_EMITTED:
            log.warning(
                "Required tools (PDBFixer/OpenMM) not installed. Skipping structure repair.")
            _REPAIR_WARNING_EMITTED = True
        _REPAIR_AVAILABLE = False
        if auto_install and _INSTALL_ATTEMPTED:
            log.debug("Repair import failure after installation attempt: %s", exc)
        return False

    globals()["PDBFile"] = _PDBFile
    globals()["PDBFixer"] = _PDBFixer
    _REPAIR_AVAILABLE = True
    return True


def fix_structure(input_path: Path, output_path: Path, auto_install: bool = False) -> bool:
    """Repair a PDB structure using PDBFixer."""
    if not _ensure_repair_dependencies(auto_install):
        return False

    try:
        try:
            from structure_analysis import parse_small_molecules  # type: ignore
        except ImportError:  # pragma: no cover - fallback when parser missing
            parse_small_molecules = None  # type: ignore

        ligands_before: list[str] = []
        if parse_small_molecules is not None:
            try:
                ligands_before = parse_small_molecules(str(input_path))
            except Exception as exc:  # pragma: no cover - defensive log
                log.debug(
                    "Failed to parse small molecules for %s: %s", input_path, exc)

        fixer = PDBFixer(str(input_path))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        fixer.removeHeterogens(keepWater=False)

        missing_residues = len(getattr(fixer, "missingResidues", []))
        missing_atoms = len(getattr(fixer, "missingAtoms", []))
        nonstandard = len(getattr(fixer, "nonstandardResidues", []))
        removed_ligands = ligands_before

        log.info(
            "Repair summary for %s :: missing_residues=%d missing_atoms=%d nonstandard=%d removed_ligands=%s",
            input_path.name,
            missing_residues,
            missing_atoms,
            nonstandard,
            ", ".join(removed_ligands) if removed_ligands else "none",
        )

        output_path = output_path.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions,
                              handle, keepIds=True)

        log.info("Repaired structure saved to %s", output_path)
        return True
    except Exception as exc:  # pragma: no cover - defensive
        log.error("Failed to repair %s with PDBFixer: %s",
                  input_path.name, exc)
        return False


HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"}


def _load_orientation_payload(pdb_path: Path):
    if PDBParser is None or PDBIO is None:
        log.error(
            "Biopython is required to orient %s but is not installed", pdb_path)
        return None
    if not pdb_path.is_file():
        log.error("PDB file not found for orientation: %s", pdb_path)
        return None

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    heavy_coords = []
    hydrophobic_coords = []
    for atom in structure.get_atoms():
        element = getattr(atom, "element", "").upper()
        name = atom.get_name().upper()
        if element == "H" or name.startswith("H"):
            continue
        coord = atom.get_coord()
        heavy_coords.append(coord)
        resname = atom.get_parent().get_resname().strip().upper()
        if resname in HYDROPHOBIC_RESIDUES:
            hydrophobic_coords.append(coord)
    return structure, np.asarray(heavy_coords, dtype=float), np.asarray(hydrophobic_coords, dtype=float)


def _principal_axes(coords: np.ndarray):
    if coords.size == 0:
        raise ValueError(
            "No coordinates available for principal axes computation")
    center = coords.mean(axis=0)
    centered = coords - center
    cov = centered.T @ centered / max(len(centered), 1)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    if np.linalg.det(eigenvectors) < 0:
        eigenvectors[:, -1] *= -1
    return center, eigenvalues, eigenvectors


def _build_rotation(eigenvectors: np.ndarray, order: tuple[int, int, int] = (0, 1, 2)) -> np.ndarray:
    rotation = np.column_stack([eigenvectors[:, idx] for idx in order])
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    return rotation


def _apply_orientation(structure, center: np.ndarray, rotation: np.ndarray, out_path: Path) -> None:
    for atom in structure.get_atoms():
        vec = atom.get_coord() - center
        atom.set_coord(vec @ rotation)
    io = PDBIO()
    io.set_structure(structure)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_path))


def _measure_thickness(coords: np.ndarray) -> float:
    if coords.size == 0:
        return 0.0
    z_vals = coords[:, 2]
    return float(z_vals.max() - z_vals.min())


def orient_with_pdbtm_cached(pdb_id: str, src_path: Path, out_path: Path) -> Optional[Dict[str, Any]]:
    """Use existing PDBTM/OPM files as reference for alignment.

    Args:
        pdb_id: PDB ID of the target structure
        src_path: Path to the PDB file to orient
        out_path: Where to save the oriented structure

    Returns:
        Dict with method, reference, and metadata
    """
    from oriented import _align_to_reference

    # Extract gene/uniprot info from the path
    parent_dir = src_path.parent.name
    if "_" not in parent_dir:
        log.warning("Cannot extract gene/uniprot from path %s", src_path)
        return None

    gene_name, uniprot_id = parent_dir.split("_", 1)
    group_dir = src_path.parent

    # Look for OPM/PDBTM references in this group's directory
    ref_file = None
    ref_source = None

    for source in ("OPM", "PDBTM"):
        for f in group_dir.glob("*.pdb"):
            if f.stem != pdb_id and source.lower() in f.name.lower():
                ref_file = f
                ref_source = source
                break
        if ref_file:
            break

    if not ref_file:
        log.info("No native reference found for %s in %s", pdb_id, group_dir)
        return None

    try:
        rms_before, rms_after = _align_to_reference(
            ref_file, src_path, out_path)
        log.info(
            "Aligned %s to %s reference %s (RMS: %.3f → %.3f)",
            pdb_id, ref_source, ref_file.stem, rms_before, rms_after
        )
        return {
            "method": f"Native_{ref_source}",
            "aligned_to": str(ref_file),
            "metadata": {
                "rms_before": rms_before,
                "rms_after": rms_after,
                "reference_source": ref_source,
            }
        }
    except Exception as exc:
        log.error("Failed to align %s to %s: %s", pdb_id, ref_file, exc)
        return None


def orient_with_memembed(pdb_path: Path, out_path: Path) -> Optional[Dict[str, Any]]:
    """Approximate Memembed orientation using hydrophobic principal axes."""
    payload = _load_orientation_payload(pdb_path)
    if payload is None:
        return None
    structure, heavy_coords, hydrophobic_coords = payload
    source_coords = hydrophobic_coords if hydrophobic_coords.shape[0] >= 3 else heavy_coords
    if source_coords.shape[0] < 3:
        log.warning(
            "Not enough hydrophobic atoms to orient %s via Memembed", pdb_path.name)
        return None

    center, eigenvalues, eigenvectors = _principal_axes(source_coords)
    rotation = _build_rotation(eigenvectors, (1, 2, 0))
    rotated = (source_coords - center) @ rotation
    thickness = _measure_thickness(rotated)

    _apply_orientation(structure, center, rotation, out_path)

    return {
        "method": "Memembed",
        "eigenvalues": eigenvalues.tolist(),
        "thickness": thickness,
    }


def orient_with_tmdet(
    pdb_path: Path,
    out_path: Path,
    thickness_window: tuple[float, float] = (25.0, 40.0),
) -> Optional[Dict[str, Any]]:
    """Heuristic TMDET-style orientation with bilayer thickness validation."""
    payload = _load_orientation_payload(pdb_path)
    if payload is None:
        return None
    structure, heavy_coords, hydrophobic_coords = payload
    source_coords = hydrophobic_coords if hydrophobic_coords.shape[0] >= 3 else heavy_coords
    if source_coords.shape[0] < 3:
        log.warning(
            "Not enough atoms to orient %s via TMDET heuristic", pdb_path.name)
        return None

    center, eigenvalues, eigenvectors = _principal_axes(source_coords)
    candidate_orders = [(1, 2, 0), (2, 0, 1), (0, 1, 2)]
    best = None
    for order in candidate_orders:
        rotation = _build_rotation(eigenvectors, order)
        rotated = (source_coords - center) @ rotation
        thickness = _measure_thickness(rotated)
        score = abs(thickness - 32.0)
        if best is None or score < best[0]:
            best = (score, rotation, thickness, order, rotated)

    if best is None:
        return None

    _, rotation, thickness, order, rotated_coords = best
    min_thk, max_thk = thickness_window
    if thickness < min_thk or thickness > max_thk:
        log.info(
            "TMDET heuristic rejected for %s: thickness %.2f ? outside [%0.1f, %0.1f]",
            pdb_path.name,
            thickness,
            min_thk,
            max_thk,
        )
        return None

    _apply_orientation(structure, center, rotation, out_path)

    midplane_offset = float(
        rotated_coords[:, 2].mean()) if rotated_coords.size else 0.0

    return {
        "method": "TMDET",
        "eigenvalues": eigenvalues.tolist(),
        "thickness": thickness,
        "axis_order": order,
        "midplane_offset": midplane_offset,
    }


def orient_with_principal_axes(pdb_path: Path, out_path: Path) -> Optional[Dict[str, Any]]:
    """Orient a structure by aligning its principal axes with the coordinate frame."""
    payload = _load_orientation_payload(pdb_path)
    if payload is None:
        return None
    structure, heavy_coords, _ = payload
    if heavy_coords.shape[0] < 3:
        log.warning("Not enough heavy atoms to orient %s", pdb_path.name)
        return None

    center, eigenvalues, eigenvectors = _principal_axes(heavy_coords)
    rotation = _build_rotation(eigenvectors, (0, 1, 2))
    _apply_orientation(structure, center, rotation, out_path)

    return {
        "method": "PrincipalAxes",
        "eigenvalues": eigenvalues.tolist(),
    }


def orient_with_ppm(pdb_path: Path, out_path: Path, ppm_executable: str) -> Optional[Dict[str, Any]]:
    """Run a local PPM calculation for membrane protein orientation."""
    executable = Path(ppm_executable)
    if not executable.exists():
        log.error("PPM executable not found at %s", ppm_executable)
        return None

    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            [str(executable), "-i", str(pdb_path), "-o", str(out_path)],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        log.error("PPM execution failed for %s: %s", pdb_path, exc)
        return None
    except Exception as exc:  # pragma: no cover - defensive
        log.error("Error running PPM for %s: %s", pdb_path, exc)
        return None

    output = result.stdout or ""
    if "error" in output.lower():
        log.error("PPM reported an error for %s: %s", pdb_path, output.strip())
        return None

    if not out_path.is_file():
        log.error("PPM did not produce an output file for %s", pdb_path)
        return None

    return {"method": "PPM_local", "stdout": output.strip()}
