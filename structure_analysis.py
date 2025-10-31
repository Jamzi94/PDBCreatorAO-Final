"""
structure_analysis.py

Parse PDB files to extract:
 - Secondary structure annotations (HELIX/SHEET)
 - SEQRES length
 - Missing-residue gaps and low-confidence coils
 - Small-molecule (HETATM) lists
"""

import re
import warnings
from typing import List, Tuple


from Bio import BiopythonDeprecationWarning

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


def parse_helix_sheet_annotations(pdb_file: str, chain_id: str = "A") -> dict:
    """Extract HELIX and SHEET annotations for a given chain.

    Args:
        pdb_file (str): Path to the PDB file.
        chain_id (str): Chain identifier to inspect, defaults to ``'A'``.

    Returns:
        dict: Mapping of residue numbers to ``'H'`` or ``'E'`` designations."""
    ss_map = {}
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("HELIX") and len(line) >= 38 and line[19] == chain_id:
                start = int(line[21:25].strip())
                end = int(line[33:37].strip())
                for res_id in range(start, end + 1):
                    ss_map[res_id] = "H"
            elif line.startswith("SHEET") and len(line) >= 37 and line[21] == chain_id:
                start = int(line[22:26].strip())
                end = int(line[33:37].strip())
                for res_id in range(start, end + 1):
                    ss_map[res_id] = "E"
    return ss_map


def get_seqres_length(pdb_file: str, chain_id: str = "A") -> int:
    """Count the number of residues declared in SEQRES records for a chain.

    Args:
        pdb_file (str): Path to the PDB file.
        chain_id (str): Chain identifier to target, defaults to ``'A'``.

    Returns:
        int: Number of residues reported in SEQRES for the chain."""
    length = 0
    with open(pdb_file, "r") as f:
        for ln in f:
            if ln.startswith("SEQRES") and ln[11] == chain_id:
                # count residue names on that line
                length += len(ln[19:].split())
    return length


def group_segments(numbers: List[int]) -> List[Tuple[int, int]]:
    """Collapse a sorted list of residue indices into contiguous ranges.

    Args:
        numbers (List[int]): Sorted sequence of residue numbers.

    Returns:
        List[Tuple[int, int]]: Inclusive ``(start, end)`` segments."""
    from itertools import count, groupby

    segments = []
    for _, grp in groupby(numbers, key=lambda n, c=count(): n - next(c)):
        grp = list(grp)
        segments.append((grp[0], grp[-1]))
    return segments


def detect_gaps_and_discard_coils(
    pdb_file: str, chain_id: str, plddt_threshold: float, coil_length_threshold: int
) -> List[Tuple[str, int, int]]:
    """Identify missing residues and low-confidence coil regions in a structure.

    Args:
        pdb_file (str): Path to the PDB file to analyse.
        chain_id (str): Chain identifier to evaluate.
        plddt_threshold (float): Minimum confidence (B-factor) to keep a residue.
        coil_length_threshold (int): Minimum consecutive residues to treat as a coil.

    Returns:
        List[Tuple[str, int, int]]: Issue records in the form ``(type, start, end)``."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("struct", pdb_file)
    except Exception:
        return []

    # find chain
    chain = next((c for c in structure[0] if c.id == chain_id), None)
    if chain is None:
        return []

    present_ids = []
    seq_plddt = []

    # record which residues are present and their avg B-factor
    for res in chain:
        hetflag, resnum, icode = res.id
        if hetflag == " " and is_aa(res):
            present_ids.append(resnum)
            bvals = [a.bfactor for a in res if a.bfactor is not None]
            avg = sum(bvals) / len(bvals) if bvals else 0.0
            seq_plddt.append((resnum, avg))

    if not present_ids:
        return []

    min_pres, max_pres = min(present_ids), max(present_ids)
    missing = []

    # parse REMARK 465 lines
    pattern = re.compile(rf"^REMARK 465\s+\w{{3}}\s+{chain_id}\s+(\d+)")
    with open(pdb_file, "r") as fh:
        for ln in fh:
            m = pattern.match(ln)
            if m:
                missing.append(int(m.group(1)))

    issues = []
    if missing:
        for start, end in group_segments(sorted(missing)):
            if end < min_pres:
                t = "GAP_NTER"
            elif start > max_pres:
                t = "GAP_CTER"
            else:
                t = "GAP"
            issues.append((t, start, end))
    else:
        # infer gaps by numbering discontinuities
        sorted_ids = sorted(present_ids)
        for a, b in zip(sorted_ids, sorted_ids[1:]):
            if b != a + 1:
                issues.append(("GAP", a + 1, b - 1))
        # Nâterminal
        if min_pres > 1:
            issues.append(("GAP_NTER", 1, min_pres - 1))
        # Câterminal
        seq_len = get_seqres_length(pdb_file, chain_id) or max_pres
        if max_pres < seq_len:
            issues.append(("GAP_CTER", max_pres + 1, seq_len))

    # now detect coil (lowâconfidence) segments
    coil_regions = []
    current = []
    for resnum, avg in seq_plddt:
        if avg < plddt_threshold:
            current.append(resnum)
        else:
            if len(current) >= coil_length_threshold:
                coil_regions.append(tuple(current))
            current = []
    if len(current) >= coil_length_threshold:
        coil_regions.append(tuple(current))

    for region in coil_regions:
        issues.append(("Coil_rejetÃ©", region[0], region[-1]))

    return issues


def parse_small_molecules(pdb_file: str) -> List[str]:
    """List unique small-molecule residue names present in a PDB file.

    Args:
        pdb_file (str): Path to the PDB file to inspect.

    Returns:
        List[str]: Sorted small-molecule residue identifiers, excluding water."""
    ligands = set()
    with open(pdb_file, "r") as f:
        for ln in f:
            # header HET records
            if ln.startswith("HET   "):
                resn = ln[7:10].strip()
                if resn and resn not in ("HOH", "WAT"):
                    ligands.add(resn)
            # actual coordinates
            elif ln.startswith("HETATM"):
                resn = ln[17:20].strip()
                if resn and resn not in ("HOH", "WAT"):
                    ligands.add(resn)
    return sorted(ligands)
