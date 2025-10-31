from __future__ import annotations

# ---------------------------------------------------------------------------
# Imports & logging
# ---------------------------------------------------------------------------
import logging
import tempfile
import warnings
from uuid import uuid4


from Bio import BiopythonDeprecationWarning

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import MDAnalysis as mda
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, PDBParser, Select, Superimposer
from Bio.PDB.Polypeptide import PPBuilder
from MDAnalysis.analysis import align, rms
from MDAnalysis.exceptions import SelectionWarning

from structure_analysis import group_segments

# Suppress specific warnings commonly encountered with MDAnalysis
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", message="Reader has no dt information")
warnings.filterwarnings("ignore", message="Unit cell dimensions not found*")
warnings.filterwarnings(
    "ignore", message="Found no information for attr: 'formalcharges'"
)
warnings.filterwarnings("ignore", message="Unknown element .*")
warnings.filterwarnings("ignore", category=SelectionWarning)

# Reduce MDAnalysis logger level if needed
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

log = logging.getLogger(__name__)

# Constants for DUM detection
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


def _get_max_resid(src: str, chain_id: str, parser: Optional[PDBParser] = None) -> int:
    """

    Parameters
    ----------
    src : str
        Path to the PDB file
    chain_id : str
        Chain identifier
    parser : PDBParser, optional
        PDB parser to reuse, default None (creates a new parser)
    Returns
    -------
    int
        Maximum residue number
    Raises
    ------
    RuntimeError
        If the chain is not found
    """
    if parser is None:
        parser = PDBParser(QUIET=True)
    struct = parser.get_structure("tmp", src)
    for ch in struct.get_chains():
        if ch.id == chain_id:
            return max(res.id[1] for res in ch if res.id[0] == " " and "CA" in res)
    raise RuntimeError(f"Chain {chain_id} not found in {src!r}")


# ---------------------------------------------------------------------------
# Chain-pairing utilities
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------


def _find_best_chain_pair(pred_dir: str, exp_path: str) -> Dict[str, Tuple[str, Dict]]:
    """Map each predicted chain to its closest experimental chain using RMSD."""
    pred_root = Path(pred_dir)
    exp_path = Path(exp_path)

    parser = PDBParser(QUIET=True)
    pp_builder = PPBuilder()  # For extracting sequences
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = 0.0
    aligner.extend_gap_score = 0.0

    exp_struct = parser.get_structure("exp", str(exp_path))
    exp_chains: Dict[str, Dict[int, np.ndarray]] = {}
    exp_sequences: Dict[str, str] = {}

    for ch in exp_struct.get_chains():
        exp_chains[ch.id] = {
            res.id[1]: res["CA"].get_coord()
            for res in ch
            if res.id[0] == " " and "CA" in res
        }
        seq = "".join(str(pp.get_sequence())
                      for pp in pp_builder.build_peptides(ch))
        if seq:
            exp_sequences[ch.id] = seq

    results: Dict[str, Tuple[str, Dict]] = {}
    used_exp_chains: set[str] = set()

    for pred_path in pred_root.glob("*.pdb"):
        pred_name = pred_path.name
        pred_struct = parser.get_structure(pred_name, str(pred_path))

        for p_chain in pred_struct.get_chains():
            p_ca = {
                res.id[1]: res["CA"].get_coord()
                for res in p_chain
                if res.id[0] == " " and "CA" in res
            }
            pred_seq = "".join(
                str(pp.get_sequence()) for pp in pp_builder.build_peptides(p_chain)
            )
            if not pred_seq:
                continue

            best_rmsd = float("inf")
            best_exp_chain = None
            best_score = float("-inf")
            pairs: List[Dict[str, Any]] = []

            for exp_id, exp_seq in exp_sequences.items():
                if exp_id in used_exp_chains:
                    continue

                score = aligner.score(pred_seq, exp_seq) if exp_seq else 0.0

                if score > best_score and exp_id in exp_chains and p_ca:
                    shared = sorted(p_ca.keys() & exp_chains[exp_id].keys())
                    if shared:
                        coords_pred = np.array([p_ca[i] for i in shared])
                        coords_exp = np.array(
                            [exp_chains[exp_id][i] for i in shared])
                        rmsd_val = rms.rmsd(
                            coords_pred, coords_exp, superposition=True)
                        pairs.append(
                            {
                                "pred_chain": p_chain.id,
                                "exp_chain": exp_id,
                                "rmsd": float(rmsd_val),
                                "seq_score": float(score),
                            }
                        )
                        best_score = score
                        if rmsd_val < best_rmsd:
                            best_rmsd = rmsd_val
                            best_exp_chain = exp_id

            if best_exp_chain is None:
                for exp_id, exp_ca in exp_chains.items():
                    if exp_id in used_exp_chains:
                        continue
                    shared = sorted(p_ca.keys() & exp_ca.keys())
                    if not shared:
                        continue
                    coords_pred = np.array([p_ca[i] for i in shared])
                    coords_exp = np.array([exp_ca[i] for i in shared])
                    rmsd_val = rms.rmsd(
                        coords_pred, coords_exp, superposition=True)
                    pairs.append(
                        {
                            "pred_chain": p_chain.id,
                            "exp_chain": exp_id,
                            "rmsd": float(rmsd_val),
                        }
                    )
                    if rmsd_val < best_rmsd:
                        best_rmsd = rmsd_val
                        best_exp_chain = exp_id

            if best_exp_chain is not None:
                used_exp_chains.add(best_exp_chain)
                key = f"{pred_name}:{p_chain.id}"
                results[key] = (
                    best_exp_chain,
                    {
                        "rmsd": best_rmsd,
                        "pairs": pairs,
                        "source": str(pred_path),
                        "pred_chain": p_chain.id,
                    },
                )

    return results


# ---------------------------------------------------------------------------
# Strict PDB reorder utility
# ---------------------------------------------------------------------------
def reorder_residues_in_pdb(pdb_input: str, pdb_output: str) -> None:
    """Rewrite a PDB file with canonical residue and atom ordering.

    Args:
        pdb_input (str): Input PDB path.
        pdb_output (str): Destination path for the reordered structure."""
    from collections import defaultdict

    with open(pdb_input, "r") as f:
        lines = f.readlines()

    atom_lines = [
        line for line in lines if line.startswith(("ATOM", "HETATM"))]
    other_lines = [
        line for line in lines if not line.startswith(("ATOM", "HETATM", "TER"))
    ]

    # Group atoms per (chain, resid, insertion_code, resname)
    residues = defaultdict(list)
    for line in atom_lines:
        chain_id = line[21]
        res_seq = int(line[22:26])
        i_code = line[26] if line[26] != " " else ""
        resname = line[17:20]
        key = (chain_id, res_seq, i_code, resname)
        residues[key].append(line)

    # Sort keys by chain, resid, insertion code
    sorted_keys = sorted(residues.keys(), key=lambda x: (x[0], x[1], x[2]))

    new_lines = []
    serial = 1
    prev_chain = None
    for key in sorted_keys:
        chain_id = key[0]
        if prev_chain and chain_id != prev_chain:
            new_lines.append("TER\n")
        for line in sorted(residues[key], key=lambda entry: entry[12:16]):
            new_line = f"{line[:6]}{serial:5d}{line[11:]}"[:80] + "\n"
            new_lines.append(new_line)
            serial += 1
        prev_chain = chain_id

    new_lines.append("TER\n")
    new_lines += [line for line in other_lines if not line.strip().startswith("END")]
    new_lines.append("END\n")

    with open(pdb_output, "w") as f:
        f.writelines(new_lines)


# ---------------------------------------------------------------------------
# PDB manipulation helpers
# ---------------------------------------------------------------------------


class _ChainSelector(Select):
    """Select only *wanted* chain up to residue *maxr* (inclusive)."""

    def __init__(self, wanted: str, maxr: int):
        super().__init__()
        self.wanted = wanted
        self.max = maxr

    def accept_residue(self, res):  # noqa: D401
        return (
            res.get_parent().id == self.wanted
            and res.id[0] == " "
            and res.id[1] <= self.max
        )


def _clean_chain(
    src: str, chain: str, maxr: int, tmpdir: str, parser: Optional[PDBParser] = None
) -> str:
    """Write *chain* from *src* into a temporary PDB and return its path."""
    if parser is None:
        parser = PDBParser(QUIET=True)

    tmp = Path(tmpdir) / f"clean_{chain}_{uuid4().hex}.pdb"
    io = PDBIO()
    io.set_structure(parser.get_structure("s", src))
    io.save(str(tmp), _ChainSelector(chain, maxr))
    return str(tmp)


def _pad_pdb(src: str, tmpdir: str) -> str:
    """Pad lines to 80 chars, remove ANISOU/TER/END, add TER between chains."""
    tmp = Path(tmpdir) / f"pad_{uuid4().hex}.pdb"
    prev_chain = None
    with open(src) as fin, open(tmp, "w") as fout:
        for ln in fin:
            if ln.startswith(("ATOM  ", "HETATM")):
                chain = ln[21]
                if prev_chain and chain != prev_chain:
                    fout.write("TER\n")
                prev_chain = chain
                fout.write((ln.rstrip("\n") + " " * 80)[:80] + "\n")
            elif ln.startswith(("ANISOU", "TER", "END")):
                continue
            else:
                fout.write(ln)
    return str(tmp)


def _trim_chain_to_universe(
    src: str, chain_id: str, tmpdir: str, parser: Optional[PDBParser] = None
) -> mda.Universe:
    """Return an *MDAnalysis Universe* containing only `chain_id`."""
    if parser is None:
        parser = PDBParser(QUIET=True)

    tmp = _clean_chain(
        src, chain_id, _get_max_resid(src, chain_id, parser), tmpdir, parser
    )
    padded = _pad_pdb(tmp, tmpdir)
    return mda.Universe(padded)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


def _safe_rmsd(u1: mda.Universe, u2: mda.Universe, sel: str) -> Tuple[float, int]:
    def _pick(u):
        ag = u.select_atoms(sel + ' and altloc ""')
        if not ag.n_atoms:
            ag = u.select_atoms(sel + ' and altloc "A"')
        return ag if ag.n_atoms else u.select_atoms(sel)

    mob, ref = _pick(u1), _pick(u2)
    mob_map = {a.resid: a.position for a in mob if a.name == "CA"}
    ref_map = {a.resid: a.position for a in ref if a.name == "CA"}
    common = sorted(mob_map.keys() & ref_map.keys())
    if not common:
        raise RuntimeError("No common residues for RMSD")
    mob_pos = np.vstack([mob_map[r] for r in common])
    ref_pos = np.vstack([ref_map[r] for r in common])
    return float(rms.rmsd(mob_pos, ref_pos, superposition=False)), len(common)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


def _reassign_chain(seg: mda.AtomGroup, chain_id: str) -> None:
    """Set both chainID (atom) and segid (segment) to *chain_id*.

    Setting just one of the two is unreliable across MDAnalysis versions.
    """
    # atoms: chainID attribute exists in topology
    for atom in seg.atoms:
        atom.chainID = chain_id
    # segments: proper place for segid
    for segment in seg.segments:
        segment.segid = chain_id


# ------ Adding new functions for handling fusion proteins and patching missing residues
# ---------------------------------------------------------------------------
# Enhanced patching helpers - from Bio.PDB
# ---------------------------------------------------------------------------


def _is_fusion_protein(pred_chain, exp_chain) -> bool:
    """
    Determine if the experimental chain is a fusion protein (longer than predicted).

    Parameters
    ----------
    pred_chain : Bio.PDB.Chain
        Predicted chain
    exp_chain : Bio.PDB.Chain
        Experimental chain
    Returns
    -------
    bool
        True if the experimental chain is a fusion protein
    """
    if not pred_chain or not exp_chain:
        return False
    # Get maximum residue ID for each chain
    try:
        pred_max = max(
            res.id[1] for res in pred_chain.get_residues() if res.has_id("CA")
        )
        exp_max = max(res.id[1]
                      for res in exp_chain.get_residues() if res.has_id("CA"))
        return exp_max > pred_max
    except ValueError:  # Empty chain
        return False


def _truncate_fusion_protein(exp_chain, max_resid) -> None:
    """
    Truncate experimental chain to match the length of the predicted chain.

    Parameters
    ----------
    exp_chain : Bio.PDB.Chain
        Experimental chain to truncate
    max_resid : int
        Maximum residue number to keep
    """
    # Remove residues beyond the maximum residue ID of the predicted chain
    for res in list(exp_chain.get_residues()):
        if res.id[1] > max_resid:
            exp_chain.detach_child(res.id)


# ---------------------------------------------------------------------------
# Alignment with reference structure
# ---------------------------------------------------------------------------


def align_with_reference(
    structure_path: str, reference_path: str, output_path: str
) -> Tuple[float, float]:
    """Superimpose a structure onto a reference using alpha-carbon atoms.

    Args:
        structure_path (str): Path to the structure being aligned.
        reference_path (str): Path to the reference structure.
        output_path (str): Destination path for the aligned structure.

    Returns:
        Tuple[float, float]: RMSD values before and after alignment.

    Raises:
        ValueError: If no CA atoms are found in either structure."""
    parser = PDBParser(QUIET=True)

    # Load structures
    ref_struct = parser.get_structure("ref", reference_path)
    tgt_struct = parser.get_structure("target", structure_path)

    # Get CA atoms
    ref_atoms = [a for a in ref_struct[0].get_atoms() if a.get_id() == "CA"]
    tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == "CA"]

    # Use the minimum number of atoms
    n = min(len(ref_atoms), len(tgt_atoms))

    if n == 0:
        raise ValueError("No CA atoms found in one or both structures")

    # Calculate RMSD before alignment
    ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
    tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
    rms_before = np.sqrt(((ref_coords - tgt_coords) ** 2).sum(axis=1).mean())

    # Superimpose
    sup = Superimposer()
    sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
    rms_after = sup.rms

    # Apply transformation
    sup.apply(tgt_struct.get_atoms())

    # Save aligned structure
    io = PDBIO()
    io.set_structure(tgt_struct)
    io.save(output_path)

    return rms_before, rms_after


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


def align_and_patch(
    af_path: str,
    exp_path: str,
    gene: str,
    prot_len: int,  # Kept for API compatibility
    window: int,  # Kept for API compatibility
    outdir: str,
    partner_ids: Optional[List[str]] = None,
    uniprot_id: Optional[str] = None,
    reference_paths: Optional[Dict[str, str]] = None,
) -> dict:
    """Align prediction files onto an experimental structure and patch missing regions.

    Args:
        af_path (str): Directory containing predicted PDB files.
        exp_path (str): Experimental PDB file path.
        gene (str): Gene name associated with the structure.
        prot_len (int): Protein length retained for backwards compatibility.
        window (int): Sliding window size retained for backwards compatibility.
        outdir (str): Output directory for generated files.
        partner_ids (Optional[List[str]]): Optional UniProt partner identifiers.
        uniprot_id (Optional[str]): Primary UniProt accession.
        reference_paths (Optional[Dict[str, str]]):
            Precomputed references by UniProt ID.

    Returns:
        dict: Metrics and file paths describing the alignment process.

    Raises:
        RuntimeError: If no predictions are found or no chain pairing is possible."""
    af_dir = Path(af_path)
    exp_file = Path(exp_path)
    out_dir_path = Path(outdir)

    all_universes = []

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # ---------------------------------------------------------------------
            # Discover prediction files
            # ---------------------------------------------------------------------
            if not any(af_dir.glob("*.pdb")):
                raise RuntimeError(f"No PDB files found in {af_dir!r}")

            mapping = _find_best_chain_pair(af_path, exp_path)
            valid = [(k, v) for k, v in mapping.items() if v[0] is not None]
            if not valid:
                raise RuntimeError(
                    "No chain could be paired with experimental model")

            pdbid_raw = exp_file.stem.split("_")[-1].upper()
            pdb_url = f"https://www.rcsb.org/structure/{pdbid_raw}"
            pdbid = f'=HYPERLINK("{pdb_url}", "{pdbid_raw}")'
            res: Dict = {
                "gene": gene,
                "uniprot_id": uniprot_id,
                "pdb_id": pdbid,
                "pdb_id_raw": pdbid_raw,
                "aligned_pdb": "",
                "patched_pdb": "",
                "patched_with_hetatm_pdb": "",
                "oriented_pdb": "",  # New field for the oriented structure
                "orientation_method": None,
                "fusion_protein": False,
                "is_heteromer": False,
                "interaction_partners": "",
                "n_CA_common_before_patch": 0,
                "n_CA_common_after_patch": 0,
                "gaps_found": "",
                "gaps_used": "",
                "gap_count": 0,
                "nb_atoms_to_patch": 0,
                "patched_atoms": 0,
                "%_gaps_patched": 0.0,
                "rmsd_exp": None,
                "rmsd_patched": None,
                "rmsd_oriented": None,  # New field for orientation RMSD
                "protonated": False,
                "rmsd_pairs": {},
                "status": "error",
            }

            # ------ Check if experimental structure is heteromeric
            # Analyze the experimental structure
            exp_parser = PDBParser(QUIET=True)
            exp_structure = exp_parser.get_structure("exp", str(exp_file))
            exp_chains = list(exp_structure[0].get_chains())
            exp_chain_ids = [ch.id for ch in exp_chains]

            log.info(
                f"Experimental structure: {pdbid} with chains {exp_chain_ids}")
            log.info(
                f"Partner IDs: {partner_ids} (type: {type(partner_ids).__name__})")
            if uniprot_id:
                log.info(f"UniProt ID: {uniprot_id}")

            # Normalize partner_ids format (regardless of chain count)
            if partner_ids is None:
                partner_ids = []
            elif isinstance(partner_ids, str):
                partner_ids = [partner_ids]
            # Process partner_ids and set metadata
            if len(exp_chains) > 1:
                # Log chain correspondences
                log.info("Chain-protein correspondences:")
                for key, (exp_chain, details) in mapping.items():
                    pred_file = Path(details["source"]).name
                    pred_chain = details["pred_chain"]
                    log.info(
                        "Mapping %s:%s -> %s:%s (RMSD: %s)",
                        pred_file,
                        pred_chain,
                        pdbid_raw,
                        exp_chain,
                        details.get("rmsd", "N/A"),
                    )

                # Associate partners with chains if possible
                if partner_ids and len(partner_ids) <= len(exp_chain_ids):
                    for i, (partner, chain) in enumerate(
                        zip(partner_ids, exp_chain_ids)
                    ):
                        log.info("Partner %s -> Chain %s", partner, chain)
            # Set metadata
            res["interaction_partners"] = partner_ids
            res["is_heteromer"] = bool(partner_ids)
            hetero_state = "HETEROMER" if res["is_heteromer"] else "HOMOMER"
            log.info("Structure defined as: %s", hetero_state)
            # ---------------------------------------------------------------------
            # ---------------------------------------------------------------------
            # Create a temporary padded version of the experimental structure
            padded_exp = _pad_pdb(str(exp_file), tmpdir)

            hybrid = mda.Universe(padded_exp)
            all_universes.append(hybrid)

            predicted: List[Tuple[str, str, mda.Universe]] = []

            for key, (exp_chain, details) in valid:
                src_pred = details["source"]
                pred_chain = details["pred_chain"]

                # ----- Check for fusion protein case
                pred_struct = exp_parser.get_structure("pred", src_pred)
                for p_chain in pred_struct.get_chains():
                    if p_chain.id == pred_chain:
                        for e_chain in exp_chains:
                            if e_chain.id == exp_chain:
                                if _is_fusion_protein(p_chain, e_chain):
                                    res["fusion_protein"] = True
                                    pred_max = max(
                                        res.id[1]
                                        for res in p_chain.get_residues()
                                        if res.has_id("CA")
                                    )
                                    _truncate_fusion_protein(e_chain, pred_max)
                                break
                        break

                uA = _trim_chain_to_universe(
                    src_pred, pred_chain, tmpdir, exp_parser)
                uE = _trim_chain_to_universe(
                    str(exp_file), exp_chain, tmpdir, exp_parser
                )

                all_universes.extend([uA, uE])

                # Calculate common residues once and reuse
                ca_resids_A = set(uA.select_atoms("name CA").resids)
                ca_resids_E = set(uE.select_atoms("name CA").resids)
                common = sorted(ca_resids_A & ca_resids_E)

                if not common:
                    continue

                sel = "protein and name CA and resid " + \
                    " ".join(map(str, common))

                try:
                    align.AlignTraj(
                        uA, uE, select=sel, in_memory=True, weights="mass"
                    ).run()
                    rms_val, n_common = _safe_rmsd(uA, uE, sel)
                    rms_val = round(rms_val, 3)
                    res["rmsd_pairs"][exp_chain] = rms_val
                    predicted.append((pred_chain, exp_chain, uA))
                    res["n_CA_common_before_patch"] += len(common)
                except (ValueError, RuntimeError) as e:
                    log.warning(
                        f"Failed alignment for chain {pred_chain}-{exp_chain}: {e}"
                    )
                    continue
            if predicted:
                gene_dir = out_dir_path / gene
                gene_dir.mkdir(parents=True, exist_ok=True)
                aligned_path = gene_dir / f"aligned_on_{exp_file.name}"

                ag_pred = [uA.select_atoms("protein")
                           for _, _, uA in predicted]

                merged = mda.Merge(*ag_pred)
                merged.atoms.write(str(aligned_path))

                res["aligned_pdb"] = str(aligned_path)

            # ---------------------------------------------------------------------
            # Detect and patch gaps per chain
            # ---------------------------------------------------------------------
            for pred_chain, exp_chain, uA in predicted:
                log.info(
                    "Checking gaps for predicted chain %s vs experimental chain %s",
                    pred_chain,
                    exp_chain,
                )

                ids_A = set(uA.select_atoms("name CA").resids)
                ids_H = set(
                    hybrid.select_atoms(
                        f"chainID {exp_chain} and name CA").resids
                )

                extras = sorted(ids_A - ids_H)
                gaps = group_segments(extras) if extras else []

                res["gap_count"] += len(gaps)
                if gaps:
                    res["gaps_found"] += (
                        ";".join(f"{exp_chain}:{s}-{e}" for s, e in gaps) + ";"
                    )

                patched_ranges = []  # For building the dynamic list per pair

                # Pre-calculate total atoms to patch
                total_atoms_to_patch = 0
                for start, end in gaps:
                    seg = uA.select_atoms(f"protein and resid {start}:{end}")
                    if seg.n_atoms:
                        total_atoms_to_patch += seg.n_atoms

                res["nb_atoms_to_patch"] += total_atoms_to_patch

                for start, end in gaps:
                    seg = uA.select_atoms(f"protein and resid {start}:{end}")
                    if not seg.n_atoms:
                        continue
                    _reassign_chain(seg, exp_chain)
                    hybrid = mda.Merge(
                        hybrid.select_atoms("protein"), seg.atoms)
                    res["patched_atoms"] += seg.n_atoms
                    patched_ranges.append(f"{exp_chain}:{start}-{end}")

                if patched_ranges:
                    res["gaps_used"] += ";".join(patched_ranges) + ";"

            # ---------------------------------------------------------------------
            # Write patched structure
            # ---------------------------------------------------------------------
            patched_path = gene_dir / f"patched_{exp_file.name}"
            # CORRECTION: For homomers, keep only one chain
            if not res["is_heteromer"] and predicted:
                # For homomers, keep only the first experimental chain
                # Get the first experimental chain
                first_exp_chain = predicted[0][1]
                log.info(
                    f"Homomer detected: keeping only chain {first_exp_chain}")
                # Select only atoms from the first chain
                hybrid.select_atoms(f"protein and chainID {first_exp_chain}").write(
                    str(patched_path)
                )
            else:
                # Original behavior for heteromers
                hybrid.select_atoms("protein").write(str(patched_path))

            # Reorder and rewrite the same file directly
            reorder_residues_in_pdb(str(patched_path), str(patched_path))
            res["patched_pdb"] = str(patched_path)

            # ---------------------------------------------------------------------
            # Unified RMSD calculation (experimental and patched)
            # ---------------------------------------------------------------------
            rmsd_exp_vals = []
            rmsd_patch_vals = []

            for pred_chain, exp_chain, uA in predicted:
                try:
                    # Experimental comparison
                    uE = _trim_chain_to_universe(
                        exp_path, exp_chain, tmpdir, exp_parser
                    )
                    all_universes.append(uE)

                    sel_common = "name CA"
                    mob = uA.select_atoms(sel_common)
                    ref = uE.select_atoms(sel_common)
                    mob_map = {a.resid: a.position for a in mob}
                    ref_map = {a.resid: a.position for a in ref}
                    common = sorted(mob_map.keys() & ref_map.keys())
                    if common:
                        rmsd = float(
                            rms.rmsd(
                                np.vstack([mob_map[r] for r in common]),
                                np.vstack([ref_map[r] for r in common]),
                                superposition=False,
                            )
                        )
                        rmsd_exp_vals.append((rmsd, len(common)))

                    # Patched comparison
                    ref = hybrid.select_atoms(
                        f"chainID {exp_chain} and name CA")
                    ref_map = {a.resid: a.position for a in ref}
                    common = sorted(mob_map.keys() & ref_map.keys())
                    if common:
                        rmsd = float(
                            rms.rmsd(
                                np.vstack([mob_map[r] for r in common]),
                                np.vstack([ref_map[r] for r in common]),
                                superposition=False,
                            )
                        )
                        rmsd_patch_vals.append((rmsd, len(common)))
                        res["n_CA_common_after_patch"] += len(common)
                except (ValueError, RuntimeError) as e:
                    log.warning(
                        "Failed RMSD computation for chain %s-%s: %s",
                        pred_chain,
                        exp_chain,
                        e,
                    )
                except Exception as e:
                    log.error(f"Unexpected error during RMSD computation: {e}")
                    # Don't re-raise to maintain compatibility

            # Final RMSD aggregation
            if rmsd_exp_vals:
                total = sum(n for _, n in rmsd_exp_vals)
                res["rmsd_exp"] = round(
                    sum(r * n for r, n in rmsd_exp_vals) / total, 3)

            if rmsd_patch_vals:
                total = sum(n for _, n in rmsd_patch_vals)
                res["rmsd_patched"] = round(
                    sum(r * n for r, n in rmsd_patch_vals) / total, 3
                )

            # ---------------------------------------------------------------------
            # HETATM reintegration (optional)
            # ---------------------------------------------------------------------
            with exp_file.open() as f:
                het_atoms = [ln for ln in f if ln.startswith("HETATM")]

            if het_atoms:
                hetout = gene_dir / f"patched_with_hetatm_{exp_file.name}"
                with patched_path.open() as p_in, hetout.open("w") as p_out:
                    p_out.writelines(
                        [ln for ln in p_in if not ln.startswith("END")]
                        + het_atoms
                        + ["END\n"]
                    )
                res["patched_with_hetatm_pdb"] = str(hetout)
            # ---------------------------------------------------------------------
            # Final metrics
            # ---------------------------------------------------------------------
            if res["gap_count"]:
                done = (
                    len(res["gaps_used"].rstrip(";").split(";"))
                    if res["gaps_used"]
                    else 0
                )
                res["%_gaps_patched"] = round(100 * done / res["gap_count"], 1)
            res["status"] = "ok"

            return res

        finally:
            for universe in all_universes:
                try:
                    close = getattr(universe, "close", None)
                    if callable(close):
                        close()
                    else:
                        trajectory = getattr(universe, "trajectory", None)
                        if trajectory is not None and callable(
                            getattr(trajectory, "close", None)
                        ):
                            trajectory.close()
                except Exception:
                    pass
