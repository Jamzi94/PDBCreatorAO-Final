#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# Module d'alignement multi-r????f????rence pour structures PDB
# ---------------------------------------------------------------------------

import logging
import os
import re
import subprocess
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

from config import PipelineConfig, SERVICE_CONFIG, ServiceConfig


from Bio import BiopythonDeprecationWarning

from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple
from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Superimposer
from networking import get_session
from MDAnalysis.exceptions import SelectionWarning


# --- Constantes & expressions régulières ---
DUM_RE = re.compile(r"\bDUM\b", re.IGNORECASE)

@dataclass
class OrientationMetadata:
    """Métadonnées d'orientation d'une structure PDB."""
    source: str  # 'OPM', 'PDBTM', 'RCSB'
    reference_id: Optional[str] = None  # ID PDB de la structure de référence utilisée
    ref_uniprot: Optional[str] = None  # UniProt ID de la référence
    convention: str = "+Z=extra"  # Convention d'orientation
    transform_matrix: Optional[np.ndarray] = None  # Matrice 4x4 de transformation
    rmsd_before: Optional[float] = None  # RMSD avant alignement
    rmsd_after: Optional[float] = None  # RMSD après alignement
    z_flipped: bool = False  # Si un flip en Z a été appliqué
    
    def to_dict(self) -> Dict[str, Any]:
        """Convertit les métadonnées en dictionnaire pour stockage."""
        d = {
            "source": self.source,
            "convention": self.convention,
            "z_flipped": self.z_flipped
        }
        if self.reference_id:
            d["reference_id"] = self.reference_id
            d["ref_uniprot"] = self.ref_uniprot
        if self.transform_matrix is not None:
            d["transform_matrix"] = self.transform_matrix.tolist()
        if self.rmsd_before is not None:
            d["rmsd_before"] = self.rmsd_before
            d["rmsd_after"] = self.rmsd_after
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "OrientationMetadata":
        """Crée une instance à partir d'un dictionnaire."""
        matrix = d.get("transform_matrix")
        if matrix:
            matrix = np.array(matrix)
        return cls(
            source=d["source"],
            reference_id=d.get("reference_id"),
            ref_uniprot=d.get("ref_uniprot"),
            convention=d.get("convention", "+Z=extra"),
            transform_matrix=matrix,
            rmsd_before=d.get("rmsd_before"),
            rmsd_after=d.get("rmsd_after"),
            z_flipped=d.get("z_flipped", False)
        )

# --- Configuration du logging ---
logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=SelectionWarning)
# ---------------------------------------------------------------------------
# Fonctions utilitaires pour l'extraction et le traitement des données
# ---------------------------------------------------------------------------


@lru_cache(maxsize=128)
def _load_structure_cached(path_str: str):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(Path(path_str).stem, path_str)


def _load_structure(path: Path):
    return deepcopy(_load_structure_cached(str(path)))


def _superimpose_worker(
    args: tuple[int, str, str, str, str, str, str],
) -> tuple[int, Optional[tuple[str, str, float, float]]]:
    from pipeline import _configure_worker_logging
    _configure_worker_logging()
    warnings.filterwarnings("ignore", category=SelectionWarning)
    warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
    
    
    idx, ref_path_str, out_root_str, pdb_id, gene, uniprot, patched_pdb = args
    try:
        pdb_path = Path(patched_pdb)
        if not pdb_path.is_file():
            logger.warning("Missing target file: %s", pdb_path)
            return idx, None

        ref_struct = _load_structure(Path(ref_path_str))
        tgt_struct = _load_structure(pdb_path)

        ref_atoms = [a for a in ref_struct[0].get_atoms() if a.get_id() == "CA"]
        tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == "CA"]
        n = min(len(ref_atoms), len(tgt_atoms))
        if n == 0:
            logger.warning("No CA atoms in common for %s", pdb_id)
            return idx, None

        ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
        tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
        rms_before = float(np.sqrt(((ref_coords - tgt_coords) ** 2).sum(axis=1).mean()))

        sup = Superimposer()
        sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
        sup.apply(tgt_struct.get_atoms())
        rms_after = float(sup.rms)

        out_root = Path(out_root_str)
        subdir = out_root / f"{gene}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        out_path = subdir / f"oriented_{gene}_{uniprot}_{pdb_id}.pdb"

        io = PDBIO()
        io.set_structure(tgt_struct)
        io.save(str(out_path))

        logger.info(
            "%s: before=%.3f, after=%.3f -> %s", pdb_id, rms_before, rms_after, out_path
        )
        return idx, (pdb_id, str(out_path), rms_before, rms_after)
    except Exception as exc:
        logger.error("Failed to orient %s: %s", pdb_id, exc)
        return idx, None



def _choose_native_reference(
    df_group: pd.DataFrame,
    source_map: Dict[str, str],
    cfg: PipelineConfig
) -> tuple[Optional[str], OrientationMetadata]:
    """Sélectionne la meilleure référence native pour un groupe UniProt.
    
    Args:
        df_group: DataFrame avec les structures d'un même groupe UniProt
        source_map: Dictionnaire {pdb_id: source}
        cfg: Configuration avec les préférences de sources natives
    
    Returns:
        tuple[str | None, OrientationMetadata]: (ID PDB de référence, métadonnées)
    """
    # Vérifie si le DataFrame est vide
    if df_group.empty:
        return None, OrientationMetadata(source="none")
        
    # Ajoute la source dans le DataFrame pour faciliter la sélection
    df_group = df_group.copy()
    df_group["source"] = df_group["pdb_id"].map(source_map)
    
    # Cherche d'abord une référence déjà orientée selon la priorité des sources
    for source in cfg.native_source_priority:
        source_structures = df_group[df_group["source"] == source]
        if not source_structures.empty:
            chosen_row = source_structures.iloc[0]
            metadata = OrientationMetadata(
                source=source,
                reference_id=None,  # C'est la référence elle-même
                ref_uniprot=chosen_row["uniprot_id"],
                convention="+Z=extra"  # Convention par défaut pour OPM/PDBTM
            )
            logger.info(f"Found native reference from {source}: {chosen_row['pdb_id']}")
            return chosen_row["pdb_id"], metadata
        for _, row in df_group.iterrows():
            pid = row["pdb_id"]
            if source_map.get(pid) == source:
                metadata = OrientationMetadata(
                    source=source,
                    reference_id=None,  # C'est la référence elle-même
                    ref_uniprot=row["uniprot_id"],
                    convention="+Z=extra"  # Convention par défaut pour OPM/PDBTM
                )
                return pid, metadata
                
    # Si on n'a pas trouvé de référence native et qu'on en exige une
    if cfg.require_native_reference:
        return None, OrientationMetadata(source="none")
        
    # Sinon on prend la première structure comme référence
    if not df_group.empty:
        pid = df_group.iloc[0]["pdb_id"]
        source = source_map.get(pid, "RCSB")
        metadata = OrientationMetadata(
            source=source,
            reference_id=None,
            ref_uniprot=df_group.iloc[0]["uniprot_id"],
            convention="+Z=extra"  # On appliquera la convention standard
        )
        return pid, metadata
        
    return None, OrientationMetadata(source="none")


def _find_tm_residues(structure) -> List[int]:
    """Identifie les résidus transmembranaires via leur position Z."""
    ca_atoms = [a for a in structure[0].get_atoms() if a.get_id() == "CA"]
    z_coords = [a.get_coord()[2] for a in ca_atoms]
    z_mean = np.mean(z_coords)
    z_std = np.std(z_coords)
    # On considère TM les résidus dans ±15Å autour du centre
    tm_indices = [i for i, z in enumerate(z_coords) if abs(z - z_mean) < 15.0]
    return tm_indices

def _align_to_reference(ref_path: Path, target_path: Path, out_path: Path) -> tuple[float, float, bool, np.ndarray]:
    """Aligne la structure cible sur la référence en utilisant les CA des résidus TM.
    
    Returns:
        tuple[float, float, bool, np.ndarray]: (RMSD avant, RMSD après, flip Z appliqué, matrice)
    """
    ref_struct = _load_structure(ref_path)
    tgt_struct = _load_structure(target_path)

    ref_tm_indices = _find_tm_residues(ref_struct)
    tgt_tm_indices = _find_tm_residues(tgt_struct)
    
    ref_atoms = [a for i, a in enumerate(ref_struct[0].get_atoms()) if a.get_id() == "CA" and i in ref_tm_indices]
    tgt_atoms = [a for i, a in enumerate(tgt_struct[0].get_atoms()) if a.get_id() == "CA" and i in tgt_tm_indices]
    n = min(len(ref_atoms), len(tgt_atoms))
    if n < 10:  # On demande au moins 10 résidus TM en commun
        ref_atoms = [a for a in ref_struct[0].get_atoms() if a.get_id() == "CA"]
        tgt_atoms = [a for a in tgt_struct[0].get_atoms() if a.get_id() == "CA"]
        n = min(len(ref_atoms), len(tgt_atoms))
        logger.warning(f"Pas assez de résidus TM en commun pour {target_path.name}, utilisation de tous les CA")
    if n == 0:
        raise ValueError(f"No CA atoms in common for {target_path.name}")

    ref_coords = np.array([a.get_coord() for a in ref_atoms[:n]])
    tgt_coords = np.array([a.get_coord() for a in tgt_atoms[:n]])
    rms_before = float(np.sqrt(((ref_coords - tgt_coords) ** 2).sum(axis=1).mean()))

    # Premier alignement
    sup = Superimposer()
    sup.set_atoms(ref_atoms[:n], tgt_atoms[:n])
    sup.apply(tgt_struct.get_atoms())
    rms_after = float(sup.rms)
    transform_matrix = sup.rotran

    # Test du flip Z
    flipped_coords = tgt_coords.copy()
    flipped_coords[:, 2] *= -1  # Flip selon Z
    rms_flipped = float(np.sqrt(((ref_coords - flipped_coords) ** 2).sum(axis=1).mean()))
    
    z_flipped = False
    if rms_flipped < rms_after:
        # Le flip Z améliore l'alignement
        logger.info(f"Applying Z-flip to {target_path.name} (RMSD: {rms_after:.3f} -> {rms_flipped:.3f})")
        for atom in tgt_struct.get_atoms():
            coords = atom.get_coord()
            coords[2] *= -1
            atom.set_coord(coords)
        rms_after = rms_flipped
        z_flipped = True
        # Met à jour la matrice pour inclure le flip Z
        flip_matrix = np.eye(4)
        flip_matrix[2, 2] = -1
        transform_matrix = np.dot(flip_matrix, transform_matrix)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(tgt_struct)
    io.save(str(out_path))

    logger.info("Aligned %s onto %s -> %s", target_path.name, ref_path.name, out_path)
    return rms_before, rms_after, z_flipped, transform_matrix


def _orient_with_fallback(
    pdb_id: str, pdb_path: Path, out_dir: Path, cfg
) -> tuple[Optional[Path], str]:
    """Cascade fallback orientation with services enabled/disabled via config."""
    raw_code = (pdb_id or "").strip()
    if not raw_code:
        logging.error("Orientation fallback called without a PDB identifier")
        return None, "Failed: No PDB ID"

    # Helper to get the correct PDB ID case for each service
    def service_code(name: str) -> str:
        conf = SERVICE_CONFIG.get(name)
        # Default to uppercase for local tools, lowercase for remote, unless specified
        if conf and conf.case_sensitive:
            return raw_code
        if name in {"pdbtm", "memprotmd"}:
            return raw_code.lower()
        return raw_code.upper()

    from structure_repair import (
        orient_with_memembed,
        orient_with_memprotmd_cached,
        orient_with_pdbtm_cached,
        orient_with_ppm,
        orient_with_principal_axes,
        orient_with_tmdet,
    )

    # Define all available services and their handlers
    all_services = {
        "pdbtm": lambda dest: orient_with_pdbtm_cached(service_code("pdbtm"), pdb_path, dest),
        "memprotmd": lambda dest: orient_with_memprotmd_cached(service_code("memprotmd"), pdb_path, dest),
        "ppm": lambda dest: orient_with_ppm(pdb_path, dest, cfg.ppm_path),
        "memembed": lambda dest: orient_with_memembed(pdb_path, dest),
        "tmdet": lambda dest: orient_with_tmdet(pdb_path, dest),
        "axes": lambda dest: orient_with_principal_axes(pdb_path, dest),
    }

    # Build the list of enabled strategies based on native source preferences
    enabled_strategies: list[str] = []
    
    # First try native sources in order of preference if enabled
    if cfg.prefer_native_sources:
        native_sources = {"pdbtm"}  # Since we've removed direct PDBTM/OPM fetching
        
        # Start with native sources according to priority
        for source in native_sources:
            conf = SERVICE_CONFIG.get(source)
            if not conf or conf.enabled:
                enabled_strategies.append(source)
                
        # If native sources are required and we found some, stop here
        if cfg.require_native_reference and enabled_strategies:
            logging.info("Using only native sources: %s", enabled_strategies)
        else:
            # Add remaining local strategies
            for name in all_services:
                if name not in native_sources:
                    conf = SERVICE_CONFIG.get(name)
                    if not conf or conf.enabled:
                        if name != "ppm" or cfg.ppm_path:  # Special PPM check
                            enabled_strategies.append(name)
    else:
        # If not preferring native sources, add all enabled strategies
        for name in all_services:
            conf = SERVICE_CONFIG.get(name)
            if not conf or conf.enabled:
                if name != "ppm" or cfg.ppm_path:  # Special PPM check
                    enabled_strategies.append(name)
                    
    # Determine the final order of execution based on config
    final_order: list[str] = []
    if bool(getattr(cfg, "enable_new_orientation_services", False)):
        # New mode: sort by priority field from config
        sorting_list = [
            (SERVICE_CONFIG.get(name, ServiceConfig(priority=99)).priority, name)
            for name in enabled_strategies
        ]
        sorting_list.sort()
        final_order = [name for _, name in sorting_list]
    else:
        # Legacy mode: hardcoded order for pdbtm/memprotmd, others after
        legacy_order = ["memprotmd", "pdbtm"]
        if getattr(cfg, "prefer_pdbtm_over_memprotmd", False):
            legacy_order = ["pdbtm", "memprotmd"]

        for name in legacy_order:
            if name in enabled_strategies:
                final_order.append(name)
        for name in enabled_strategies:
            if name not in legacy_order:
                final_order.append(name)

    # Execute strategies in the determined order
    for service_name in final_order:
        handler = all_services[service_name]
        oriented_path = out_dir / f"{pdb_path.stem}_oriented_{service_name}.pdb"
        try:
            meta = handler(oriented_path)
        except Exception as exc:
            logging.error("%s orientation failed for %s: %s", service_name.upper(), raw_code.upper(), exc)
            continue

        if meta and oriented_path.is_file() and oriented_path.stat().st_size > 0:
            method = (
                meta.get("method", service_name.upper())
                if isinstance(meta, dict)
                else service_name.upper()
            )
            logging.info(
                "Orientation by %s succeeded for %s -> %s",
                method,
                raw_code.upper(),
                oriented_path,
            )
            return oriented_path, method

    return None, "Failed: All methods exhausted"


def superimpose_structures(
    ref_path: Path, df_targets: pd.DataFrame, out_root: Path
) -> list[tuple[str, str, float, float]]:
    """Superimpose patched structures onto a reference PDB.

    Args:
        ref_path (Path): Reference structure used as alignment target.
        df_targets (pd.DataFrame): Candidate structures to orient.
        out_root (Path): Root directory for generated files.

    Returns:
        list[tuple[str, str, float, float]]:
            (pdb_id, path, rms_before, rms_after) tuples.
    """
    tasks = []
    required_cols = {"pdb_id", "gene_name", "uniprot_id", "patched_pdb"}
    if not required_cols.issubset(df_targets.columns):
        missing = required_cols - set(df_targets.columns)
        logger.error(
            "Missing required columns in target DataFrame for superimposition: %s",
            missing,
        )
        return []

    # We pass the original index to the worker to be able to log errors correctly
    for idx, row in df_targets.iterrows():
        args = (
            idx,
            str(ref_path),
            str(out_root),
            row["pdb_id"],
            row["gene_name"],
            row["uniprot_id"],
            row["patched_pdb"],
        )
        tasks.append(args)

    if not tasks:
        return []

    results = []
    max_workers = os.cpu_count() or 1
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_map = {executor.submit(_superimpose_worker, task): task for task in tasks}
        for future in as_completed(future_map):
            try:
                _, res = future.result()
                if res:
                    results.append(res)
            except Exception as exc:
                task_info = future_map.get(future)
                pdb_id = task_info[3] if task_info else "unknown"
                logger.error(
                    "A superimposition worker failed for PDB ID %s: %s", pdb_id, exc
                )

    return results


def group_structures_by_attribute(
    df_det: pd.DataFrame, attribute: str
) -> Dict[str, List[str]]:
    """Return a mapping {attribute_value: [pdb ids]} using pandas groupby.

    Args:
        df_det (pd.DataFrame): PDB details table.
        attribute (str): Column to group by (e.g., 'gene_name').

    Returns:
        Dict[str, List[str]]: A dictionary mapping attribute values to lists of PDB IDs."""
    if attribute not in df_det.columns:
        logger.error(f"Attribut '{attribute}' non trouvé dans le DataFrame")
        return {}
    if "pdb_id" not in df_det.columns:
        logger.error("Colonne 'pdb_id' manquante dans PDB_Details")
        return {}

    subset = df_det.dropna(subset=[attribute, "pdb_id"])
    if subset.empty:
        return {}
    return subset.groupby(attribute)["pdb_id"].apply(list).to_dict()


def find_dum_references(
    df_det: pd.DataFrame, groups: Dict[str, List[str]]
) -> Dict[str, str]:
    """Select a DUM-containing reference structure for each group.

    Args:
        df_det (pd.DataFrame): Detail records containing small-molecule metadata.
        groups (Dict[str, List[str]]): Mapping from group to candidate PDB identifiers.

    Returns:
        Dict[str, str]: Reference PDB identifier per group when available."""
    references = {}

    if "small_molecules" not in df_det.columns:
        logger.warning("Colonne 'small_molecules' non trouvée, impossible de détecter DUM")
        return references

    # Pour chaque groupe, trouve une référence avec DUM
    for group, pdb_ids in groups.items():
        # Filtre le DataFrame pour ce groupe
        group_df = df_det[df_det["pdb_id"].isin(pdb_ids)].copy()

        # Cherche les structures avec DUM
        sm = group_df["small_molecules"].fillna("").str.upper()
        dum_mask = sm.str.contains(DUM_RE)

        if dum_mask.any():
            # Prend la première structure avec DUM comme référence
            ref_id = group_df.loc[dum_mask, "pdb_id"].iloc[0]
            references[group] = ref_id            
            logger.info(f"Groupe {group}: référence avec DUM trouvée: {ref_id}")
        else:
            logger.info(f"Groupe {group}: aucune référence avec DUM trouvée")
            
    return references


def get_oriented_path(
    out_root: Path,
    group: str,
    uniprot: str,
    pdb_id: str,
    metadata: OrientationMetadata
) -> Path:
    """Génère le chemin pour une structure orientée."""
    subdir = out_root / f"{group}_{uniprot}"
    subdir.mkdir(parents=True, exist_ok=True)
    suffix = "_native" if metadata.source in ("OPM", "PDBTM") else ""
    return subdir / f"oriented_{group}_{uniprot}_{pdb_id}{suffix}.pdb"

def align_structures_by_groups(
    excel_file: Optional[Path | str] = None,
    group_by: str = "gene_name",
    out_root: Optional[Path | str] = None,
    data_frames: Optional[Mapping[str, pd.DataFrame]] = None,
    cfg: Optional[PipelineConfig] = None,
    orientation_fallback_fn: Optional[Callable[[str, Path, Path, PipelineConfig], Tuple[Optional[Path], str]]] = None,
) -> tuple[Dict[str, List[tuple]], pd.DataFrame]:
    """Align patched structures within groups using DUM references.

    Args:
        excel_file (Optional[Path | str]):
            Source Excel workbook, if not supplying DataFrames.
        group_by (str): Column used to build groups (``gene_name`` or ``uniprot_id``).
        out_root (Optional[Path | str]): Destination directory for oriented structures.
        data_frames (Optional[Mapping[str, pd.DataFrame]]):
            Pre-loaded tables keyed by sheet name.
        cfg (Optional[PipelineConfig]): Pipeline configuration for feature flags.

    Returns:
        tuple[Dict[str, List[tuple]], pd.DataFrame]:
            Raw results per group and a consolidated DataFrame.

    Raises:
        ValueError: If neither an Excel file nor data frames are provided."""
    if cfg is None:
        cfg = PipelineConfig.from_args()  # Default config if not provided

    # This function is now defined in pipeline.py, but we need a way to load sheets
    # if we are running this module standalone.
    from pipeline import _prepare_dataframe, load_sheet, filter_and_cascade

    if data_frames:
        try:
            df_det_all = _prepare_dataframe(data_frames["PDB_Details"])
            df_ali_all = _prepare_dataframe(data_frames["Alignment_Summary"])
        except KeyError as exc:
            raise ValueError(
                "data_frames must contain 'PDB_Details' and 'Alignment_Summary'"
            ) from exc
        excel_path = Path(excel_file) if excel_file else None
    else:
        if excel_file is None:
            raise ValueError("excel_file or data_frames must be provided")
        excel_path = Path(excel_file)
        df_det_all = _prepare_dataframe(load_sheet(excel_path, "PDB_Details"))
        df_ali_all = _prepare_dataframe(load_sheet(excel_path, "Alignment_Summary"))

    out_root_path = (
        Path(out_root) if out_root is not None else Path("pdb_structures_oriented")
    )
    out_root_path.mkdir(parents=True, exist_ok=True)

    if "pdb_id" in df_det_all.columns:
        df_det_all["pdb_id"] = df_det_all["pdb_id"].astype(str).str.upper()
    if "pdb_id" in df_ali_all.columns:
        df_ali_all["pdb_id"] = df_ali_all["pdb_id"].astype(str).str.upper()

    df_det_all, df_ali_all = filter_and_cascade(df_det_all, df_ali_all)

    required_cols = {"pdb_id", group_by, "small_molecules"}
    missing = required_cols - set(df_det_all.columns)
    if missing:
        logger.error(f"Colonnes requises manquantes dans PDB_Details: {missing}")
        return {}, pd.DataFrame()

    if "patched_pdb" not in df_ali_all.columns:
        logger.error("Colonne 'patched_pdb' non trouvée dans Alignment_Summary")
        return {}, pd.DataFrame()

    groups = group_structures_by_attribute(df_det_all, group_by)
    logger.info(f"Grouped structures by {group_by}: {len(groups)} groups")

    references = find_dum_references(df_det_all, groups)
    logger.info(f"Found DUM references for {len(references)}/{len(groups)} groups")

    source_map: Dict[str, str] = {}
    if "pdb_source" in df_det_all.columns:
        source_map = {
            str(row["pdb_id"]): str(row["pdb_source"]).upper()
            for _, row in df_det_all.iterrows()
            if pd.notna(row.get("pdb_id"))
        }
        
    # Ensure we have required columns with valid data
    missing_data = df_det_all[["pdb_id", "gene_name", "uniprot_id"]].isna().any(axis=1)
    if missing_data.any():
        logger.warning(f"Missing required data for {missing_data.sum()} structures")
        # Fill missing gene_name with uniprot_id if available
        mask = df_det_all["gene_name"].isna() & df_det_all["uniprot_id"].notna()
        df_det_all.loc[mask, "gene_name"] = df_det_all.loc[mask, "uniprot_id"]

    all_results: Dict[str, List[tuple]] = {}
    all_aligned_results: List[tuple] = []
    orientation_metadata: Dict[str, OrientationMetadata] = {}

    # Helper function for consistent file naming
    def get_oriented_path(group: str, uniprot: str, pdb_id: str, metadata: OrientationMetadata) -> Path:
        """Génère le chemin pour une structure orientée."""
        subdir = out_root_path / f"{group}_{uniprot}"
        subdir.mkdir(parents=True, exist_ok=True)
        suffix = "_native" if metadata.source in ("OPM", "PDBTM") else ""
        return subdir / f"oriented_{group}_{uniprot}_{pdb_id}{suffix}.pdb"

    for group, pdb_ids in groups.items():
        # Get group-specific data
        group_df = df_det_all[df_det_all["pdb_id"].isin(pdb_ids)].copy()
        
        # Find native references
        native_ids = [pid for pid in pdb_ids if source_map.get(pid) in {"OPM", "PDBTM"}]
        preferred_native, native_metadata = _choose_native_reference(group_df, source_map, cfg) if native_ids else (None, OrientationMetadata(source="none"))

        native_done = False

        if preferred_native:
            ref_row = df_ali_all[df_ali_all["pdb_id"] == preferred_native]
            if ref_row.empty:
                logger.warning("Native reference %s missing from Alignment_Summary; skipping group %s", preferred_native, group)
            else:
                det_row_native = df_det_all[df_det_all["pdb_id"] == preferred_native]
                if det_row_native.empty:
                    logger.warning("Native reference %s missing from PDB_Details; skipping group %s", preferred_native, group)
                else:
                    gene_native = det_row_native.iloc[0].get("gene_name", "")
                    uniprot_native = det_row_native.iloc[0].get("uniprot_id", "") or ""
                    if gene_native:
                        ref_path = Path(ref_row.iloc[0]["patched_pdb"])
                        out_dir = out_root_path / f"{gene_native}_{uniprot_native}"
                        out_dir.mkdir(parents=True, exist_ok=True)
                        ref_source = source_map.get(preferred_native, "Native")
                        native_done = True
                        for pid in pdb_ids:
                            row = df_ali_all[df_ali_all["pdb_id"] == pid]
                            if row.empty:
                                continue
                            det_row = df_det_all[df_det_all["pdb_id"] == pid]
                            gene = det_row.iloc[0].get("gene_name", gene_native) if not det_row.empty else gene_native
                            uniprot = det_row.iloc[0].get("uniprot_id", uniprot_native) if not det_row.empty else uniprot_native
                            if not gene:
                                continue
                            target_path = Path(row.iloc[0]["patched_pdb"])
                            out_path = out_dir / f"oriented_{gene}_{uniprot}_{pid}_native.pdb"
                            if pid == preferred_native:
                                shutil.copy(ref_path, out_path)
                                rms_before = rms_after = None
                            else:
                                try:
                                    rms_before, rms_after, z_flipped, transform = _align_to_reference(ref_path, target_path, out_path)
                                except Exception as exc:
                                    logger.error("Native alignment failed for %s: %s", pid, exc)
                                    native_done = False
                                    break
                            method_source = source_map.get(pid, ref_source)
                            method = f"Native_{method_source}" if method_source else "Native"
                            all_aligned_results.append((pid, str(out_path), rms_before, rms_after, method))
                        if native_done:
                            continue


        ref_id = references.get(group)

        if ref_id:
            # Try to find a native reference (OPM/PDBTM) first
            source_map = {}
            for pid in pdb_ids:
                source = df_det_all.loc[df_det_all["pdb_id"] == pid, "pdb_source"].iloc[0] if not df_det_all[df_det_all["pdb_id"] == pid].empty else ""
                if source:
                    source_map[pid] = source
                    
            # Get group-specific data for native reference selection
            group_df = df_det_all[df_det_all["pdb_id"].isin(pdb_ids)].copy()
            native_ref_id, native_metadata = _choose_native_reference(group_df, source_map, cfg)
            
            if native_ref_id:
                # Use native reference first if available
                logger.info(f"Processing group {group} with native reference {native_ref_id} (source: {native_metadata.source})")
                ref_id = native_ref_id
            
            ref_row = df_ali_all[df_ali_all["pdb_id"] == ref_id]
            if ref_row.empty:
                logger.warning(
                    "Reference %s missing from Alignment_Summary; skipping group %s",
                    ref_id,
                    group,
                )
                continue
            ref_path = Path(ref_row.iloc[0]["patched_pdb"])
            if not ref_path.is_file():
                logger.error(
                    f"Reference file missing: {ref_path}; skipping group {group}"
                )
                continue

            target_pdb_ids = [pid for pid in pdb_ids if pid != ref_id]
            if not target_pdb_ids:
                logger.warning(f"No target structures for group {group}")
                # Still copy the reference
            else:
                df_targets = pd.merge(
                    df_ali_all[df_ali_all["pdb_id"].isin(target_pdb_ids)],
                    df_det_all[["pdb_id", "gene_name", "uniprot_id"]],
                    on="pdb_id",
                    how="left",
                )
                uniprot_cols = [c for c in df_targets.columns if c.startswith("uniprot_id")]
                if uniprot_cols:
                    df_targets["uniprot_id"] = (
                        df_targets[uniprot_cols]
                        .bfill(axis=1)
                        .iloc[:, 0]
                    )
                    extra_cols = [c for c in uniprot_cols if c != "uniprot_id"]
                    if extra_cols:
                        df_targets = df_targets.drop(columns=extra_cols)
                df_targets = df_targets.dropna(subset=["patched_pdb"], how="any")
                if not df_targets.empty:
                    results = superimpose_structures(ref_path, df_targets, out_root_path)
                    for res_tuple in results:
                        # Record whether it was aligned to a native reference or DUM reference
                        method = f"Aligned_to_{source_map.get(ref_id)}" if ref_id in source_map else "DUM_Ref"
                        all_aligned_results.append(res_tuple + (method,))

            # Copy reference and add to results
            ref_det_row = df_det_all[df_det_all["pdb_id"] == ref_id]
            if not ref_det_row.empty:
                gene = ref_det_row.iloc[0].get("gene_name", "")
                uniprot = ref_det_row.iloc[0].get("uniprot_id", "") or ""
                if gene:
                    subdir = out_root_path / f"{gene}_{uniprot}"
                    subdir.mkdir(parents=True, exist_ok=True)
                    out_name = f"oriented_{gene}_{uniprot}_{ref_id}.pdb"
                    out_path = subdir / out_name
                    shutil.copy(str(ref_path), str(out_path))
                    logger.info(f"Copied reference {ref_id} -> {out_path}")
                    # Record the orientation method based on source
                    method = f"Native_{source_map.get(ref_id)}" if ref_id in source_map else "DUM_Ref"
                    all_aligned_results.append(
                        (ref_id, str(out_path), None, None, method)
                    )

        elif cfg.use_orientation_fallback and orientation_fallback_fn:
            # --- No DUM Reference: Fallback Orientation ---
            logger.info(f"No DUM reference for group {group}. Attempting to create one with fallback.")

            new_ref_path: Optional[Path] = None
            new_ref_id: Optional[str] = None
            method_used = "Failed"

            # Find the first valid PDB in the group to orient and use as a reference
            for candidate_ref_id in pdb_ids:
                row = df_ali_all[df_ali_all["pdb_id"] == candidate_ref_id]
                if row.empty:
                    continue

                pdb_path = Path(row.iloc[0]["patched_pdb"])
                if not pdb_path.is_file():
                    continue

                det_row = df_det_all[df_det_all["pdb_id"] == candidate_ref_id]
                if det_row.empty:
                    continue

                gene = det_row.iloc[0].get("gene_name", "")
                uniprot = det_row.iloc[0].get("uniprot_id", "") or ""
                if not gene:
                    continue
                
                out_dir = out_root_path / f"{gene}_{uniprot}"
                out_dir.mkdir(parents=True, exist_ok=True)

                # Use the provided fallback function to orient this candidate reference
                oriented_path, method = orientation_fallback_fn(candidate_ref_id, pdb_path, out_dir, cfg)

                if oriented_path:
                    new_ref_path = oriented_path
                    new_ref_id = candidate_ref_id
                    method_used = method
                    logger.info(f"Group {group}: created new reference {new_ref_id} using {method_used}.")
                    all_aligned_results.append(
                        (new_ref_id, str(new_ref_path), None, None, method_used)
                    )
                    break  # We have our reference, stop searching.

            if new_ref_path and new_ref_id:
                # Now superimpose the rest of the group onto this new reference
                target_pdb_ids = [pid for pid in pdb_ids if pid != new_ref_id]
                if target_pdb_ids:
                    df_targets = pd.merge(
                        df_ali_all[df_ali_all["pdb_id"].isin(target_pdb_ids)],
                        df_det_all[["pdb_id", "gene_name", "uniprot_id"]],
                        on="pdb_id",
                        how="left",
                    )
                    uniprot_cols = [c for c in df_targets.columns if c.startswith("uniprot_id")]
                    if uniprot_cols:
                        df_targets["uniprot_id"] = df_targets[uniprot_cols].bfill(axis=1).iloc[:, 0]
                        extra_cols = [c for c in uniprot_cols if c != "uniprot_id"]
                        if extra_cols:
                            df_targets = df_targets.drop(columns=extra_cols)
                    df_targets = df_targets.dropna(subset=["patched_pdb"], how="any")
                    if not df_targets.empty:
                        results = superimpose_structures(new_ref_path, df_targets, out_root_path)
                        for res_tuple in results:
                            all_aligned_results.append(res_tuple + (f"Fallback_Ref_{method_used}",))
            else:
                logger.warning(
                    f"Could not create a fallback reference for group {group}. "
                    "All enabled orientation methods (e.g., PDBTM, MemProtMD, Memembed) failed for all candidate structures in this group. "
                    "Check service configurations in config.py and network connectivity."
                )

    df_res = pd.DataFrame(
        all_aligned_results,
        columns=[
            "pdb_id",
            "aligned_path",
            "rmsd_before",
            "rmsd_after",
            "orientation_method",
        ],
    )

    if "repaired" in df_det_all.columns:
        df_res = pd.merge(
            df_res,
            df_det_all[["pdb_id", "repaired"]],
            on="pdb_id",
            how="left",
        )

    write_back = excel_path is not None and data_frames is None
    if write_back and not df_res.empty:
        with pd.ExcelWriter(
            excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace"
        ) as writer:
            df_res.to_excel(writer, sheet_name="Alignment_RMSD_MultiRef", index=False)
        logger.info(f"Results saved to 'Alignment_RMSD_MultiRef' in {excel_path}")

    return all_results, df_res


# ---------------------------------------------------------------------------
# Fonction principale pour l'utilisation en ligne de commande
# ---------------------------------------------------------------------------


def main():
    """Command-line entry point for multi-reference alignment.

    Raises:
        SystemExit: Propagates parsing errors from ``argparse``."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Aligne et oriente les fichiers PDB par groupes "
            "sur des références contenant DUM"
        )
    )
    parser.add_argument(
        "-e", "--excel", required=True, help="Fichier Excel contenant les données"
    )
    parser.add_argument(
        "-g",
        "--group-by",
        choices=["gene_name", "uniprot_id"],
        default="gene_name",
        help="Attribut pour le regroupement (gene_name ou uniprot_id)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="pdb_structures_oriented",
        help="Répertoire racine pour les structures alignées",
    )

    args = parser.parse_args()

    # Configure le logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Ex??cute l'alignement
    excel_file = Path(args.excel)
    out_root = Path(args.output)

    align_structures_by_groups(excel_file, args.group_by, out_root)


if __name__ == "__main__":
    main()
