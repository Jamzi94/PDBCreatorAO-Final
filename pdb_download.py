import logging
from pathlib import Path

from networking import get_session

log = logging.getLogger(__name__)


def download_pdb_structure(
    pdb_id: str, output_dir: Path | str, filename: str | None = None
) -> tuple[str, str]:
    """Download a PDB file from OPM/PDBTM/RCSB and report its provenance.

    Args:
        pdb_id (str): Four-character PDB identifier to download.
        output_dir (Path | str): Directory where the file will be saved.
        filename (str | None): Optional filename; defaults to `{pdb_id}.pdb`.

    Returns:
        tuple[str, str]: `(absolute_path, source)` where `source` indicates the
            repository that provided the file (`OPM`, `PDBTM`, `RCSB`).
            Returns ("", "") when every source fails.
    """
    pdb_id = pdb_id.lower()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not filename:
        filename = f"{pdb_id}.pdb"
    output_path = output_dir / filename

    # 1) Try OPM mirror (already membrane-oriented)
    opm_url = f"https://biomembhub.org/shared/opm-assets/pdb/{pdb_id}.pdb"
    try:
        resp = get_session().get(opm_url, timeout=10)
        if resp.ok:
            output_path.write_bytes(resp.content)
            log.info(f"Downloaded PDB {pdb_id} from OPM to {output_path}")
            return str(output_path), "OPM"
        log.debug(f"OPM returned status {resp.status_code} for {pdb_id}")
    except Exception as exc:
        log.warning(f"OPM download failed for {pdb_id}: {exc}")

    # 2) Fallback to PDBTM (UniTmp) oriented PDB (.trpdb is PDB text)
    pdbtm_url = f"https://pdbtm.unitmp.org/api/v1/entry/{pdb_id}.trpdb"
    try:
        resp = get_session().get(pdbtm_url, timeout=15)
        if resp.ok:
            head = resp.content[:4096]
            if any(token in head for token in (b"ATOM", b"HETATM", b"HEADER")):
                output_path.write_bytes(resp.content)
                log.info(f"Downloaded oriented PDB {pdb_id} from PDBTM to {output_path}")
                return str(output_path), "PDBTM"
            log.debug(f"PDBTM content for {pdb_id} does not look like PDB text")
        else:
            log.debug(f"PDBTM returned status {resp.status_code} for {pdb_id}")
    except Exception as exc:
        log.info(f"PDBTM retrieval failed for {pdb_id}: {exc}")

    # 3) Fallback to RCSB files server
    rcsb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        resp = get_session().get(rcsb_url, timeout=10)
        resp.raise_for_status()
        output_path.write_bytes(resp.content)
        log.info(f"Downloaded PDB {pdb_id} from RCSB to {output_path}")
        return str(output_path), "RCSB"
    except Exception as exc:
        log.error(f"Failed to download PDB {pdb_id} from RCSB: {exc}")
        return "", ""

