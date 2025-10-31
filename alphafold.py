import logging
import os

from pathlib import Path
from networking import get_session

log = logging.getLogger(__name__)

# If you have an API key, set it as environment variable ALPHAFOLD_API_KEY;
# otherwise you can omit the key parameter.
API_KEY = os.getenv("ALPHAFOLD_API_KEY", "")


def _prediction_url(uniprot_id: str) -> str:
    base_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    return f"{base_url}?key={API_KEY}" if API_KEY else base_url


def _fetch_predictions(uniprot_id: str) -> list:
    """Fetch the AlphaFold prediction list for a UniProt accession.

    Args:
        uniprot_id (str): Target UniProt accession.

    Returns:
        list: API response containing zero or more prediction dictionaries.
    """
    url = _prediction_url(uniprot_id)
    try:
        resp = get_session().get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        if not isinstance(data, list):
            log.error(f"Unexpected AlphaFold response format for {uniprot_id}")
            return []
        return data
    except Exception as e:
        log.error(f"Error fetching AlphaFold prediction for {uniprot_id}: {e}")
        return []


def get_alphafold_pdb_info(uniprot_id: str) -> str:
    """Return the AlphaFold PDB identifier for a UniProt accession.

    Args:
        uniprot_id (str): UniProt accession identifier.

    Returns:
        str: AlphaFold entry ID such as ``AF-P14859-F1``.
            Returns an empty string when no prediction is available.
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds:
        log.info(f"No AlphaFold prediction available for {uniprot_id}")
        return ""
    entry_id = preds[0].get("entryId", "")
    return entry_id


def download_alphafold_pdb(uniprot_id: str, output_path: str) -> str:
    """Download the AlphaFold PDB file to ``output_path``.

    Args:
        uniprot_id (str): UniProt accession to download.
        output_path (str): Destination path on disk.

    Returns:
        str: Absolute path to the downloaded file, or an empty string on failure.
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds or not preds[0].get("pdbUrl"):
        log.warning(f"No PDB URL for AlphaFold model of {uniprot_id}")
        return ""
    pdb_url = preds[0]["pdbUrl"]
    destination = Path(output_path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        resp = get_session().get(pdb_url, timeout=10)
        resp.raise_for_status()
        destination.write_bytes(resp.content)
        log.info(f"Downloaded AlphaFold PDB for {uniprot_id} to {destination}")
        return str(destination)
    except Exception as e:
        log.error(f"Error downloading AlphaFold PDB for {uniprot_id}: {e}")
        return ""


def get_alphafold_details(uniprot_id: str, prot_name: str) -> dict:
    """Summarise AlphaFold metadata for the provided accession.

    Args:
        uniprot_id (str): UniProt accession identifier.
        prot_name (str): Human-readable protein name used when descriptions are absent.

    Returns:
        dict: Details describing the model methodology, counts, chains and organisms.
    """
    preds = _fetch_predictions(uniprot_id)
    if not preds:
        return {
            "method": "Prediction",
            "resolution": "N/A",
            "methodology": "AlphaFold",
            "conformation": "Not specified",
            "environment": "NA",
            "count": 0,
            "descs": prot_name,
            "chains": "",
            "orgs": "",
        }

    p = preds[0]
    method = "Prediction"
    resolution = "N/A"
    methodology = "AlphaFold"
    conformation = "Not specified"
    environment = "NA"

    # You could decide to infer conformation from templateModelQuality or other fields
    # but by default AlphaFold predictions have no explicit 'conformation'.
    entry_id = p.get("entryId", "")
    count = 1 if entry_id else 0

    # Description: use the protein name passed in
    descs = prot_name

    # Chains: AlphaFold models are single-chain, chain ID 'A'
    chains = "A"

    # Organism name from the returned JSON
    orgs = p.get("organismScientificName", "")

    return {
        "method": method,
        "resolution": resolution,
        "methodology": methodology,
        "conformation": conformation,
        "environment": environment,
        "count": count,
        "descs": descs,
        "chains": chains,
        "orgs": orgs,
    }
