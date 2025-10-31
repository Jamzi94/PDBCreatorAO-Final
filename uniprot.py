import json
import logging
from pathlib import Path
from typing import Dict, List

from networking import get_session

BASE_URL = "https://rest.uniprot.org/uniprotkb"
FIELDS = ",".join(
    [
        "gene_primary",
        "protein_name",
        "xref_pdb",
        "ft_variant",
        "cc_pharmaceutical",
        "xref_tcdb",
        "xref_drugcentral",
        "cc_tissue_specificity",
        "annotation_score",
        "cc_cofactor",
        "cc_polymorphism",
        "sequence",
        "cc_interaction",
    ]
)

log = logging.getLogger(__name__)

# Directory for caching UniProt JSON results
CACHE_DIR = Path.cwd() / "cache" / "uniprot"
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def fetch_uniprot_json(uniprot_id: str) -> dict:
    """Retrieve the UniProt JSON payload for an accession, using on-disk caching.

    Args:
        uniprot_id (str): UniProt accession to download.

    Returns:
        dict: Parsed JSON document describing the UniProt entry.

    Raises:
        requests.HTTPError: If the remote API responds with an error status.
    """
    uniprot_id = uniprot_id.strip()
    cache_file = CACHE_DIR / f"{uniprot_id}.json"
    if cache_file.is_file():
        try:
            with cache_file.open("r", encoding="utf-8") as handle:
                data = json.load(handle)
                return data
        except Exception as exc:
            log.warning(
                f"Failed to load cache for {uniprot_id}: {exc}. Fetching fresh."
            )
    url = f"{BASE_URL}/{uniprot_id}"
    resp = get_session().get(
        url,
        params={"fields": FIELDS},
        headers={"Accept": "application/json"},
        timeout=10,
    )
    resp.raise_for_status()
    data = resp.json()
    try:
        with cache_file.open("w", encoding="utf-8") as handle:
            json.dump(data, handle)
    except Exception as exc:
        log.warning(f"Failed to write cache for {uniprot_id}: {exc}")
    return data


def get_first_comment(comments: List[dict], comment_type: str) -> str:
    """Return the first comment matching ``comment_type`` from a UniProt payload.

    Args:
        comments (List[dict]): Comment entries from the UniProt JSON payload.
        comment_type (str): Targeted UniProt comment type (e.g. ``PHARMACEUTICAL``).

    Returns:
        str: Extracted comment text, or an empty string when not present.
    """
    for c in comments:
        if c.get("commentType") == comment_type:
            texts = c.get("texts", [])
            if texts:
                val = texts[0]
                if isinstance(val, dict) and "value" in val:
                    return val["value"]
                elif isinstance(val, str):
                    return val
    return ""


def get_interactions(uniprot_id: str) -> dict:
    """
    Récupère uniquement le bloc cc_interaction.
    """
    url = f"{BASE_URL}/{uniprot_id}"
    resp = get_session().get(
        url,
        params={"fields": "cc_interaction"},
        headers={"accept": "application/json"},
        timeout=10,
    )
    resp.raise_for_status()
    return resp.json()


def find_related_uniprots(entry_uniprot: str, data: dict) -> List[Dict[str, str]]:
    """Discover interaction partners that share the accession prefix.

    Args:
        entry_uniprot (str): Reference UniProt accession.
        data (dict): Parsed UniProt JSON payload.

    Returns:
        List[Dict[str, str]]: Partner records containing UniProt and gene identifiers.
    """
    related = []
    prefix = entry_uniprot[:5]
    for c in data.get("comments", []):
        if c.get("commentType") != "INTERACTION":
            continue
        for inter in c.get("interactions", []):
            p = inter.get("interactantTwo", {})
            pid = p.get("uniProtKBAccession", "")
            pg = p.get("geneName", "")
            # Only add entries with both UniProt ID and gene name
            if not pid or "-" in pid or not pg:
                continue
            if pid.startswith(prefix):
                related.append({"uniprot_id": pid, "gene_name": pg})
    return related


def parse_uniprot_data(data: dict) -> dict:
    """Flatten UniProt JSON into the structure consumed by the pipeline.

    Args:
        data (dict): JSON document returned by the UniProt REST API.

    Returns:
        dict: Key fields for downstream processing.
            Includes sequence and interaction information.
    """
    # Gene name
    gene = ""
    try:
        gene = data.get("genes", [{}])[0].get("geneName", {}).get("value", "")
    except Exception:
        log.debug("No gene name found in UniProt data")
    # Protein name
    prot = ""
    try:
        prot = (
            data.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        )
    except Exception:
        log.debug("No recommended protein name found")
    # PDB cross-references
    pdb_ids = [
        x.get("id")
        for x in data.get("uniProtKBCrossReferences", [])
        if x.get("database") == "PDB" and x.get("id")
    ]
    # Natural variant (SNP) features
    variants_list = []
    for feat in data.get("features", []):
        if feat.get("type") == "Natural variant":
            vid = feat.get("featureId", "")
            loc = feat.get("location", {}).get("start", {}).get("value", "")
            alt = feat.get("alternativeSequence", {})
            orig = alt.get("originalSequence", "")
            new = ""
            if alt.get("alternativeSequences"):
                new = alt["alternativeSequences"][0]
            if vid and loc:
                variants_list.append(f"{vid}({loc}{orig}>{new})")
    variants = "; ".join(variants_list) if variants_list else ""
    # Comments
    pharm = get_first_comment(data.get("comments", []), "PHARMACEUTICAL")
    tissue = get_first_comment(data.get("comments", []), "TISSUE SPECIFICITY")
    cofactor = get_first_comment(data.get("comments", []), "COFACTOR")
    polymorphism = get_first_comment(data.get("comments", []), "POLYMORPHISM")
    # Cross-refs TCDB / DrugCentral
    tcdb_id, drugcentral_id = "", ""
    for x in data.get("uniProtKBCrossReferences", []):
        if x.get("database") == "TCDB":
            tcdb_id = x.get("id", "")
        if x.get("database") == "DrugCentral":
            drugcentral_id = x.get("id", "")
    # Annotation score & sequence
    annotation_score = data.get("annotationScore", "")
    sequence = data.get("sequence", {}).get("value", "")
    # Add interactions
    related = find_related_uniprots(data.get("primaryAccession", ""), data)
    related_ids = ""
    related_genes = ""
    if related:
        # If partners found, concatenate gene names
        related_genes = "-".join([gene] + [r["gene_name"] for r in related])
        related_ids = ";".join([r["uniprot_id"] for r in related])
    else:
        related_genes = gene
    return {
        "gene": related_genes,
        "protein_name": prot,
        "pdb_ids": pdb_ids,
        "variants": variants,
        "pharmaceutical": pharm,
        "tcdb_id": tcdb_id,
        "drugcentral_id": drugcentral_id,
        "tissue": tissue,
        "annotation_score": annotation_score,
        "cofactor": cofactor,
        "polymorphism": polymorphism,
        "sequence": sequence,
        "interaction_partners": related_ids,
    }
