# Protein Structure Analysis Pipeline (PDBCreatorAO)

## English

### Overview
The **PDBCreatorAO** pipeline automates end-to-end protein structure analysis. For a list of UniProt accessions it
retrieves protein metadata, downloads experimental and AlphaFold models, assesses structural quality, aligns and
patches AlphaFold predictions onto experimental structures, orients the resulting models against DUM-containing
references, and produces a curated Excel report plus organized output directories.

### Installation
1. Clone the repository.
2. *(Optional)* Create and activate a virtual environment.
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. *(Optional)* Export `ALPHAFOLD_API_KEY` if you have access to the AlphaFold API.

### Usage
1. Fill `uniprot_template.xlsx` with UniProt IDs in the first column (header row is optional).
2. Run the pipeline from the project root:
   ```bash
   python pipeline.py
   ```
   or supply accessions directly:
   ```bash
   python pipeline.py --ids P12345 Q8ZIN0
   ```
3. Generated artefacts appear under the chosen base directory (default: current working directory).

### CLI Arguments
| Flag | Description |
| ### Repair Toolchain
- Prefer the provided conda/micromamba environment (setup_env.ps1 / environment.yml) to obtain compatible pdbfixer/openmm wheels.
- The pipeline logs a single warning and skips repairs when run on unsupported interpreters (Python ?3.13).

### Orientation Pipeline & Documentation
- Orientation cascade: see ORIENT.md for the fallback decision tree and metadata.
- PDBFixer/OpenMM guidance: see MD_SETUP.md.
- Release notes: see CHANGELOG.md.

### Testing
```bash
pytest
```
--- | ### Repair Toolchain
- Prefer the provided conda/micromamba environment (setup_env.ps1 / environment.yml) to obtain compatible pdbfixer/openmm wheels.
- The pipeline logs a single warning and skips repairs when run on unsupported interpreters (Python ?3.13).

### Orientation Pipeline & Documentation
- Orientation cascade: see ORIENT.md for the fallback decision tree and metadata.
- PDBFixer/OpenMM guidance: see MD_SETUP.md.
- Release notes: see CHANGELOG.md.

### Testing
```bash
pytest
```
--- |
| `--base PATH` | Base working directory for cache, logs, outputs. |
| `--input PATH` | Excel template containing UniProt IDs (default: `uniprot_template.xlsx`). |
| `--output PATH` | Output Excel filename (default: `uniprot_results.xlsx`). |
| `--window INT` | Sliding window length used during patching (default: 5). |
| `--ids ID ...` | Process the listed UniProt IDs instead of reading the Excel template. |
| `-v/--verbose` | Enable verbose logging (DEBUG to console). |

### Outputs
- `pipeline.log`: Complete execution log.
- `cache/uniprot/`: JSON cache of UniProt API responses.
- `pdb_structures_exp/GENE_UNIPROT/`: Experimental PDB downloads.
- `pdb_structures_pred/GENE_UNIPROT/`: AlphaFold predictions and partner models.
- `pdb_structures_af_aligned/GENE/`: Patched AlphaFold?experimental hybrids.
- `pdb_structures_oriented/GENE/`: Structures reoriented against DUM references.
- `uniprot_results.xlsx`: Multi-sheet report containing:
  - **UniProt Results** ? primary metadata with hyperlinks.
  - **PDB_Details** ? experimental structure information (method, resolution, environment, ligands, literature).
  - **Structural_Integrity** ? detected gaps/coils with residue ranges.
  - **Alignment_Summary** ? alignment quality metrics and generated file paths.
  - **Alignment_RMSD_MultiRef** ? RMSD values computed during multi-reference orientation.

### Repair Toolchain
- Prefer the provided conda/micromamba environment (setup_env.ps1 / environment.yml) to obtain compatible pdbfixer/openmm wheels.
- The pipeline logs a single warning and skips repairs when run on unsupported interpreters (Python ?3.13).

### Orientation Pipeline & Documentation
- Orientation cascade: see ORIENT.md for the fallback decision tree and metadata.
- PDBFixer/OpenMM guidance: see MD_SETUP.md.
- Release notes: see CHANGELOG.md.

### Testing
```bash
pytest
```
---

## Fran?ais

### Aper?u
**PDBCreatorAO** automatise l'analyse structurale des prot?ines. ? partir d'identifiants UniProt, il r?cup?re les
m?tadonn?es, t?l?charge les structures exp?rimentales et AlphaFold, ?value la qualit?, aligne et raccommode les
mod?les AlphaFold sur les structures exp?rimentales, oriente les mod?les via des r?f?rences contenant DUM et produit
un rapport Excel accompagn? de dossiers de sortie structur?s.

### Installation
1. Cloner ce d?p?t.
2. *(Optionnel)* Cr?er et activer un environnement virtuel.
3. Installer les d?pendances :
   ```bash
   pip install -r requirements.txt
   ```
4. *(Optionnel)* D?finir `ALPHAFOLD_API_KEY` si vous disposez d'une cl? AlphaFold.

### Utilisation
1. Compl?ter `uniprot_template.xlsx` avec les identifiants UniProt dans la premi?re colonne (ligne d'ent?te optionnelle).
2. Ex?cuter le pipeline depuis la racine du projet :
   ```bash
   python pipeline.py
   ```
   ou fournir les identifiants directement :
   ```bash
   python pipeline.py --ids P12345 Q8ZIN0
   ```
3. Les artefacts sont g?n?r?s dans le r?pertoire de base choisi (par d?faut : dossier courant).

### Arguments CLI
| Option | Description |
| ### Repair Toolchain
- Prefer the provided conda/micromamba environment (setup_env.ps1 / environment.yml) to obtain compatible pdbfixer/openmm wheels.
- The pipeline logs a single warning and skips repairs when run on unsupported interpreters (Python ?3.13).

### Orientation Pipeline & Documentation
- Orientation cascade: see ORIENT.md for the fallback decision tree and metadata.
- PDBFixer/OpenMM guidance: see MD_SETUP.md.
- Release notes: see CHANGELOG.md.

### Testing
```bash
pytest
```
--- | ### Repair Toolchain
- Prefer the provided conda/micromamba environment (setup_env.ps1 / environment.yml) to obtain compatible pdbfixer/openmm wheels.
- The pipeline logs a single warning and skips repairs when run on unsupported interpreters (Python ?3.13).

### Orientation Pipeline & Documentation
- Orientation cascade: see ORIENT.md for the fallback decision tree and metadata.
- PDBFixer/OpenMM guidance: see MD_SETUP.md.
- Release notes: see CHANGELOG.md.

### Testing
```bash
pytest
```
--- |
| `--base PATH` | R?pertoire de travail pour le cache, les journaux et les sorties. |
| `--input PATH` | Mod?le Excel contenant les identifiants UniProt (d?faut : `uniprot_template.xlsx`). |
| `--output PATH` | Fichier Excel de sortie (d?faut : `uniprot_results.xlsx`). |
| `--window INT` | Taille de fen?tre utilis?e pendant le raccommodage (d?faut : 5). |
| `--ids ID ...` | Traiter les identifiants UniProt fournis au lieu de lire le mod?le Excel. |
| `-v/--verbose` | Active un journal d?taill? (niveau DEBUG sur la console). |

### Sorties
- `pipeline.log` : journal complet de l'ex?cution.
- `cache/uniprot/` : cache JSON des r?ponses UniProt.
- `pdb_structures_exp/GENE_UNIPROT/` : structures PDB exp?rimentales.
- `pdb_structures_pred/GENE_UNIPROT/` : pr?dictions AlphaFold et mod?les partenaires.
- `pdb_structures_af_aligned/GENE/` : structures exp?rimentales raccommod?es avec AlphaFold.
- `pdb_structures_oriented/GENE/` : structures r?orient?es via une r?f?rence contenant DUM.
- `uniprot_results.xlsx` : rapport multi-feuilles comprenant :
  - **UniProt Results** ? m?tadonn?es principales avec hyperliens.
  - **PDB_Details** ? informations exp?rimentales (m?thode, r?solution, environnement, ligands, r?f?rences).
  - **Structural_Integrity** ? gaps/coils d?tect?s avec plages de r?sidus.
  - **Alignment_Summary** ? m?triques d'alignement et chemins des fichiers g?n?r?s.
  - **Alignment_RMSD_MultiRef** ? valeurs de RMSD obtenues lors de l'orientation multi-r?f?rence.


### Cha?ne de r?paration
- Utilisez de pr?f?rence l'environnement conda/micromamba fourni (`setup_env.ps1` / `environment.yml`) afin d'obtenir des versions compatibles de `pdbfixer`/`openmm`.
- Le pipeline journalise un avertissement unique et ignore la r?paration si l'interpr?teur n'est pas support? (Python ?3.13).

### Orientation & documentation
- Cascade d'orientation : consultez `ORIENT.md` pour le d?tail des ?tapes et des m?tadonn?es.
- Guide PDBFixer/OpenMM : voir `MD_SETUP.md`.
- Notes de version : voir `CHANGELOG.md`.

### Tests
```bash
pytest
```
