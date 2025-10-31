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
|------|-------------|
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
- `pdb_structures_af_aligned/GENE/`: Patched AlphaFold–experimental hybrids.
- `pdb_structures_oriented/GENE/`: Structures reoriented against DUM references.
- `uniprot_results.xlsx`: Multi-sheet report containing:
  - **UniProt Results** — primary metadata with hyperlinks.
  - **PDB_Details** — experimental structure information (method, resolution, environment, ligands, literature).
  - **Structural_Integrity** — detected gaps/coils with residue ranges.
  - **Alignment_Summary** — alignment quality metrics and generated file paths.
  - **Alignment_RMSD_MultiRef** — RMSD values computed during multi-reference orientation.

### Repair Toolchain
- Prefer the provided conda/micromamba environment (setup_env.ps1 / environment.yml) to obtain compatible pdbfixer/openmm wheels.
- The pipeline logs a single warning and skips repairs when run on unsupported interpreters (Python ≥3.13).

### Orientation Pipeline & Documentation
- Orientation cascade: see ORIENT.md for the fallback decision tree and metadata.
- PDBFixer/OpenMM guidance: see MD_SETUP.md.
- Release notes: see CHANGELOG.md.

### Testing
```bash
pytest
```
---

## Français

### Aperçu
**PDBCreatorAO** automatise l'analyse structurale des protéines. À partir d'identifiants UniProt, il récupère les
métadonnées, télécharge les structures expérimentales et AlphaFold, évalue la qualité, aligne et raccommode les
modèles AlphaFold sur les structures expérimentales, oriente les modèles via des références contenant DUM et produit
un rapport Excel accompagné de dossiers de sortie structurés.

### Installation
1. Cloner ce dépôt.
2. *(Optionnel)* Créer et activer un environnement virtuel.
3. Installer les dépendances :
   ```bash
   pip install -r requirements.txt
   ```
4. *(Optionnel)* Définir `ALPHAFOLD_API_KEY` si vous disposez d'une clé AlphaFold.

### Utilisation
1. Compléter `uniprot_template.xlsx` avec les identifiants UniProt dans la première colonne (ligne d'entête optionnelle).
2. Exécuter le pipeline depuis la racine du projet :
   ```bash
   python pipeline.py
   ```
   ou fournir les identifiants directement :
   ```bash
   python pipeline.py --ids P12345 Q8ZIN0
   ```
3. Les artefacts sont générés dans le répertoire de base choisi (par défaut : dossier courant).

### Arguments CLI
| Option | Description |
|--------|-------------|
| `--base PATH` | Répertoire de travail pour le cache, les journaux et les sorties. |
| `--input PATH` | Modèle Excel contenant les identifiants UniProt (défaut : `uniprot_template.xlsx`). |
| `--output PATH` | Fichier Excel de sortie (défaut : `uniprot_results.xlsx`). |
| `--window INT` | Taille de fenêtre utilisée pendant le raccommodage (défaut : 5). |
| `--ids ID ...` | Traiter les identifiants UniProt fournis au lieu de lire le modèle Excel. |
| `-v/--verbose` | Active un journal détaillé (niveau DEBUG sur la console). |

### Sorties
- `pipeline.log` : journal complet de l'exécution.
- `cache/uniprot/` : cache JSON des réponses UniProt.
- `pdb_structures_exp/GENE_UNIPROT/` : structures PDB expérimentales.
- `pdb_structures_pred/GENE_UNIPROT/` : prédictions AlphaFold et modèles partenaires.
- `pdb_structures_af_aligned/GENE/` : structures expérimentales raccommodées avec AlphaFold.
- `pdb_structures_oriented/GENE/` : structures réorientées via une référence contenant DUM.
- `uniprot_results.xlsx` : rapport multi-feuilles comprenant :
  - **UniProt Results** — métadonnées principales avec hyperliens.
  - **PDB_Details** — informations expérimentales (méthode, résolution, environnement, ligands, références).
  - **Structural_Integrity** — gaps/coils détectés avec plages de résidus.
  - **Alignment_Summary** — métriques d'alignement et chemins des fichiers générés.
  - **Alignment_RMSD_MultiRef** — valeurs de RMSD obtenues lors de l'orientation multi-référence.


### Chaîne de réparation
- Utilisez de préférence l'environnement conda/micromamba fourni (`setup_env.ps1` / `environment.yml`) afin d'obtenir des versions compatibles de `pdbfixer`/`openmm`.
- Le pipeline journalise un avertissement unique et ignore la réparation si l'interpréteur n'est pas supporté (Python ≥3.13).

### Orientation & documentation
- Cascade d'orientation : consultez `ORIENT.md` pour le détail des étapes et des métadonnées.
- Guide PDBFixer/OpenMM : voir `MD_SETUP.md`.
- Notes de version : voir `CHANGELOG.md`.

### Tests
```bash
pytest
```
