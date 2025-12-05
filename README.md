# uniprot-blast-viewer

A small toolkit + GUI for:
- submitting BLAST queries to EBI,
- downloading BLAST JSON,
- interactively exploring HSPs with a PyQt GUI,
- querying UniProt, ChEBI, etc.

## Installation (from GitHub)

```bash
pip install git+https://github.com/niladriranjandas/uniprot_blast_viewer.git

# for linux

python3 -m venv /tmp/testenv
source /tmp/testenv/bin/activate
pip install dist/uniprot_blast_viewer-0.1.0-py3-none-any.whl
uniprot-blast-gui
