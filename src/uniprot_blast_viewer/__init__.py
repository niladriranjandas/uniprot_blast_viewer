# __init__.py
# Initializes the package and exposes key functions and globals.

from .uniprot_querry import getFasta_from_accession, doBlast_from_accession, doBlast_from_sequence, see_Blast_status, see_Blast_result, see_Blast_result_json, get_isoforms, get_isoforms_list, align_all_accessions, see_alignJob_status, see_alignJob_result, get_aligned_results_accessions, get_PTMs_from_accession, get_PTMs_from_accession_list, import_BindingInfo_uniProt, get_mass_chEBI
from .general_utils import aa_list, aa_check, mass_shift, mw, enumerate_terminal_ptm_masses,fasta_to_sequence, calculate_peptide_masses, get_url, get_url_blast, list_to_dataframe, sequence_to_fasta, tsv_to_file, tsv_to_list, df_to_list, df_to_list, parse_chebi_json
from .view_blast_jason_web import run_blast_hsp_viewer
from .globals import WEBSITE_API, PROTEINS_API, EBI_API

__version__ = "1.1.0"
__author__ = "Your Name"

# Optional: convenience dictionary to show all exports

__version__ = "0.1.0"
__author__ = "Your Name"

__all__ = [
    "getFasta_from_accession",
    "doBlast_from_accession",
    "get_url",
    "get_url_blast",
    "run_blast_hsp_viewer",
    "WEBSITE_API",
    "PROTEINS_API",
    "EBI_API",
]
