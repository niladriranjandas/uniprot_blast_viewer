import sys
import json
import requests

import pandas as pd
import re

import json

from pyteomics import mass
from textwrap import wrap

# Helper function to download data
#def get_url(url, **kwargs):
#  response = requests.get(url, **kwargs);

#  if not response.ok:
#    print(response.text)
#    response.raise_for_status()
#    #sys.exit()

#  return response


def aa_list():
    """
    Returns the list of allowed amino acid one-letter codes
    as used by the Peptides package.
    """
    return list("ACDEFGHIKLMNPQRSTVWYUO")

def aa_check(seq):
    """
    Accept either a string sequence or a list of amino acids.
    Converts the sequence into a validated list of uppercase AA letters.
    """
    valid = set(aa_list())

    # If already a list (e.g. ['A','C','D']) → just validate
    if isinstance(seq, list):
        cleaned = []
        for aa in seq:
            aa = aa.upper()
            if aa not in valid:
                raise ValueError(f"Invalid amino acid: {aa}")
            cleaned.append(aa)
        return cleaned

    # Otherwise treat as a string
    if isinstance(seq, str):
        seq = seq.strip().upper()
        for aa in seq:
            if aa not in valid:
                raise ValueError(f"Invalid amino acid: {aa}")
        return list(seq)

    # Wrong datatype
    raise TypeError(f"seq must be a string or list, not {type(seq)}")



def mass_shift(seq, label="none", aaShift=None, monoisotopic=True):
    """
    Python equivalent of the R function massShift() from the Peptides package.

    Parameters
    ----------
    seq : str
        Amino-acid sequence.
    label : str
        One of: "none", "silac_13c", "silac_13c15n", "15n".
        Overrides aaShift if not "none".
    aaShift : dict or None
        Custom amino-acid mass shifts. Keys must be amino-acid letters or
        "Nterm" / "Cterm". Values are mass shifts in Dalton.
    monoisotopic : bool
        Whether to use monoisotopic adjustments (TRUE/FALSE).
    """

    # --- Normalize inputs ---
    label = label.lower()

    valid_labels = {"none", "silac_13c", "silac_13c15n", "15n"}
    if label not in valid_labels:
        raise ValueError("Unknown label. Use one of: 'none', 'silac_13c', 'silac_13c15n', '15n'")

    # Check aaShift correctness
    allowed = set(aa_list()) | {"Nterm", "Cterm"}

    if aaShift is not None:
        if not all(isinstance(k, str) for k in aaShift.keys()):
            raise ValueError("aaShift must have string names (amino-acid codes).")

        if not all(k in allowed for k in aaShift.keys()):
            good = ", ".join(sorted(allowed))
            raise ValueError(f"Unknown names in aaShift. Allowed: {good}")

    # R logic: predefined isotope labels override aaShift
    if label == "silac_13c":
        # R code: 6.020129 - 0.064229 * !monoisotopic
        corr = -0.064229 * (not monoisotopic)
        aaShift = {"K": 6.020129 + corr, "R": 6.020129 + corr}

    elif label == "silac_13c15n":
        # R code: K = 8.014199 - 0.071499 * !mono
        #         R = 10.008269 - 0.078669 * !mono
        aaShift = {
            "K": 8.014199 - 0.071499 * (not monoisotopic),
            "R": 10.008269 - 0.078669 * (not monoisotopic),
        }

    elif label == "15n":
        # Exact vector from R code
        base = {
            "A": 1, "R": 4, "N": 2, "D": 1, "C": 1, "E": 1, "Q": 2, "G": 1,
            "H": 3, "I": 1, "L": 1, "K": 2, "M": 1, "F": 1, "P": 1, "S": 1,
            "T": 1, "W": 2, "Y": 1, "V": 1
        }

        # 0.997035 * base - 0.003635 * !monoisotopic
        shift_factor = 0.997035
        offset = -0.003635 * (not monoisotopic)

        aaShift = {aa: base[aa] * shift_factor + offset for aa in base}

    # If label == "none" and no aaShift given: aaShift = empty dict
    if aaShift is None:
        aaShift = {}

    # Ensure Nterm and Cterm keys exist (as in R)
    if "Nterm" not in aaShift:
        aaShift["Nterm"] = 0.0
    if "Cterm" not in aaShift:
        aaShift["Cterm"] = 0.0

    # --- Validate and split sequence ---
    seq_list = aa_check(seq)

    # --- Compute mass shift exactly like R ---
    total_shift = 0.0
    for aa in seq_list:
        aa_mass = aaShift.get(aa, 0.0)
        total_shift += aa_mass

    # Add Nterm and Cterm (R does this inside each lapply)
    total_shift += aaShift.get("Nterm", 0.0)
    total_shift += aaShift.get("Cterm", 0.0)

    return total_shift

def mw(seq, monoisotopic=False, avgScale="expasy",
       label="none", aaShift=None):
    """
    Python equivalent of the R function mw().

    Parameters
    ----------
    seq : str
        Amino acid sequence (1-letter codes).
    monoisotopic : bool
        Use monoisotopic masses if True.
    avgScale : str
        "expasy" (default) or "mascot" if monoisotopic=False.
    label : str
        Predefined isotope label: "none", "silac_13c", 
        "silac_13c15n", "15n".
    aaShift : dict
        Custom mass shifts (e.g., {"C": 57.02}).
    """

    seq = aa_check(seq)

    # --- Define amino acid weight scales ---
    mono = {
        "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
        "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
        "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
        "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
        "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
        "U": 150.95363, "O": 237.14772, "H2O": 18.01056
    }

    avg_expasy = {
        "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
        "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
        "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
        "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
        "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
        "U": 150.0388, "O": 237.3018, "H2O": 18.01524
    }

    avg_mascot = {
        "A": 71.0779, "R": 156.1857, "N": 114.1026, "D": 115.0874,
        "C": 103.1429, "E": 129.1140, "Q": 128.1292, "G": 57.0513,
        "H": 137.1393, "I": 113.1576, "L": 113.1576, "K": 128.1723,
        "M": 131.1961, "F": 147.1739, "P": 97.1152, "S": 87.0773,
        "T": 101.1039, "W": 186.2099, "Y": 163.1733, "V": 99.1311,
        "U": 150.0379, "O": 237.2982, "H2O": 18.01528
    }

    # Choose the mass scale
    if monoisotopic:
        weight = mono
    else:
        if avgScale == "expasy":
            weight = avg_expasy
        elif avgScale == "mascot":
            weight = avg_mascot
        else:
            raise ValueError("avgScale must be 'expasy' or 'mascot'.")

    # --- Compute molecular weight (sum of AA masses + H2O) ---
    mass = sum(weight[aa] for aa in seq) + weight["H2O"]

    # --- Apply isotope mass shift ---
    mass += mass_shift(seq, label=label, aaShift=aaShift,
                       monoisotopic=monoisotopic)

    return mass


# -----------------------
# PTM enumeration and mass calculation
# -----------------------

# Define a small library of N- and C-terminal PTMs.
# Mass shifts are given as monoisotopic values; average values can differ slightly.
# For more exhaustive and precise PTMs refer to Unimod (www.unimod.org).

NTERM_PTMS = {
    "none": 0.0,
    # Example: Acetyl (Unimod:1) ~ +42.010565 Da
    "Acetyl": 42.010565,
    # Example: Pyro-glutamate from N-terminal Q (~ -17.026549, but context-specific)
    "PyroGlu_from_Q": -17.026549,
}

CTERM_PTMS = {
    "none": 0.0,
    # Example: Amide (amidated C-terminus), ~ -0.984016 Da vs free acid
    "Amidated": -0.984016,
    # Example: Methyl ester (approx +14.01565), simplified
    "MethylEster": 14.01565,
}


def enumerate_terminal_ptm_masses(seq):
    """
    For a given peptide, enumerate all combinations of N-term and C-term PTMs
    and compute monoisotopic and average masses.

    Returns
    -------
    list[dict]
        Each dict has:
          - 'nterm'
          - 'cterm'
          - 'mono_mass'
          - 'avg_mass_expasy'
          - 'avg_mass_mascot'
    """
    results = []

    for n_name, n_shift in NTERM_PTMS.items():
        for c_name, c_shift in CTERM_PTMS.items():
            # Construct aaShift with terminal modifications
            aa_shift = {
                "Nterm": n_shift,
                "Cterm": c_shift,
            }

            mono = mw(seq, monoisotopic=True, label="none", aaShift=aa_shift)
            avg_expasy = mw(seq, monoisotopic=False, avgScale="expasy",
                            label="none", aaShift=aa_shift)
            avg_mascot = mw(seq, monoisotopic=False, avgScale="mascot",
                            label="none", aaShift=aa_shift)

            results.append({
                "nterm": n_name,
                "cterm": c_name,
                "mono_mass": mono,
                "avg_mass_expasy": avg_expasy,
                "avg_mass_mascot": avg_mascot,
            })

    return results


def fasta_to_sequence(fasta_text, clean_non_aa=True):
    """
    Convert a FASTA-formatted sequence into a single-line amino acid string.

    Parameters
    ----------
    fasta_text : str
        A string containing FASTA data (can be multi-line, may include header).
    clean_non_aa : bool (default=True)
        If True, remove any characters not valid one-letter AA codes.

    Returns
    -------
    str
        The amino acid sequence in a single line, uppercase.
    """
    # Split into lines and remove header lines starting with '>'
    lines = fasta_text.strip().splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]

    # Join into a single string
    sequence = "".join(seq_lines).upper()

    # Optionally clean invalid characters
    if clean_non_aa:
        valid = set("ACDEFGHIKLMNPQRSTVWYUO")
        sequence = "".join([aa for aa in sequence if aa in valid])

    return sequence


def calculate_peptide_masses(sequence, modifications=None):
    """
    Calculates average mass, monoisotopic mass, and termini masses.

    Args:
        sequence (str): The peptide sequence (e.g., "PEPTIDE").
        modifications (dict): A dictionary mapping amino acids to their
                              modification mass (e.g., {'C': 57.021}).

    Returns:
        dict: A dictionary containing the calculated masses.
    """
    # Calculate average and monoisotopic masses
    mono_mass = mass.calculate_mass(sequence, average_mass=False)
    avg_mass = mass.calculate_mass(sequence, average_mass=True)

    # Calculate termini masses
    # N-terminal mass (monoisotopic)
    n_term_mass = mass.calculate_mass(sequence[1:], average_mass=False)

    # C-terminal mass (monoisotopic)
    c_term_mass = mass.calculate_mass(sequence[:-1], average_mass=False)

    # You can also add modifications using the 'modification' argument
    if modifications:
        mono_mass_mod = mass.calculate_mass(sequence, average_mass=False, modification=modifications)
        avg_mass_mod = mass.calculate_mass(sequence, average_mass=True, modification=modifications)
        
        n_term_mass_mod = mass.calculate_mass(sequence[1:], average_mass=False, modification=modifications)
        c_term_mass_mod = mass.calculate_mass(sequence[:-1], average_mass=False, modification=modifications)
    else:
        mono_mass_mod = 0
        avg_mass_mod = 0
        
        n_term_mass_mod = 0
        c_term_mass_mod = 0


    return {
        "monoisotopic_mass": mono_mass,
        "average_mass": avg_mass,
        "n_terminal_mass": n_term_mass,
        "c_terminal_mass": c_term_mass,
        "monoisotopic_mass_modified": mono_mass_mod,
        "average_mass_modified": avg_mass_mod,
        "n_terminal_mass_modified": n_term_mass_mod,
        "c_terminal_mass_modified": c_term_mass_mod,
    }

def get_url(url, **kwargs):
    """
    Safe version of get_url():
    - Never quits the program
    - Catches HTTP errors (404, 500, etc.)
    - Catches timeout, connection errors, SSL errors
    - Returns None on failure
    """
    try:
        response = requests.get(url, **kwargs)

        # Explicitly raise HTTP errors (raises for 4xx, 5xx)
        response.raise_for_status()

        return response

    except requests.exceptions.HTTPError as e:
        print(f"[HTTP ERROR] {e}")
        print(f"URL: {url}")
        return None

    except requests.exceptions.RequestException as e:
        # Catches: Timeout, ConnectionError, SSLError, too many redirects, etc.
        print(f"[REQUEST ERROR] Could not retrieve URL: {url}")
        print(f"Details: {e}")
        return None

def get_url_blast(url, **kwargs):
  params = {"format": "tsv"}
  response = requests.get(url, **kwargs)

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response

def list_to_dataframe(data: list, headers: list) -> pd.DataFrame:
    """
    Convert a list of lists into a pandas DataFrame with user-specified headers.

    Parameters
    ----------
    data : list
        A list of lists, where each inner list is a row of data.
    headers : list
        A list of column names (strings). Must match the number of columns in data.

    Returns
    -------
    df : pandas.DataFrame
        A DataFrame with the provided headers.
    """
    # Validate that header length matches number of columns
    if len(headers) != len(data[0]):
        raise ValueError(
            f"Number of headers ({len(headers)}) does not match number of columns ({len(data[0])})"
        )

    # Create DataFrame
    df = pd.DataFrame(data, columns=headers)
    return df


def df_to_list(df: pd.DataFrame):
   """
   """
   return df.values.tolist()

# Example: Convert a raw amino acid sequence string into FASTA format
# and save it as a .fasta file for use in BLAST (e.g., UniProt BLAST).

def sequence_to_fasta(seq, identifier="query_protein", line_length=60, filename=''):
    """
    Convert a raw amino acid sequence string to FASTA format and write to a file.

    Parameters
    ----------
    seq : str
        The amino acid sequence (string of residues, e.g., 'MTEYKLVVVG...')
    identifier : str
        A header identifier for the sequence (appears after '>' in the FASTA file)
    line_length : int
        Number of characters per line in the FASTA sequence block (default=60, standard practice)
    filename : str
        The name of the output FASTA file
    
    Returns
    -------
    str
        The FASTA-formatted string (also written to file)
    """
    
    # Clean the sequence (remove whitespace, make uppercase)
    seq = "".join(seq.split()).upper()
    
    # Wrap the sequence into lines of fixed length (default 60 chars)
    wrapped_seq = "\n".join(wrap(seq, line_length))
    
    # Build the FASTA string
    fasta_str = f">{identifier}\n{wrapped_seq}\n"
    
    # Write to file
    if filename:
        with open(filename, "w") as f:
             f.write(fasta_str)

    return fasta_str


def save_tsv_to_csv(tsv_text, output_csv):
    """
    Convert TSV text into a CSV file.
    """
    lines = tsv_text.strip().split("\n")
    reader = csv.reader(lines, delimiter="\t")

    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(reader)

    print(f"[INFO] Results saved to {output_csv}")

# reading given tsv file
def tsv_to_file(tsv_text, filename='tsv_to_csv.csv'):
   """ 
   
   """
  
   with open(filename, 'w') as csv_file:
      for line in tsv_text:
      
         # Replace every tab with comma
         fileContent = re.sub("\t", ",", line)
      
         # Writing into csv file
         csv_file.write(fileContent)

   # output
   print("Successfully made csv file")

# reading given tsv file
def tsv_to_list(tsv_text: str, delimiter=None) -> list:
   """ 
   
   """
   res = []
   lines = tsv_text.strip().split('\n')
   res = [line.strip().split(delimiter) for line in lines if line.strip()]
   # output
   return res

def df_to_list(df: pd.DataFrame)-> list:
   """
   
   """
   res = []
   for index, row in df.iterrows():
     res.append(row)

   return res

def parse_chebi_json(json_data):
    """
    Parse a ChEBI JSON entry and flatten nested data into a structured DataFrame.

    Parameters
    ----------
    json_data : dict or str
        JSON object or JSON string from the ChEBI REST API.

    Returns
    -------
    pd.DataFrame
        A DataFrame with separate fields for key information such as:
        - Basic info (id, name, definition)
        - Chemical data (mass, formula, charge)
        - Synonyms and names
        - Roles and ontology relations
        - Species origins
    """

    # --- Step 1: Parse string input if necessary ---
    if isinstance(json_data, str):
        json_data = json.loads(json_data)

    # --- Step 2: Extract simple top-level fields ---
    record = {
        "ChEBI_ID": json_data.get("chebi_accession"),
        "Name": json_data.get("name"),
        "ASCII_Name": json_data.get("ascii_name"),
        "Definition": json_data.get("definition"),
        "Stars": json_data.get("stars"),
        "Modified_On": json_data.get("modified_on"),
        "Formula": json_data.get("chemical_data", {}).get("formula"),
        "Mass": json_data.get("chemical_data", {}).get("mass"),
        "Monoisotopic_Mass": json_data.get("chemical_data", {}).get("monoisotopic_mass"),
        "Charge": json_data.get("chemical_data", {}).get("charge"),
        "SMILES": json_data.get("default_structure", {}).get("smiles"),
        "InChI": json_data.get("default_structure", {}).get("standard_inchi"),
        "InChI_Key": json_data.get("default_structure", {}).get("standard_inchi_key"),
    }

    # --- Step 3: Extract names and synonyms ---
    names_dict = json_data.get("names", {})
    record["Synonyms"] = ", ".join(
        [entry["ascii_name"] for entry in names_dict.get("SYNONYM", [])]
    )
    record["IUPAC_Name"] = ", ".join(
        [entry["ascii_name"] for entry in names_dict.get("IUPAC NAME", [])]
    )
    record["UniProt_Name"] = ", ".join(
        [entry["ascii_name"] for entry in names_dict.get("UNIPROT NAME", [])]
    )

    # --- Step 4: Extract species origin info ---
    origins = json_data.get("compound_origins", [])
    record["Species"] = ", ".join(
        [o.get("species_text", "") for o in origins if o.get("species_text")]
    )

    # --- Step 5: Extract ontology relations (incoming/outgoing) ---
    rels = json_data.get("ontology_relations", {})
    incoming = rels.get("incoming_relations", [])
    outgoing = rels.get("outgoing_relations", [])

    record["Incoming_Relations"] = "; ".join(
        [f"{r['init_name']} ({r['relation_type']})" for r in incoming]
    )
    record["Outgoing_Relations"] = "; ".join(
        [f"{r['final_name']} ({r['relation_type']})" for r in outgoing]
    )

    # --- Step 6: Extract role classifications ---
    roles = json_data.get("roles_classification", [])
    record["Roles"] = ", ".join([r["name"] for r in roles])

    # --- Step 7: Convert to pandas DataFrame ---
    df = pd.DataFrame([record])
    return df
