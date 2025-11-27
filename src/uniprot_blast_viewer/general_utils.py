import sys
import json
import requests

import pandas as pd
import re

import json

from textwrap import wrap

# Helper function to download data
#def get_url(url, **kwargs):
#  response = requests.get(url, **kwargs);

#  if not response.ok:
#    print(response.text)
#    response.raise_for_status()
#    #sys.exit()

#  return response

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