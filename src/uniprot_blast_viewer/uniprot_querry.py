import sys
import json
import requests
from prettytable import PrettyTable

import pandas as pd
import re
import time

from textwrap import wrap

from .globals import WEBSITE_API, PROTEINS_API, EBI_API
from .general_utils import get_url, get_url_blast, list_to_dataframe, sequence_to_fasta, save_tsv_to_csv, tsv_to_file, tsv_to_list, df_to_list, parse_chebi_json

def getFasta_from_accession(accession: str)-> str:
   """
   
   """
   # get FASTA file
   r = get_url(f"{WEBSITE_API}/uniprotkb/{accession}?format=fasta")
   return r.text


def doBlast_from_accession(accession: str)-> str:
   """
   
   """

   fasta = getFasta_from_accession(accession)
   # submit blast job
   r = requests.post("https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run", data={
                     "email": "example@example.com",
                     "program": "blastp",
                     "matrix": "BLOSUM62",
                     "alignments": 250,
                     "scores": 250,
                     "exp": 10,
                     "filter": "F",
                     "gapalign": "true",
                     "stype": "protein",
                     "database": "uniprotkb_refprotswissprot",
                    "sequence": fasta
    })
    # documentation here https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?pageId=94147939#NCBIBLAST+HelpandDocumentation-RESTAPI

   job_id = r.text
   return job_id

def doBlast_from_sequence(sequence: str)-> str:
   """
   
   """

   fasta = sequence_to_fasta(sequence)
   # submit blast job
   r = requests.post("https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run", data={
                     "email": "example@example.com",
                     "program": "blastp",
                     "matrix": "BLOSUM62",
                     "alignments": 250,
                     "scores": 250,
                     "exp": 10,
                     "filter": "F",
                     "gapalign": "true",
                     "stype": "protein",
                     "database": "uniprotkb_refprotswissprot",
                    "sequence": fasta
    })
    # documentation here https://www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?pageId=94147939#NCBIBLAST+HelpandDocumentation-RESTAPI

   job_id = r.text
   return job_id


def see_Blast_status(jobid: str)-> str:
   """
   
   """
   # get job status
   r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{jobid}")
   return r.text

def see_Blast_result(jobid: str)-> str:
   """
   """
   # Run the following again to check the status until finished
   #job_id = "ncbiblast-R20250925-172248-0014-13953597-p1m"

   # get job status
   #r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{jobid}")
   #print(r.text)

   r = get_url_blast(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{jobid}/out")
   return r

def see_Blast_result_json(jobid: str)-> str:
   """
   """
   # Run the following again to check the status until finished
   #job_id = "ncbiblast-R20250925-172248-0014-13953597-p1m"

   # get job status
   #r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{jobid}")
   #print(r.text)

   #r = get_url_blast(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{jobid}/json")
   return(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{jobid}/json")

# Download BLAST results as TSV
# ---------------------------------------------------------
def download_results(job_id):
    """
    Download results in TSV format.

    API endpoint:
        https://rest.uniprot.org/blast/results/{job_id}?format=tsv
    """

    result_url = f"https://rest.uniprot.org/blast/results/{job_id}"

    params = {"format": "tsv"}

    r = requests.get(result_url, params=params)

    print("[INFO] Results downloaded successfully.")
    return r.text

def do_Blast(accession: str)-> str:
    """

    """

    r = doBlast_from_accession(accession)

    status = 'NO'
    while status != 'FINISHED':
        time.sleep(0.5)
        status = see_Blast_status(r)

    res = see_Blast_result(r)
    return res

def do_seq_Blast_get_json(sequence: str)-> str:
    """

    """

    r = doBlast_from_sequence(sequence)

    status = 'NO'
    while status != 'FINISHED':
        time.sleep(0.5)
        print(r)
        status = see_Blast_status(r)

    json_url = see_Blast_result_json(r)
    return json_url

def get_isoforms(accession: str, csv_file: str):
   """
   
   """
   # natural variants info for O60260 / PRKN_HUMAN
   r = get_url(f"{WEBSITE_API}/uniprotkb/search?query={accession}&&includeIsoform=true&fields=accession,cc_function,cc_subcellular_location,cc_ptm,sequence&format=tsv")

   #print(r.text)
   tsv_to_file(r.text, csv_file)

def get_isoforms_list(accession: str) -> list:
   """
   
   """
   # natural variants info for O60260 / PRKN_HUMAN
   r = get_url(f"{WEBSITE_API}/uniprotkb/search?query={accession}&&includeIsoform=true&fields=accession,cc_function,cc_subcellular_location,cc_ptm,sequence&format=tsv")

   #print(r.text)
   return tsv_to_list(r.text)


def align_all_accessions(accessions: list)-> str:
   """
   
   """
   # manually selected accessions
   str_accessions = ",".join(accessions)

   r = get_url(f"{WEBSITE_API}/uniprotkb/accessions?accessions={str_accessions}&format=fasta")
   fasta = r.text
   #print(fasta)

   # submit align job using clustalo
   r = requests.post("https://www.ebi.ac.uk/Tools/services/rest/clustalo/run", data={
                     "email": "example@example.com",
                     "iterations": 0,
                     "outfmt": "clustal_num",
                     "order": "aligned",
                     "sequence": fasta
                    })

   # documentation here https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation#ClustalOmegaHelpandDocumentation-RESTAPI
   job_id = r.text


def see_alignJob_status(jobid: str)->str:
   """
   
   """
   # get job status
   r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{jobid}")
   return r.text

def see_alignJob_result(jobid: str)-> str:
   """
   
   """
   # * : Fully conserved residues.
   # : : Conservation between groups of strongly similar properties (Gonnet PAM 250 score > 0.5).
   # . : Conservation between groups of weakly similar properties (Gonnet PAM 250 score ≤ 0.5).
   #   : Non-conserved residues.
   r = get_url(f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{jobid}/aln-clustal_num")
   return r.text

def get_aligned_results_accessions(accessions: list)-> str:
   """
   
   """

   r = align_all_accessions(accessions)

   status = 'NO'
   while status != 'FINISHED':
     time.sleep(0.5)
     status = see_alignJob_status(r)

     #res = see_alignJob_result('clustalo-R20251008-151813-0476-53561636-p1m')

   res = see_alignJob_result(r)     
   return res

def get_PTMs_from_accession(accession: str, csv_file: str):
   """
   
   """

   # documentation https://www.ebi.ac.uk/proteins/api/doc/#!/proteomics-ptm/getByAccession
   #r = get_url(f"{PROTEINS_API}/proteomics-ptm/{accession}")
   r = get_url(f"{PROTEINS_API}/proteomics/ptm/{accession}")

   if r is None:
       return None

   data = r.json()

   t = PrettyTable(['name', 'position', 'sources', 'id', 'confidence'])
   table_data = []
   for feature in data['features']:
       for ptm in feature['ptms']:
           for dbRef in ptm['dbReferences']:
                table_data.append(
                                  (ptm['name'],
                                   int(feature['begin']) + int(ptm['position']) - 1,','.join(ptm['sources']),
                                   dbRef['id'],
                                   dbRef['properties']['Confidence score']))

   # sort by "position" column
   table_data = sorted(table_data, key=lambda x: x[1])
   print(table_data)
   # filter by "confidence" column, only "Gold" values
   #table_data = filter(lambda x: x[4] == 'Gold', table_data)
   #print(table_data)

   #
   df_table = list_to_dataframe(table_data, ['name', 'position', 'sources', 'id', 'confidence'])
   df_table.to_csv(csv_file)

   t.add_rows(table_data)

def get_PTMs_from_accession_list(accession: str) -> list:
   """
   
   """

   # documentation https://www.ebi.ac.uk/proteins/api/doc/#!/proteomics-ptm/getByAccession
   #r = get_url(f"{PROTEINS_API}/proteomics-ptm/{accession}")
   r = get_url(f"{PROTEINS_API}/proteomics/ptm/{accession}")

   if r is None:
       return None

   data = r.json()

   t = PrettyTable(['name', 'position', 'sources', 'id', 'confidence'])
   table_data = []
   for feature in data['features']:
       for ptm in feature['ptms']:
           for dbRef in ptm['dbReferences']:
                table_data.append(
                                  [ptm['name'],
                                   int(feature['begin']) + int(ptm['position']) - 1,','.join(ptm['sources']),
                                   dbRef['id'],
                                   dbRef['properties']['Confidence score']])

   # sort by "position" column
   table_data = sorted(table_data, key=lambda x: x[1])
   table_data.insert(0,['name', 'position', 'sources', 'id', 'confidence'])
   # filter by "confidence" column, only "Gold" values
   #table_data = filter(lambda x: x[4] == 'Gold', table_data)
   #print(table_data)

   #
   df_table = list_to_dataframe(table_data, ['name', 'position', 'sources', 'id', 'confidence'])
   t.add_rows(table_data)

   return table_data


def import_BindingInfo_uniProt(accession: str)-> list:
    """

    """
    r = get_url(f"{WEBSITE_API}/uniprotkb/{accession}?fields=ft_binding,sequence&format=json")
    json_input = r.text
    
    # --- Step 1: Parse JSON string if needed ---
    if isinstance(json_input, str):
        try:
            data = json.loads(json_input)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON input: {e}")
    else:
        data = json_input  # Already a dict or list

    features = data['features']

    #df = pd.DataFrame(columns=['type', 
    #                           'start_pos','start_modifier', 
    #                           'end_pos','end_modifier',
    #                           'ligand_name','ligand_id',
    #                           'mass', 'monisotopic_mass', 'charge'])
    
    res_table = [['type', 
                  'start_pos','start_modifier', 
                  'end_pos','end_modifier',
                  'ligand_name','ligand_id',
                  'mass', 'monisotopic_mass', 'charge']]

    for items in features:
        type      = items['type']
        
        start_val = items['location']['start']['value']
        start_mod = items['location']['start']['modifier']

        end_val = items['location']['end']['value']
        end_mod = items['location']['end']['modifier']
        
        lig_name = items['ligand']['name']
        lig_id   = items['ligand']['id']

        if lig_id.startswith('ChEBI'):
            parts = lig_id.split(':')
            res = get_mass_chEBI(f"{parts[1]}:{parts[2]}")

            #print(res['Monoisotopic_Mass'])

            mass = res['Mass']
            monoisotopic_mass = res['Monoisotopic_Mass']
            charge = res['Charge']
        else:
            mass = ''
            monoisotopic_mass = '' 
            charge = '' 
        
        tmp = [type, 
               start_val, start_mod,
               end_val,   end_mod,
               lig_name,  lig_id,
               mass, monoisotopic_mass, charge]
        
        res_table.append(tmp)
       

    return res_table

#querry chEBI database
def get_mass_chEBI(chebi_accno: str)-> dict:
    """
    
    """
    r = get_url(f"{EBI_API}/compound/{chebi_accno}/?only_ontology_parents=false&only_ontology_children=false")

    df = parse_chebi_json(r.text)

    new_dict = {"Mass": float(df['Mass'].values[0]), 
                "Monoisotopic_Mass": float(df['Monoisotopic_Mass'].values[0]), 
                "Charge": int(df['Charge'].values[0])}
    return new_dict