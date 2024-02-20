from Bio import Entrez, SeqIO
from collections import defaultdict
from tqdm import tqdm
import time

class TaxonomyManager:
    def __init__(self, email, api_key="", chunk_size=250):
        self.email = email
        self.api_key = api_key
        self.chunk_size = chunk_size
        Entrez.email = email
        if api_key!="":
        	Entrez.api_key = api_key

    def get_taxid_for_protein(self, protein_ids):
        taxid_map = {}
        assert len(protein_ids) == len(set(protein_ids))

        for i in tqdm(range(0, len(protein_ids), self.chunk_size)):
            chunk_ids = protein_ids[i:i + self.chunk_size]
            handle = Entrez.efetch(db="protein", id=chunk_ids, rettype="gb", retmode="text")
            records = SeqIO.parse(handle, "genbank")

            for r in records:
                protein_id = r.id
                f = r.features[0]
                assert f.type == 'source'
                taxid = f.qualifiers['db_xref'][0].split(':')[1]
                assert taxid.isdigit()
                assert protein_id not in taxid_map
                taxid_map[protein_id] = taxid
        return taxid_map

    def get_lineage(self, taxids):
        lineage_map = {}
        taxid_chunks = [list(taxids.values())[i:i + self.chunk_size] for i in range(0, len(taxids), self.chunk_size)]
        for chunk in tqdm(taxid_chunks):
            id_list = ','.join(chunk)
            handle = Entrez.efetch(db="taxonomy", id=id_list, retmode="xml")
            records = Entrez.read(handle)
            time.sleep(1)
            assert len(records) == len(chunk)

            for record, taxid_ver in zip(records, chunk):
                taxid = record['TaxId']
                assert taxid_ver == taxid

                lineage = record['LineageEx']
                full_lineage = {}
                for item in lineage:
                    full_lineage[item['Rank']] = item['ScientificName']
                lineage_map[taxid] = full_lineage
        return lineage_map

