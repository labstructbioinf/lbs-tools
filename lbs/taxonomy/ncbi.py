import pandas as pd
import numpy as np


'''
	This is an old script for dealing with NCBI taxonomy files
	downloaded from NCBI ftp
'''

class TaxDB:
    def __init__(self, ncbidir):
        self.ncbidir = ncbidir.rstrip('/')  # Remove trailing slash if present

        # Define column names for the DataFrame
        column_names = [
            'taxid', 'parent_taxid', 'rank', 'embl_code', 'divid', 'divflag',
            'gencode', 'gencode_flag', 'mitgencode', 'mitgencode_flag',
            'gb_flag', 'subtree_flag', 'comments'
        ]

        # Read the nodes data file into a pandas DataFrame
        self.nodes = pd.read_csv(
            f"{ncbidir}/nodes.dmp",
            sep='|',
            index_col=False,
            dtype={'taxid': np.int, 'parent_taxid': np.int},
            names=column_names
        )

        # Clean certain columns by removing unwanted characters
        for col in ['embl_code', 'rank', 'comments']:
            self.nodes[col] = self.nodes[col].str.replace(r'[\t\f]', '')

    def keep_genomes(self, genome_taxids):
        """
        Remove species (rank) that aren't in the genome_taxids list.
        """
        self.nodes = self.nodes[
            ~((self.nodes['rank'] == 'species') & (~self.nodes['taxid'].isin(genome_taxids)))
        ]

    def get_children(self, taxid):
        """
        Get all leaf nodes (nodes without children) that are descendants of the given taxid.
        """
        children = self.nodes[self.nodes['parent_taxid'] == taxid]['taxid'].tolist()

        if not children:
            yield taxid
        else:
            for child in children:
                yield from self.get_children(child)
