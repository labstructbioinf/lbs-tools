from lbs.localpdb.utils import multiprocess, os_cmd, mkdir, get_previous_local_version
from lbs.localpdb.plugins.Plugin import Plugin
from Bio.PDB import Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import os
import pickle
import numpy as np


class PDBChain(Plugin):

    def __init__(self, pdb):
        self.plugin_col_name = 'pdb_chain'
        self.pdb = pdb
        self.required_plugins = []

    def update(self):
        new_pdbs = self.pdb.get_new_pdbs()
        self.pdb.chains['_new'] = self.pdb.chains.index.map(lambda x: True if x[0:4] in new_pdbs else False)
        tmp = self.pdb.chains[self.pdb.chains['_new'] == True]
        cmds = {}
        for key, value in tmp.iterrows():
            pdb, chain = key.split('_')
            out_fn = "{}/pdb_chain/{}/{}.pdb".format(self.pdb.abs_db_path, pdb[1:3], key)
            in_fn = "{}/structures/{}/{}.pdb".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmds[key] = 'timeout 10 pdb_extract_chain -i {} -o {} -chain {}'.format(in_fn, out_fn, chain)
        failed = multiprocess(os_cmd, cmds, return_type='failed', np=20, print_progress=False)

    def load(self):
        data_dict = {}
        for pdb in self.pdb.chains.index.tolist():
            out_fn = "{}/pdb_chain/{}/{}.pdb".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if os.path.isfile(out_fn):
                data_dict[pdb] = out_fn
        self.pdb.chains = self.add_col(self.pdb.chains, data_dict, self.plugin_col_name)

    def setup(self):
        cmds = {}
        for key, value in self.pdb.chains.iterrows():
            pdb, chain = key.split('_')
            out_fn = "{}/pdb_chain/{}/{}.pdb".format(self.pdb.abs_db_path, pdb[1:3], key)
            in_fn = "{}/structures/{}/{}.pdb".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmds[key] = 'timeout 10 pdb_extract_chain -i {} -o {} -chain {}'.format(in_fn, out_fn, chain)
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)

    def prepare_paths(self):
        mkdir("{}/pdb_chain".format(self.pdb.db_path))
        for pdb in self.pdb.structures.index.tolist():
            mkdir("{}/pdb_chain/{}".format(self.pdb.db_path, pdb[1:3]))
