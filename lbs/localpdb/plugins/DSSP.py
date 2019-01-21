from lbs.localpdb.utils import multiprocess, os_cmd, mkdir
from lbs.localpdb.plugins.Plugin import Plugin
import pandas as pd
from Bio.PDB import *
import os


class DSSP(Plugin):

    def __init__(self, pdb):
        self.pdb = pdb
        self.req_columns = {'chains': ['pdb_chain'], 'structures': ['pdb_merge_bio']}

    def load(self):
        data_dict = {}
        for pdb in self.pdb.structures.index.tolist():
            out_fn = "{}/dssp/{}/{}.dssp".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if os.path.isfile(out_fn):
                data_dict[pdb] = out_fn
        self.pdb.structures = self.add_col(self.pdb.structures, data_dict, 'dssp')
        data_dict = {}
        for pdb in self.pdb.chains.index.tolist():
            out_fn = "{}/dssp/{}/{}.dssp".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if os.path.isfile(out_fn):
                data_dict[pdb] = out_fn
        self.pdb.chains = self.add_col(self.pdb.chains, data_dict, 'dssp')

    def update(self):
        cmds = {}
        new_pdbs = self.pdb.get_new_pdbs()

        update_df = self.pdb.structures.loc[new_pdbs]
        for pdb, value in update_df.iterrows():
            in_fn = value['pdb_merge_bio']
            out_fn = "{}/dssp/{}/{}.dssp".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "dssp -i {} -o {}".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)

        self.pdb.chains['_new'] = self.pdb.chains.index.map(lambda x: True if x[0:4] in new_pdbs else False)
        tmp = self.pdb.chains[self.pdb.chains['_new'] == True]
        cmds = {}
        for pdb, value in tmp.iterrows():
            in_fn = value['pdb_chain']
            out_fn = "{}/dssp/{}/{}.dssp".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "dssp -i {} -o {}".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)

    def setup(self):
        cmds = {}
        for pdb, value in self.pdb.structures.iterrows():
            in_fn = value['pdb_merge_bio']
            out_fn = "{}/dssp_biounit/{}/{}.dssp".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "dssp -i {} -o {}".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)
        cmds = {}
        for pdb, value in self.pdb.chains.iterrows():
            in_fn = value['pdb_chain']
            out_fn = "{}/dssp_chain/{}/{}.dssp".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "dssp -i {} -o {}".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)

    def prepare_paths(self):
        mkdir("{}/dssp".format(self.pdb.db_path))
        for pdb in self.pdb.chains.index.tolist():
            mkdir("{}/dssp/{}".format(self.pdb.db_path, pdb[1:3]))
