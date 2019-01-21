from lbs.localpdb.utils import multiprocess, os_cmd, mkdir
from lbs.localpdb.plugins.Plugin import Plugin
import pandas as pd
from Bio.PDB import *
import os


class Master(Plugin):
    def __init__(self, pdb):
        self.pdb = pdb
        self.req_columns = {'chains': ['pdb_chain'], 'structures': ['pdb_merge_bio']}

    def load(self):
        data_dict = {}
        for pdb in self.pdb.structures.index.tolist():
            out_fn = "{}/master/{}/{}.pds".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if os.path.isfile(out_fn):
                data_dict[pdb] = out_fn
        self.pdb.structures = self.add_col(self.pdb.structures, data_dict, 'master')

        data_dict = {}
        for pdb in self.pdb.chains.index.tolist():
            out_fn = "{}/master/{}/{}.pds".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if os.path.isfile(out_fn):
                data_dict[pdb] = out_fn
        self.pdb.chains = self.add_col(self.pdb.chains, data_dict, 'master')

    def update(self):
        new_pdbs = self.pdb.get_new_pdbs()

        cmds = {}
        update_df = self.pdb.structures.loc[new_pdbs]
        for pdb, value in update_df.iterrows():
            in_fn = value['pdb_merge_bio']
            out_fn = "{}/master/{}/{}.pds".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "timeout 60 /opt/apps/master-1.3.1/createPDS --pdb {} --pds {} --type target".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)

        self.pdb.chains['_new'] = self.pdb.chains.index.map(lambda x: True if x[0:4] in new_pdbs else False)
        tmp = self.pdb.chains[self.pdb.chains['_new'] == True]
        cmds = {}
        for pdb, value in tmp.iterrows():
            in_fn = value['pdb_chain']
            out_fn = "{}/master/{}/{}.pds".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "timeout 60 /opt/apps/master-1.3.1/createPDS --pdb {} --pds {} --type target".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)

    def setup(self):
        cmds = {}
        for pdb, value in self.pdb.chains.iterrows():
            in_fn = value['pdb_chain']
            out_fn = "{}/master/{}/{}.pds".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "timeout 60 /opt/apps/master-1.3.1/createPDS --pdb {} --pds {} --type target".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed')

        cmds = {}
        for pdb, value in self.pdb.structures.iterrows():
            in_fn = value['pdb_merge_bio']
            out_fn = "{}/master/{}/{}.pds".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd = "timeout 60 /opt/apps/master-1.3.1/createPDS --pdb {} --pds {} --type target".format(in_fn, out_fn)
            cmds[pdb] = cmd
        failed = multiprocess(os_cmd, cmds, return_type='failed')

    def prepare_paths(self):
        mkdir("{}/master".format(self.pdb.db_path))
        for pdb in self.pdb.structures.index.tolist():
            mkdir("{}/master/{}".format(self.pdb.db_path, pdb[1:3]))