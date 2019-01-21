from lbs.localpdb.utils import multiprocess, os_cmd, mkdir
from lbs.localpdb.plugins.Plugin import Plugin
import pandas as pd
from Bio.PDB import *
import os


class PDBMergeBio(Plugin):

    def __init__(self, pdb):
        self.plugin_col_name = 'pdb_merge_bio'
        self.pdb = pdb
        self.required_plugins = []

    def update(self):
        cmds_xray = {}
        cmds_nmr = {}
        for pdb, value in self.pdb.structures.loc[self.pdb.get_new_pdbs()].iterrows():
            in_fn = value['fn_bio']
            out_fn = "{}/pdb_merge_bio/{}/{}.pdb1".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if value['method'] != 'NMR':
                cmd = "timeout 60 pdb_merge_bio -i {} -o {}".format(in_fn, out_fn)
                cmds_xray[pdb] = cmd
            else:
                cmd = [in_fn, out_fn]
                cmds_nmr[pdb] = cmd
        failed = multiprocess(os_cmd, cmds_xray, return_type='failed', print_progress=False)
        failed = multiprocess(self.__nmr_model_to_pdb, cmds_nmr, return_type='failed', print_progress=False)

    def load(self):
        data_dict = {}
        for pdb in self.pdb.structures.index.tolist():
            out_fn = "{}/pdb_merge_bio/{}/{}.pdb1".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if os.path.isfile(out_fn):
                data_dict[pdb] = out_fn
        self.pdb.structures = self.add_col(self.pdb.structures, data_dict, self.plugin_col_name)

    @staticmethod
    def __nmr_model_to_pdb(data):
        in_fn, out_fn = data
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pdb', in_fn)
        try:
            io = PDBIO()
            io.set_structure(structure[0])
            io.save(out_fn)
            return 0
        except KeyError:
            return -1

    def setup(self):
        cmds_xray = {}
        cmds_nmr = {}
        for pdb, value in self.pdb.structures.iterrows():
            in_fn = value['fn_bio']
            out_fn = "{}/pdb_merge_bio/{}/{}.pdb1".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            if value['method'] != 'NMR':
                cmd = "timeout 60 pdb_merge_bio -i {} -o {}".format(in_fn, out_fn)
                cmds_xray[pdb] = cmd
            else:
                cmd = [in_fn, out_fn]
                cmds_nmr[pdb] = cmd

        print("Merging biounits for the x-ray PDB entries...")
        failed = multiprocess(os_cmd, cmds_xray, return_type='failed')
        print("Accounting for NMR entries...")
        failed = multiprocess(self.__nmr_model_to_pdb, cmds_nmr, return_type='failed')

    def prepare_paths(self):
        mkdir("{}/pdb_merge_bio".format(self.pdb.db_path))
        for pdb in self.pdb.structures.index.tolist():
            mkdir("{}/pdb_merge_bio/{}".format(self.pdb.db_path, pdb[1:3]))
