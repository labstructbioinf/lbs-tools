from lbs.localpdb.utils import multiprocess, os_cmd, mkdir, check_socket_output, get_previous_local_version
from lbs.localpdb.plugins.Plugin import Plugin
import pandas as pd
import os
import pickle





class ECOD(Plugin):

    def __init__(self, pdb):
        self.pdb = pdb

    def load_ecod_data(self, ecod_fn=''):
        tmp_df = pd.read_csv(ecod_fn, sep='\t', skiprows=4, index_col=1)
        tmp_df = tmp_df[tmp_df['chain'] != '.']
        tmp_df['ident'] = tmp_df['pdb'] + '_' + tmp_df['chain']
        del tmp_df['asm_status']
        del tmp_df['manual_rep']
        del tmp_df['#uid']
        cols = tmp_df.columns.tolist()
        tmp_df = tmp_df[[cols[-1]] + [cols[0]] + cols[3:-2] + cols[1:3]]
        self.pdb.ecod = tmp_df

    def update(self):
        pass

    def load(self):
        self.pdb.load_ecod_data = self.load_ecod_data

    def setup(self):
        pass

    def prepare_paths(self):
        pass
