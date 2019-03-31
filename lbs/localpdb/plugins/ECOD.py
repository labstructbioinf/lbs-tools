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
        tmp_df['pdb_id'] = tmp_df['pdb'] + '_' + tmp_df['chain']
        del tmp_df['pdb']
        del tmp_df['chain']
        del tmp_df['asm_status']
        del tmp_df['manual_rep']
        del tmp_df['#uid']
        cols = tmp_df.columns.tolist()
        tmp_df = tmp_df[[cols[-1]] + [cols[0]] + cols[3:-1] + cols[1:3]]
        cols = [x if x != 'f_id' else 'ecod_number' for x in tmp_df.columns]
        tmp_df.columns = cols
        tmp_df['x_id'] = tmp_df['ecod_number'].apply(lambda x: x.split('.')[0])
        tmp_df['h_id'] = tmp_df['ecod_number'].apply(lambda x: x.split('.')[1])
        tmp_df['t_id'] = tmp_df['ecod_number'].apply(lambda x: x.split('.')[2])
        tmp_df['f_id'] = tmp_df['ecod_number'].apply(lambda x: x.split('.')[3] if len(x.split('.')) == 4 else 0)
        #TODO: add this as option...
        tmp_df = tmp_df[tmp_df['pdb_id'].isin(self.pdb.chains.index)]
        self.pdb.ecod = tmp_df

    def update(self):
        pass

    def load(self):
        self.pdb.load_ecod_data = self.load_ecod_data

    def setup(self):
        pass

    def prepare_paths(self):
        pass
