from lbs.localpdb.utils import multiprocess, os_cmd, mkdir, check_socket_output, get_previous_local_version
from lbs.localpdb.plugins.Plugin import Plugin
from lbs.coiledcoils.socket import parse_socket_output
import numpy as np
import os
import pickle


class Socket(Plugin):

    def __init__(self, pdb):
        self.pdb = pdb
        self.pickle_fn = '{}/socket/{}.p'.format(self.pdb.db_path, self.pdb.version)

    def update(self):
        previous_version = get_previous_local_version(self.pdb.db_path, self.pdb.version)
        prev_pickle_fn = '{}/socket/{}.p'.format(self.pdb.db_path, previous_version)
        cmds = {}
        update_df = self.pdb.structures.loc[self.pdb.get_new_pdbs()]
        for pdb, value in update_df.iterrows():
            in_dssp = value['dssp']
            in_fn = value['pdb_merge_bio']
            out_70 = "{}/socket/{}/{}.cc_70".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd_70 = "/opt/apps/socket3.03/socket -f {} -s {} -c 7.0 > {}".format(in_fn, in_dssp, out_70)
            out_72 = "{}/socket/{}/{}.cc_72".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd_72 = "/opt/apps/socket3.03/socket -f {} -s {} -c 7.2 > {}".format(in_fn, in_dssp, out_72)
            out_74 = "{}/socket/{}/{}.cc_74".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd_74 = "/opt/apps/socket3.03/socket -f {} -s {} -c 7.4 > {}".format(in_fn, in_dssp, out_74)
            cmds['{}_70'.format(pdb)] = cmd_70
            cmds['{}_72'.format(pdb)] = cmd_72
            cmds['{}_74'.format(pdb)] = cmd_74
        failed = multiprocess(os_cmd, cmds, return_type='failed', print_progress=False)
        up = {1: 'python3.5 {}/update_digest.py'.format(os.path.dirname(os.path.realpath(__file__)))}
        failed = multiprocess(os_cmd, up, return_type='failed', print_progress=False)

        data_dict = pickle.load(open(prev_pickle_fn, 'rb'))
        for pdb in update_df.index.tolist():
            fns = {}
            for variant in ['70', '72', '74']:
                fn = "{}/socket/{}/{}.cc_{}".format(self.pdb.abs_db_path, pdb[1:3], pdb, variant)
                cc = check_socket_output(fn)
                if cc:
                    fns[variant] = fn
                else:
                    fns[variant] = 0
            data_dict[pdb] = fns
        pickle.dump(data_dict, open(self.pickle_fn, 'wb'))

    def load(self):
        data = pickle.load(open(self.pickle_fn, 'rb'))
        pdbs = set(self.pdb.structures.index.tolist())
        data_clean = {key: value for key, value in data.items() if key in pdbs}
        self.pdb.structures = self.add_col(self.pdb.structures, data_clean, 'socket')

    def setup(self):
        cmds = {}
        for pdb, value in self.pdb.structures[~self.pdb.structures['dssp'].isnull()].iterrows():
            in_dssp = value['dssp']
            in_fn = value['pdb_merge_bio']
            out_70 = "{}/socket/{}/{}.cc_70".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd_70 = "/opt/apps/socket3.03/socket -f {} -s {} -c 7.0 > {}".format(in_fn, in_dssp, out_70)
            out_72 = "{}/socket/{}/{}.cc_72".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd_72 = "/opt/apps/socket3.03/socket -f {} -s {} -c 7.2 > {}".format(in_fn, in_dssp, out_72)
            out_74 = "{}/socket/{}/{}.cc_74".format(self.pdb.abs_db_path, pdb[1:3], pdb)
            cmd_74 = "/opt/apps/socket3.03/socket -f {} -s {} -c 7.4 > {}".format(in_fn, in_dssp, out_74)
            cmds['{}_70'.format(pdb)] = cmd_70
            cmds['{}_72'.format(pdb)] = cmd_72
            cmds['{}_74'.format(pdb)] = cmd_74
        failed = multiprocess(os_cmd, cmds, return_type='failed')
        data_dict = {}
        for pdb in self.pdb.structures.index.tolist():
            fns = {}
            for variant in ['70', '72', '74']:
                fn = "{}/socket/{}/{}.cc_{}".format(self.pdb.abs_db_path, pdb[1:3], pdb, variant)
                cc = check_socket_output(fn)
                if cc:
                    fns[variant] = fn
                else:
                    fns[variant] = 0
            data_dict[pdb] = fns
        pickle.dump(data_dict, open(self.pickle_fn, 'wb'))

    def prepare_paths(self):
        mkdir("{}/socket".format(self.pdb.db_path))
        for pdb in self.pdb.structures.index.tolist():
            mkdir("{}/socket/{}".format(self.pdb.db_path, pdb[1:3]))
