from lbs.localpdb.utils import multiprocess, os_cmd, mkdir, get_previous_local_version
from lbs.localpdb.plugins.Plugin import Plugin
import os
import pickle
import numpy as np


class PSSM(Plugin):

    def __init__(self, pdb):
        self.plugin_col_name = 'pssm'
        self.pdb = pdb
        self.required_plugins = []
        self.pssm_path = '/home/db/psiblast/PSSM/'
        self.out_path = '/home/db/psiblast/OUT/'
        self.pickle_fn = '{}/pssm/{}.p'.format(self.pdb.db_path, self.pdb.version)
        self.psiblast_db = '/home/users/jludwiczak/nr_vikram/nr90/nr90'
        self.psiblast_exe = '/opt/ncbi-blast-2.7.1+/bin/psiblast'
        self.previous_version = get_previous_local_version(self.pdb.db_path, self.pdb.version)
        self.prev_pickle_fn = '{}/pssm/{}.p'.format(self.pdb.db_path, self.previous_version)

    @staticmethod
    def __pssm_to_seq(fn):
        dat = np.genfromtxt(fn, skip_header=3, skip_footer=5, usecols=(1), dtype=(str))
        return ''.join(dat)

    @staticmethod
    def __seq_to_fn(seq, dat):
        try:
            fn = dat[seq]
        except KeyError:
            fn = np.nan
        return fn

    def update(self):
        new_pdbs = self.pdb.get_new_pdbs()
        data = pickle.load(open(self.prev_pickle_fn, 'rb'))
        self.pdb.chains[self.plugin_col_name] = self.pdb.chains['seq'].map(lambda x: self.__seq_to_fn(x, data))
        self.pdb.chains['_new'] = self.pdb.chains.index.map(lambda x: True if x[0:4] in new_pdbs else False)
        tmp = self.pdb.chains[(self.pdb.chains['_new'] == True) & (self.pdb.chains['pssm'].isnull())]
        data = {}
        for key, value in tmp.iterrows():
            if value['seq'] in data.keys():
                data[value['seq']].append(key)
            else:
                data[value['seq']] = [key]
        cmds = {}
        for key, value in data.items():
            ident = value[0]
            in_fn = "{}/pssm/.tmp/{}.fasta".format(self.pdb.db_path, ident)
            f = open(in_fn, 'w')
            f.write('>{}\n{}\n'.format(ident, key))
            f.close()
            cmd = '{} -query {} -db {} -evalue 0.001 -num_iterations 3 -num_threads 8 -out_ascii_pssm {}/{}.pssm -outfmt 5 -out {}/{}.out'.format(
                self.psiblast_exe, in_fn, self.psiblast_db, self.pssm_path, ident, self.out_path, ident)
            cmds[ident] = cmd
        multiprocess(os_cmd, cmds, np=1)
        self.setup()

    def load(self):
        data = pickle.load(open(self.pickle_fn, 'rb'))
        self.pdb.chains[self.plugin_col_name] = self.pdb.chains['seq'].map(lambda x: self.__seq_to_fn(x, data))

    def setup(self):
        fns = os.listdir(self.pssm_path)
        data = {}
        for fn in fns:
            fn = '/home/db/psiblast/PSSM/{}'.format(fn)
            seq = self.__pssm_to_seq(fn)
            data[seq] = fn
        pickle.dump(data, open(self.pickle_fn, 'wb'))

    def prepare_paths(self):
        mkdir("{}/pssm".format(self.pdb.db_path))
        mkdir("{}/pssm/.tmp".format(self.pdb.db_path))
