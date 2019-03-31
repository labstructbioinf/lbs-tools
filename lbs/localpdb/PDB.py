import pandas as pd
import os
import numpy as np
from Bio import SeqIO
from lbs.localpdb.utils import parse_cluster_data, get_current_local_version, is_nucl_seq, is_x_seq
import importlib


class PDB(object):

    @staticmethod
    def __parse_pdb_data(entries_fn, bundles_fn, res_fn, seqres_fn):
        """
        Parses the downloaded PDB files and returns dataframe with merged data
        :param entries_fn: PDB entries file
        :param bundles_fn: PDB non-standard entries file
        :param res_fn: PDB resolution data
        :return: merged dataframe with all data
        """
        # Parse PDB entries file
        f = open(entries_fn)
        pdb_data = {key: (_type, method) for (key, _type, method) in list(map(str.split, f.readlines()))}
        f.close()
        df_bio = pd.DataFrame.from_dict(pdb_data, orient='index')
        df_bio.columns = ['type', 'method']
        # Parse PDB non-standard entries file
        f = open(bundles_fn)
        pdb_nonstd_data = {ident for ident in map(str.rstrip, f.readlines())}
        f.close()
        df_bio['nonstd'] = df_bio.index.isin(pdb_nonstd_data)
        # Parse PDB resolution data
        f = open(res_fn)
        pdb_res_data = {data[0].lower(): float(data[2]) for data in list(map(str.split, f.readlines()[6:])) if
                        len(data) == 3}
        f.close()
        # Join data into one dataframe
        df_bio = pd.concat([df_bio, pd.DataFrame.from_dict(pdb_res_data, orient='index')], axis=1)
        # Filter non-standard entries that do not have PDB files (large macromolecules)
        df_bio = df_bio[df_bio['nonstd'] == False]
        del df_bio['nonstd']
        df_bio.columns = ['type', 'method', 'resolution']
        df_bio = df_bio[df_bio['type'].isin(['prot', 'prot-nuc'])]
        df_bio = df_bio[df_bio['method'].isin(['diffraction', 'NMR'])]
        # Create dataframe with data in a 'per-chain' format
        pdb_ids = set(df_bio.index.tolist())
        id_seq = {str(entry.id): str(entry.seq) for entry in SeqIO.parse(seqres_fn, 'fasta') if
                  entry.id[0:4] in pdb_ids}
        id_data = {key: [pdb_data[key[0:4]][0], pdb_data[key[0:4]][1], pdb_res_data[key[0:4]]] for key, value in
                   id_seq.items()}
        df_chain = pd.concat(
            [pd.DataFrame.from_dict(id_data, orient='index'), pd.DataFrame.from_dict(id_seq, orient='index')], axis=1,
            sort=True)
        df_chain.columns = ['type', 'method', 'resolution', 'seq']
        df_chain = df_chain[df_chain['type'].isin(['prot', 'prot-nuc'])]
        df_chain = df_chain[df_chain['method'].isin(['diffraction', 'NMR'])]
        # Filter chains with nucleic acids only
        df_chain['is_nucl'] = df_chain['seq'].map(lambda x: is_nucl_seq(x))
        df_chain['is_x'] = df_chain['seq'].map(lambda x: is_x_seq(x))
        df_chain = df_chain[df_chain['is_nucl'] == False]
        df_chain = df_chain[df_chain['is_x'] == False]
        del df_chain['is_nucl']
        del df_chain['is_x']
        # Remove structures (biounits) with nucleic acids only
        ok_ids = {ident.split('_')[0] for ident in df_chain.index.tolist()}
        df_bio = df_bio.loc[ok_ids]
        # Fix resolution of NMR structures
        df_chain['resolution'] = df_chain['resolution'].map(lambda x: x if x > 0 else np.nan)
        df_bio['resolution'] = df_bio['resolution'].map(lambda x: x if x > 0 else np.nan)
        # Return results
        return df_bio, df_chain

    def __init__(self, db_path='', version='latest', plugins=[], verify=False, common_filters=False):
        self.db_path = db_path
        self.abs_db_path = os.path.realpath(self.db_path)
        if version == 'latest':
            self.version = get_current_local_version(db_path)
        else:
            self.version = version
        self._working_path = "{}/data/{}".format(self.db_path, self.version)
        pdb_bundles_fn = '{}/pdb_bundles.txt'.format(self._working_path)
        pdb_entries_fn = "{}/pdb_entries.txt".format(self._working_path)
        pdb_res_fn = "{}/pdb_resolution.txt".format(self._working_path)
        pdb_seqres_fn = "{}/pdb_seqres.txt".format(self._working_path)
        # Create dataframe with original data, "biological units" and per chain data
        self.structures, self.chains = self.__parse_pdb_data(pdb_entries_fn, pdb_bundles_fn, pdb_res_fn,
                                                             pdb_seqres_fn)
        if common_filters:
            self.structures = self.structures[self.structures['resolution'] <= 4.0]
            self.structures = self.structures[self.structures['method'] == 'diffraction']
            self.chains = self.chains[self.chains['resolution'] <= 4.0]
            self.chains = self.chains[self.chains['method'] == 'diffraction']
            self.chains = self.chains[self.chains['seq'].str.len() >= 25]
        if verify:
            self.__verify_integrity()
        # Point to proper file locations
        self.structures['fn_struct'] = self.structures.index.map(
            lambda x: '{}/structures/{}/{}.pdb'.format(self.abs_db_path, x[1:3], x))
        self.structures['fn_bio'] = self.structures.index.map(
            lambda x: '{}/biounits/{}/{}.pdb1'.format(self.abs_db_path, x[1:3], x))
        for plugin in plugins:
            Pl = getattr(importlib.import_module('lbs.localpdb.plugins.{}'.format(plugin)), plugin)
            pl = Pl(self)
            pl.load()

    def __repr__(self):
        return "PDB database (v{}) holding {} entries ({} chains)".format(self.version, len(self.structures),
                                                                          len(self.chains))

    def __verify_integrity(self):
        # Check original PDB files
        for index in self.structures.index.tolist():
            if not os.path.isfile('{}/pdb/{}/{}.pdb'.format(self.db_path, index[1:3], index)):
                print('{}/pdb/{}/{}.pdb'.format(self.db_path, index[1:3], index))
            if not os.path.isfile('{}/pdb_biounit/{}/{}.pdb1'.format(self.db_path, index[1:3], index)):
                print(self.structures.loc[index])

    def __add_col_structures(self, data, added_col_name):
        cols = list(self.structures.columns)
        cols.append(added_col_name)
        if type(data) == dict:
            data = pd.DataFrame.from_dict(data, orient='index')
        self.structures = pd.concat((self.structures, data), axis=1, sort=True)
        self.structures.columns = cols

    def __add_col_chains(self, data, added_col_name):
        cols = list(self.chains.columns)
        cols.append(added_col_name)
        if type(data) == dict:
            data = pd.DataFrame.from_dict(data, orient='index')
        self.chains = pd.concat((self.chains, data), axis=1, sort=True)
        self.chains.columns = cols

    def load_clustering_data(self, identity=50, method='blast'):
        indexes = set(self.chains.index.tolist())
        valid_identities = {'blast': ['30', '40', '50', '70', '90', '95', '100'],
                            'mmseqs': ['30', '40', '50', '70', '90', '95', '100']}
        if str(method) in valid_identities.keys():
            if str(identity) in valid_identities[method]:
                if method == 'blast':
                    fn = '{}/clusters/bc-{}.out'.format(self._working_path, identity)
                    cluster_data = parse_cluster_data(fn)
                    cluster_data_filt = {key: value for key, value in cluster_data.items() if key in indexes}
                else:
                    raise Warning('Clustering method \'{}\' is not implemented yet!'.format(method))
                self.__add_col_chains(cluster_data_filt, '{}-{}'.format(method, identity))
            else:
                raise Warning('Identitity cutoff {} is not available for \'{}\' method!'.format(identity, method))

    def get_new_pdbs(self):
        added_fn = "{}/data/{}/added.pdb".format(self.db_path, self.version)
        if os.path.isfile(added_fn):
            pdbs = set(self.structures.index.tolist())
            f = open(added_fn)
            new_entries = [line.rstrip() for line in f.readlines() if line.rstrip() in pdbs]
            f.close()
            return new_entries
        else:
            return []
