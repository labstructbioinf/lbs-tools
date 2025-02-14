import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from Bio import SeqIO
from lbs.rosetta.utils import parse_score_and_silent
from lbs.utils.multiprocess import multiprocess, os_cmd
import warnings




ff_terms = ['SCORE:', 'score', 'dslf_fa13', 'exph', 'fa_atr',
            'fa_dun_dev', 'fa_dun_rot', 'fa_dun_semi', 'fa_elec',
            'fa_intra_atr_xover4', 'fa_intra_elec', 'fa_intra_rep_xover4',
            'fa_intra_sol_xover4', 'fa_rep', 'fa_sol', 'hbond_bb_sc', 'hbond_lr_bb',
            'hbond_sc', 'hbond_sr_bb', 'hxl_tors', 'lk_ball', 'lk_ball_bridge',
            'lk_ball_bridge_uncpl', 'lk_ball_iso', 'omega', 'p_aa_pp',
            'pro_close', 'rama_prepro', 'ref', 'sc_filter',
            'sc_filter_int_area', 'fa_dun', 'fa_intra_rep', 'rama', 'yhh_planarity']



class RosettaResults(object):

    """Object for processing and manipulating Rosetta results
    Parameters
    ----------
    scorefile : string or list of strings
    Filename(s) of the scorefiles

    silentfile : string or list of strings
    Filename(s) of the corresponding silentfiles

    drop_ff_terms : bool
    Remove additional force field terms


    Attributes
    ----------

    results : Pandas dataframe
    Holds the processed Rosetta results

    """
    def __init__(self, scorefile, silentfile, drop_ff_terms=False):
        """

        :param scorefile: filename of Rosetta scorefile
        :param silentfile: filename of Rosetta silentfiles
        :param drop_ff_terms: determines whether detailed Rosetta information is removed from scorefile
        """

        self.scorefiles = []
        self.silentfiles = []
        self.results = pd.DataFrame()
        if type(scorefile) == str and type(silentfile) == str:
            self.scorefiles.append(scorefile)
            self.silentfiles.append(silentfile)
        elif type(scorefile) == list and type(silentfile) == list:
            self.scorefiles = scorefile
            self.silentfiles = silentfile
        else:
            raise TypeError('Silent and score filename must be both either strings or lists of strings!')
        assert len(self.scorefiles) == len(self.silentfiles)
        
        c = 0
        for scorefile, silentfile in zip(self.scorefiles, self.silentfiles):
            if not (os.path.isfile(scorefile) and os.path.isfile(silentfile)):
                raise FileNotFoundError()
                
            result = parse_score_and_silent(scorefile, silentfile)
            result['source'] = os.path.splitext(scorefile)[0]
            
            # debug
            #print(result)
            
            result.drop('SCORE:', inplace=True, axis=1)
            
            if drop_ff_terms:
                for term in ff_terms:
                    try:
                        result.drop(term, inplace=True, axis=1)
                    except KeyError:
                        pass
            if c == 0:
                self.results = result
            else:
                self.results = pd.concat([self.results, result], ignore_index=True)
            c += 1
            
        seq_lens = {len(seq) for seq in self.results['sequence'].tolist()}
        if len(seq_lens) != 1:
            warnings.warn('Sequences are of different lengths! This may results in an unexpected behavior.')
        model_ids = {model_id for model_id in self.results['description'].tolist()}
        if len(model_ids) != len(self.results['description'].tolist()):
            warnings.warn('Model identificators are not unique. This may results in an unexpected behavior')

    def __repr__(self):
        return("Rosetta results: storing %s models with sequences." % (self.results.shape[0]))

    def __add__(self, other):
        if len(set(self.results.columns) - set(other.results.columns)) != 0:
            raise Error('Cannot combine Rosetta results with different columns!')
        self.results = pd.concat([self.results, other.results])
        return self

    def sort_by_energy(self):
        """
        Sorts results by decoy energy
        :return:
        """
        self.results = self.results.sort_values(by='total_score')

    def get_top_n_percent(self, n):
        """
        Gets n percent of top models from the Rosetta run
        :param n:
        :return:
        """
        self.results = self.results[(self.results.total_score <= self.results.total_score.quantile(n / 100.0))]

    def save_fasta(self, filename):
        """
        Saves fasta file with sequences from Rosetta run
        :param filename: output fasta filename
        :return:
        """
        f = open(filename, 'w')
        for index, row in self.results.iterrows():
            f.write('>%s\n%s\n' % (row['description'], row['sequence']))
            # if len(row['sequence']) != 50:
            # print "A"
        f.close()

    def get_model_list(self, filename):
        """
        Saves the list of model names from Rosetta run
        :param filename: output model list filename
        :return:
        """
        f = open(filename, 'w')
        for model in self.results['description'].tolist():
            f.write('%s\n' % (model))
        f.close()

    def filter_sequences(self, identity, exe='hhfilter', tmp_dir='/tmp/'):
        """
        Filters the sequences from Rosetta run based on their similarity
        :param identity: identity cutoff for filtering
        :param exe: location of the hhfiler binary
        :param tmp_dir: temporary files directory
        :return:
        """
        f_in = tempfile.NamedTemporaryFile(mode='w+b', delete=True)
        fn_in = f_in.name
        f_out = tempfile.NamedTemporaryFile(mode='w+b', delete=True)
        fn_out = f_out.name
        self.save_fasta(fn_in)
        cmd = '%s -id %s -i %s -o %s' % (exe, identity, fn_in, fn_out)
        subprocess.call([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        filtered_sequences = []
        for record in SeqIO.parse(fn_out, "fasta"):
            filtered_sequences.append(record.id)
        for index, row in self.results.iterrows():
            if not row['description'] in filtered_sequences:
                self.results.drop(index, inplace=True)
        f_in.close()
        f_out.close()

    def extract_pdbs(self, out_path='', n_cores=10, exe=''):
        """
        Extracts decoys from Rosetta run in PDB format
        :param out_path: output directory
        :param n_cores: number of cores for parallel processing
        :param exe: location of the extract_pdbs binary
        :return:
        """
        k = np.array(self.results.description.tolist())
        chunks = np.array_split(k, n_cores)
        silent_str = ' '.join(self.silentfiles)
        cmds = []
        for chunk in chunks:
            model_str = ' '.join(chunk)
            cmd = '%s -in:file:silent %s -in:file:tags %s -out:prefix %s' % (exe, silent_str, model_str, out_path)
            cmds.append(cmd)
        multiprocess(os_cmd, cmds, n_cores=n_cores)

    def get_sequence_entropy(self):
        """
        Calculates the entropy of the sequences from Rosetta run
        :return: entropies for each position in the sequence
        """
        canonical_aas = 'ACDEFGHIKLMNPQRSTVWY'
        aa_matrix_column_indices = {}
        aas_ = sorted(list(canonical_aas))
        for x in range(len(aas_)):
            aa_matrix_column_indices[aas_[x]] = x
        sequence_matrix = np.array([list(seq) for seq in self.results['sequence'].tolist()])
        num_sequences, sequence_length = np.shape(sequence_matrix)
        positions = range(sequence_length)
        with np.errstate(divide='ignore', invalid='ignore'):
            count_matrix = np.zeros((sequence_length, len(canonical_aas)))
            for i in positions:
                position_aas = sequence_matrix[:,i]
                unique, counts = np.unique(position_aas, return_counts=True)
                for cx in range(len(unique)):
                        count_matrix[i][aa_matrix_column_indices[unique[cx]]] = counts[cx]
            count_matrix = count_matrix / num_sequences
            count_matrix = -1.0 * count_matrix * (np.log(count_matrix)/np.log(20))
            entropies = np.nansum(count_matrix, axis = 1)
            entropies = dict(zip(positions, entropies.tolist()))
        return entropies