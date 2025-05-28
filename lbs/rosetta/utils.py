import re
import pandas as pd
import numpy as np

# def silent_to_seq_old(silent_file):  # Kod Staszka do wyciagania fasta z silent file
#     id2seq = {}
#     for l in open(silent_file):
#         if l.startswith("ANNOTATED_SEQUENCE:"):
#             _, seq, id = l.split(' ')
#             seq = seq.replace('[HIS_D]', '')
#             seq = re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", seq)
#             seq = seq.split('[]')
#             final_seq = ''
#             for i in range(0, len(seq)):
#                 if (len(seq[i]) == 1 and len(seq[i + 1]) > 1) or (len(seq[i - 1]) == 1 and len(seq[i]) > 1):
#                     final_seq += seq[i]
#             id2seq[id[:-1]] = {'seq': final_seq}
#     return id2seq

# def silent_to_seq_old(silent_file):
#     id2seq = {}
#     for l in open(silent_file):
#         if l.startswith("ANNOTATED_SEQUENCE:"):
#             _, seq, id = l.split(' ')
#             trail_x = True
#             trail_x_count = 0
#             for c in reversed(seq):
#                 if c == 'X' and trail_x == True:
#                     trail_x_count += 1
#                 else:
#                     trail_x = False
# 
#             par_open = False
#             par_text = []
#             seq_monomeric = []
#             start_indices = []
#             stop_indices = []
#             for i in range(0, len(seq) - trail_x_count):
#                 if par_open:
#                     par_text.append(seq[i])
#                 else:
#                     seq_monomeric.append(seq[i])
#                 try:
#                     if seq[i + 1] == '[':
#                         par_open = True
#                         start_indices.append(i + 1)
#                     elif seq[i] == ']':
#                         par_open = False
#                         stop_indices.append(i)
#                         par_text = []
#                 except IndexError:
#                     pass
#             if trail_x_count == 0:
#                 seq_monomeric = ''.join(seq_monomeric)
#                 id2seq[id[:-1]] = {'seq': seq_monomeric}
#             else:
#                 seq_heteromeric = []
#                 c_indices = []
#                 n_indices = []
#                 del_text = set()
#                 for start, stop in zip(start_indices, stop_indices):
#                     text = seq[start + 1:stop]
#                     full_text = seq[start:stop + 1]
#                     if text.endswith(':NtermProteinFull'):
#                         n_indices.append(start)
#                     elif text.endswith(':CtermProteinFull'):
#                         c_indices.append(start)
#                     else:
#                         del_text.add(full_text)
#                 for start, stop in zip(n_indices, c_indices):
#                     temp_seq = re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", seq[start - 1:stop])
#                     temp_seq = temp_seq.replace(']', '')
#                     temp_seq = temp_seq.replace('[', '')
#                     seq_heteromeric.append(temp_seq)
#                 temp_set = set(seq_heteromeric)
#                 if len(temp_set) == 1:
#                     id2seq[id[:-1]] = {'seq': seq_heteromeric[0]}
#                 else:
#                     id2seq[id[:-1]] = {'seq': seq_heteromeric}
#     return id2seq


def silent_to_seq(silent_file):
    seq_list = []
    
    for l in open(silent_file):
        if l.startswith("ANNOTATED_SEQUENCE:"):
            _, seq, id = l.split(' ')
            trail_x = True
            trail_x_count = 0
            for c in reversed(seq):
                if c == 'X' and trail_x == True:
                    trail_x_count += 1
                else:
                    trail_x = False

            par_open = False
            par_text = []
            seq_monomeric = []
            start_indices = []
            stop_indices = []
            for i in range(0, len(seq) - trail_x_count):
                if par_open:
                    par_text.append(seq[i])
                else:
                    seq_monomeric.append(seq[i])
                try:
                    if seq[i + 1] == '[':
                        par_open = True
                        start_indices.append(i + 1)
                    elif seq[i] == ']':
                        par_open = False
                        stop_indices.append(i)
                        par_text = []
                except IndexError:
                    pass
            if trail_x_count == 0:
                seq_monomeric = ''.join(seq_monomeric)
                seq_list.append(seq_monomeric)
            else:
                seq_heteromeric = []
                c_indices = []
                n_indices = []
                del_text = set()
                for start, stop in zip(start_indices, stop_indices):
                    text = seq[start + 1:stop]
                    full_text = seq[start:stop + 1]
                    if text.endswith(':NtermProteinFull'):
                        n_indices.append(start)
                    elif text.endswith(':CtermProteinFull'):
                        c_indices.append(start)
                    else:
                        del_text.add(full_text)
                for start, stop in zip(n_indices, c_indices):
                    temp_seq = re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", seq[start - 1:stop])
                    temp_seq = temp_seq.replace(']', '')
                    temp_seq = temp_seq.replace('[', '')
                    seq_heteromeric.append(temp_seq)
                temp_set = set(seq_heteromeric)
                if len(temp_set) == 1:
                    seq_list.append(seq_heteromeric[0])
                else:
                    seq_list.append(seq_heteromeric)
                
    return seq_list

# def parse_score_and_silent_OLD(score_file, silent_file):
#     results = pd.read_csv(score_file, sep='\s+', header=0,
#                                 skiprows=1, index_col=False, engine='c', dtype={"total_score": np.float64})
#     id2seq = silent_to_seq(silent_file)
#     sequences = []
#     # Pandas merge
#     for index, row in results.iterrows():
#         try:
#             sequences.append(id2seq[row['description']]['seq'])
#         except KeyError:
#             raise Warning("Model %s present in score file but not present in silent file." % (row['description']))
#             results.drop(index, inplace=True)
#     #Concat sequences and prefixes/suffixes with results dataframe
#     results = pd.concat([results, pd.DataFrame(sequences, columns=['sequence'])], axis=1)
#     return results


def parse_score_and_silent(score_file, silent_file):
    results = pd.read_csv(score_file, sep='\s+', header=0,
                                skiprows=1, index_col=False, engine='c', dtype={"total_score": np.float64})
       

    seq_list = silent_to_seq(silent_file)
    
    assert len(results) == len(seq_list)
    results['sequence'] = seq_list

    return results