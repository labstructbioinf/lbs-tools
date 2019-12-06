import numpy as np

radius = {
    'A': 88.6,
    'R': 173.4,
    'G': 60.1,
    'S': 89.0,
    'C': 108.5,
    'D': 111.1,
    'P': 112.7,
    'N': 114.1,
    'T': 116.1,
    'E': 138.4,
    'Q': 143.8,
    'H': 153.2,
    'M': 162.9,
    'I': 166.7,
    'L': 166.7,
    'K': 168.6,
    'F': 189.9,
    'Y': 193.6,
    'W': 227.8,
    'V': 140.0
}
phobos = {
    'I': 4.5,
    'V': 4.2,
    'L': 3.8,
    'F': 2.8,
    'C': 2.5,
    'M': 1.9,
    'A': 1.8,
    'G': -0.4,
    'T': -0.7,
    'S': -0.8,
    'W': -0.9,
    'Y': -1.3,
    'P': 1.6,
    'H': -3.2,
    'D': -3.5,
    'E': -3.5,
    'N': -3.5,
    'Q': -3.5,
    'K': -3.9,
    'R': -4.5
}


def enc_seq_onehot(seq, pad_length=None, pad_left=0):
    """
    Encodes amino acid sequence with one-hot encoding scheme
    :param seq: amino acid sequence
    :param pad_length: zero padding length (assures equal length of variable length sequences)
    :return: encoded sequence
    """
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa_indices = {aa1[k]: k for k in range(0, len(aa1))}
    enc_seq = []
    for aa in seq:
        enc_aa = np.zeros(20)
        enc_aa[aa_indices[aa]] = 1
        enc_seq.append(enc_aa)
    matrix = np.asarray(enc_seq)
    if pad_length:
        pad_matrix = np.zeros((pad_length, 20))
        pad_matrix[pad_left:matrix.shape[0]+pad_left, 0:matrix.shape[1]] = matrix
        return pad_matrix
    return matrix


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


phobos_norm = {key: 4*(value-min(phobos.values()))/(max(phobos.values())-min(phobos.values())) for key, value in phobos.items()}
radius_norm = {key: 4*(value-min(radius.values()))/(max(radius.values())-min(radius.values())) for key, value in radius.items()}


def enc_seq_sizehydro(seq, pad_length=None, pad_left=0):
    """
    Encodes amino acid sequence based on the size and hydrophobicity
    :param seq: amino acid sequence
    :param pad_length: zero padding length (assures equal length of variable length sequences)
    :return: encoded sequence
    """
    enc_seq = []
    for aa in seq:
        enc_aa_size = np.zeros(5)
        enc_aa_hydro = np.zeros(5)
        enc_aa_size[int(radius_norm[aa])] = 1
        enc_aa_hydro[int(phobos_norm[aa])] = 1
        enc_aa = np.concatenate((enc_aa_size, enc_aa_hydro))
        enc_seq.append(enc_aa)
    matrix = np.asarray(enc_seq)
    if pad_length:
        pad_matrix = np.zeros((pad_length, 10))
        pad_matrix[pad_left:matrix.shape[0]+pad_left, 0:matrix.shape[1]] = matrix
        return pad_matrix
    return matrix


def enc_pssm(pssm_file, pad_length=None, pad_left=0, skip_footer=5):
    """
    Encodes evolutionary information stored in PSSM file
    :param pssm_file: PSSM filename
    :param pad_length: zero padding length (assures equal length of variable length sequences)
    :param skip_footer: PSSM footer may have different lengths
    :return: encoded PSSM
    """
    pssm_matrix = sigmoid(np.genfromtxt(pssm_file, skip_header=3, skip_footer=skip_footer, usecols=(i for i in range(2, 22))))
    if pad_length:
        pad_matrix = np.zeros((pad_length, 20))
        pad_matrix[pad_left:pssm_matrix.shape[0]+pad_left, 0:pssm_matrix.shape[1]] = pssm_matrix
        return pad_matrix
    return pssm_matrix


def enc_pipred(pipred_file, pad_length=None, pad_left=0):
    """
    Encodes results of the PiPred secondary structure prediction
    :param pipred_file: PiPred results filename
    :param pad_length: zero padding length (assures equal length of variable length sequences)
    :return: encoded PiPred results
    """
    pipred_matrix = np.genfromtxt(pipred_file, usecols=(1, 2, 3, 4))
    if pad_length:
        pad_matrix = np.zeros((pad_length, 4))
        pad_matrix[pad_left:pipred_matrix.shape[0]+pad_left, 0:pipred_matrix.shape[1]] = pipred_matrix
        return pad_matrix
    return pipred_matrix


def enc_label(labels, alphabet='01', pad_length=None, pad_left=0):
    """
    Encodes labels based on the one-hot scheme based on the provided alphabet
    :param labels: string with labels
    :param alphabet: alphabet used for labeling
    :param pad_length: zero padding length (assures equal length of variable length labels)
    :return:
    """
    label1 = list(alphabet)
    label_indices = {label1[k]: k for k in range(0, len(label1))}
    enc_labels = []
    for label in labels:
        enc_label = np.zeros(len(label1))
        enc_label[label_indices[label]] = 1
        enc_labels.append(enc_label)
    matrix = np.asarray(enc_labels)
    if pad_length:
        pad_matrix = np.zeros((pad_length, len(label1)))
        pad_matrix[pad_left:matrix.shape[0]+pad_left, 0:matrix.shape[1]] = matrix
        return pad_matrix
    return matrix
