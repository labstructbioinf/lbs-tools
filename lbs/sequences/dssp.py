from Bio.PDB.DSSP import _make_dssp_dict


def get_dssp_seq(fn):
    """
    Extracts sequence from DSSP output file.
    :param fn: input filename
    :return: sequence in PDB structure
    # TODO: Connectivity info??
    """
    f = open(fn, 'r')
    out_dict, keys = _make_dssp_dict(f)
    seq = [(out_dict[key][0]) for key in keys]
    f.close()
    return ''.join(seq)