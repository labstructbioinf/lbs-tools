import re
import itertools


def parse_socket_output(filename, method='overlap'):
    """
    Parser of the Socket (J. Mol. Biol., 307 (5), 1427-1450) coiled coil assignment program
    :param filename: Socket output (short version) filename
    :param method: Assignment method
                    'heptads' - based on the heptad assignment
                    'knobs' - based on the knobs assignment
                    'overlap' - sum of the ranges of 'heptads' and 'knobs'
    :return: dict with the parsed coiled coils informations (if present) or 0 (if coiled coils absent)
    """
    if type(filename) is not str:
        return 0
    assert method in ['overlap', 'heptads', 'knobs']

    # Read all lines from file
    try:
        f = open(filename, 'r')
        lines = [line.rstrip() for line in f.readlines()]
        f.close()
    except (OSError, FileNotFoundError):
        return 0

    # Get indices (numbers of first line in file) for each CC assignment
    start_indices = []
    for i in range(0, len(lines)):
        line = lines[i].rstrip()
        if line.startswith('coiled coil') and line.endswith(':'):
            start_indices.append(i)
        if 'Finished' in line:
            start_indices.append(i)

    # Iterate over all CC assignments
    all_coils = {}
    relations = []
    for i in range(0, len(start_indices) - 1):  # Last one is 'Finished' line
        start = start_indices[i]  # Starting line for CC assignment
        stop = start_indices[i + 1]  # End line for CC assignment (next assignment starts or 'Finished')
        cc_id = 'cc_{}'.format(lines[start].split()[2].replace(':', ''))
        # Data
        coil_info = {'helix_ids': [], 'indices': {}, 'sequences': {}, 'heptads': {}, 'ambigous': False, 'relations': []}
        for k in range(start, stop):
            line = lines[k].rstrip()
            if line.startswith('assigning heptad to helix'):
                temp = line.split(' ')
                # Get start and end residues for helix
                try:
                    if temp[6][0] == '-':  # First residue is negative
                        if temp[6].count('-') == 3:  # Both residues are negative
                            start_res = -int(temp[6].split('-')[1])
                            end_res = -int(temp[6].split('-')[3].split(':')[0])
                        else:
                            start_res = -int(temp[6].split('-')[1])
                            end_res = int(temp[6].split('-')[2].split(':')[0])
                    else:
                        start_res = int(temp[6].split('-')[0])
                        end_res = int(temp[6].split('-')[1].split(':')[0])
                    helix_id = 'helix_{}'.format(int(temp[4]))
                    chain = line.split(':')[1]
                except ValueError:
                    coil_info['ambigous'] = True
                coil_info['helix_ids'].append(helix_id)
                seq_line = lines[k + 2].rstrip("\n")
                register_line = lines[k + 3].rstrip("\n")
                knobtype_line = lines[k + 5].rstrip("\n")
                seq = seq_line[9:]
                if len(seq) != (end_res - start_res + 1):
                    coil_info['ambigous'] = True
                register = register_line[9:]
                knobtype = knobtype_line[9:]
                if method == 'heptads':
                    left_marg = len(register) - len(register.lstrip())
                    right_marg = len(register) - len(register.rstrip())
                elif method == 'knobs':
                    left_marg = len(knobtype) - len(knobtype.lstrip('-'))
                    right_marg = len(knobtype) - len(knobtype.rstrip('-'))
                elif method == 'overlap':
                    left_marg1 = len(register) - len(register.lstrip())
                    left_marg2 = len(knobtype) - len(knobtype.lstrip('-'))
                    right_marg1 = len(register) - len(register.rstrip())
                    right_marg2 = len(knobtype) - len(knobtype.rstrip('-'))
                    left_marg = min(left_marg1, left_marg2)
                    right_marg = min(right_marg1, right_marg2)
                start_res = start_res + left_marg
                end_res = end_res - right_marg
                fseq = seq[left_marg:len(seq) - right_marg]
                register = register[left_marg:len(seq) - right_marg].replace(' ', '-')
                coil_info['heptads'][helix_id] = register
                coil_info['sequences'][helix_id] = fseq
                coil_info['indices'][helix_id] = {'start': int(start_res), 'end': int(end_res), 'chain': chain}
            if 'length max' in line and 'PRESENT' not in line and 'REPEATS' not in line:
                inf = re.findall('\((.*?)\)', line)
                data = inf[1].split(' ')
                oligomerization = data[1][0]
                orientation = data[0]
                coil_info['oligomerization'] = oligomerization
                coil_info['orientation'] = orientation
            if line.startswith('	angle between helices'):
                rel_data = line.split()
                first_helix, second_helix, orientation = 'helix_{}'.format(rel_data[3]), 'helix_{}'.format(rel_data[5]), rel_data[8]
                coil_info['relations'].append([first_helix, second_helix, orientation])
        all_coils[cc_id] = coil_info
    if not all_coils:
        return {}
    else:
        return all_coils


def check_socket_output(filename):
    """
    Checks Socket output to determine whether CC domain is present.
    Runs faster than parse_socket_output() and returns only boolean output.
    :param filename: Socket output (short version) filename
    :return bool indicating whether CC domain is present or not
    """
    try:
        f = open(filename, 'rb')
        f.seek(-1024, 2)
        last = f.readlines()[-2].decode()
        f.close()
        if 'NO COILED COILS' in last:
            return False
        else:
            return True
    except (TypeError, KeyError, OSError):
        return False


def __find_index(s, ch, remove_ends=False):
    ind = [i for i, ltr in enumerate(s) if ltr == ch]
    if remove_ends:
        indices = [[match.start(), match.end()] for match in re.finditer('{}+'.format(ch), s)]
        del_ind = set()
        for indice in indices:
            if indice[0] == 0 or indice[1] == len(s):
                for k in range(indice[0], indice[1] + 1):
                    del_ind.add(k)
        ind_corr = [_ind for _ind in ind if _ind not in del_ind]
        return ind_corr
    else:
        return ind


def __change_X_seq(indexes, values, seq):
    seq = list(seq)
    for index, value in zip(indexes, values):
        seq[index] = value
    return ''.join(seq)


def map_socket_to_sequence(seq, chain, socket_output, min_length=7, verbose=False):
    assignment = len(seq)*'0'
    if type(socket_output) != tuple:
        return assignment
    if type(socket_output[0]) != list:
        return assignment
    data = socket_output[0]
    for cc in data:
        for indice, cc_seq in zip(cc['indices'], cc['sequences']):
            if indice[3] == chain:
                matched = False
                if 'X' in cc_seq:
                    if not all([aa == 'X' for aa in cc_seq]):
                        indexes = __find_index(cc_seq, 'X')
                        for perm in itertools.product('MCKHX', repeat=len(indexes)):
                            cc_seq2 = __change_X_seq(indexes, perm, cc_seq)
                            starts = [match.start() for match in re.finditer(re.escape(cc_seq2), seq)]
                            if starts:
                                matched = True
                                if len(cc_seq) >= min_length:
                                    for start in starts:
                                        temp_list = list(assignment)
                                        for i in range(0, len(cc_seq)):
                                            temp_list[start + i] = '1'
                                        assignment = ''.join(temp_list)
                else:
                    starts = [match.start() for match in re.finditer(re.escape(cc_seq), seq)]
                    if starts:
                        matched = True
                        if len(cc_seq) >= min_length:
                            for start in starts:
                                temp_list = list(assignment)
                                for i in range(0, len(cc_seq)):
                                    temp_list[start + i] = '1'
                            assignment = ''.join(temp_list)
                if not matched:
                    if verbose:
                        print(cc_seq, seq, socket_output)
    return assignment
