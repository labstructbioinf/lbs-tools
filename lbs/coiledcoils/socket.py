import re
import pickle
def parse_socket_output(filename, entry, method = 'cc'):
    
    """
    filename: Socket output filename
    entry: PDB entry name
    method: Detection method
        cc - take only the coiled coil assignment
        kih - take the knobs into holes assignment
    """
    ### Read all lines from file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    
    ### Get indices (first line numbers) for each CC assignment
    start_indices = [] 
    for i in range(0, len(lines)):
        line = lines[i].rstrip()
        if line.startswith('coiled coil') and line.endswith(':'):
            start_indices.append(i)
            #print(line)
        if 'Finished' in line:
            start_indices.append(i)
    
    # Iterate over all CC assignments
    all_coils = []
    relations = []
    for i in range(0, len(start_indices)-1): ### Last one is 'Finished' line
        start = start_indices[i] # Starting line for CC assignment
        stop = start_indices[i+1] # End line for CC assignment (next assignment starts or 'Finished')
        
        ### Data
        coil_info = {}
        coil_info['indices'] = []
        coil_info['sequences'] = []
        coil_info['heptads'] = []
        coil_info['orientation'] = []
        for k in range(start, stop):
            line = lines[k].rstrip()
            if entry in line and 'PRESENT' not in line and 'REPEATS' not in line:
                inf = re.findall('\((.*?)\)', line)
                data = inf[1].split(' ')
                oligomerization = data[1][0]
                orientation = data[0]
                coil_info['oligomerization'] = oligomerization
                coil_info['orientation'] =  orientation
            if line.startswith('extent of coiled coil packing:'):
                temp = line.split(':')
                chain = temp[2]
                start_res = temp[1].replace(' ', '').split('-')[0]
                end_res = temp[1].replace(' ', '').split('-')[1]
                line_prev = lines[k-1].rstrip()
                data_prev = line_prev.split(' ')
                helix_id = data_prev[4]
            if line.startswith('sequence'):
                seq2 = line[9:]
                hept = lines[k+1][9:]
                partner = lines[k+2][9:]
                if method == 'cc':
                    ls = len(hept) - len(hept.lstrip())
                    seq2 = seq2[ls:]
                    hept = hept[ls:].rstrip()
                    ts = len(seq2) - len(hept)
                    if ts > 0:
                        seq2 = seq2[:-ts]
                    else:
                        pass
                    assert len(seq2) == len(hept)
                    coil_info['heptads'].append(hept)
                    coil_info['sequences'].append(seq2)
                    coil_info['indices'].append((helix_id, int(start_res), int(end_res), chain))
                elif method == 'kih':
                    st = False
                    partner_start = 0
                    partner_end = 0
                    for z in range(0, len(partner)):
                        if partner[z] != '-' and partner[z] != '\n':
                            #print(partner[z], z)
                            if not st:
                                st = True
                                partner_start = z
                            partner_end = z
                    st = False
                    hept_start = 0
                    hept_end = 0
                    for z in range(0, len(hept)):
                        if hept[z] != ' ' and hept[z] != '\n':
                            if not st:
                                st = True
                                hept_start = z
                            hept_end = z
                    fin_start = 0
                    fin_end = 0
                    if partner_start >= hept_start:
                        fin_start = hept_start
                    else:
                        fin_start = partner_start
                        offset = partner_start - hept_start
                        print(partner_start, hept_start)
                        start_res = int(start_res) + offset
                    if partner_end >= hept_end:
                        fin_end = partner_end
                        offset = partner_end - hept_end
                        end_res = int(end_res) + offset
                    else:
                        fin_end = hept_end
                    coil_info['sequences'].append(seq2[fin_start:fin_end+1])
                    fin_hept = ['-'for y in range(0, len(hept))]
                    coil_info['heptads'].append(''.join(fin_hept))
                    coil_info['indices'].append((helix_id, int(start_res), int(end_res), chain))
            if line.startswith('	angle between helices'):
                rel_data = line.split()
                first_helix, second_helix, orientation = rel_data[3], rel_data[5], rel_data[8]
                relations.append([first_helix, second_helix, orientation])
        all_coils.append(coil_info)
    return all_coils, relations
