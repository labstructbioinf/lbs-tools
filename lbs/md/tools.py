import os
from functools import partial
from typing import Union, List, Tuple

import atomium
import numpy as np
from scipy.spatial.distance import cdist



def get_atom_xyz(residue, atom_name):
    for a in residue.atoms():
        if a.name == atom_name:
            return a.location
    return (np.nan, np.nan, np.nan)

get_CA_xyz = partial(get_atom_xyz, atom_name='CA')

def parse_xyz(path: Union[str, atomium.structures.Model],
              chain: Union[str, None] = None,
              get_pdb_ss:bool=False)->(np.ndarray, List[str]):
    '''
    params:
        chain (str, None) default None if str read chain if None read all
    return:
        ca_xyz (torch.Tensor), sequence (list)
        
    '''
    if isinstance(path, str):
        if not os.path.isfile(path):
            raise KeyError(f'path: {path} doesn\'t exist')
        else:
            data = atomium.open(path).model
            chain = data.chain(chain) if chain is not None else data
    elif isinstance(path, atomium.structures.Model):
        chain = path.chain(chain) if chain is not None else path
    else:
        raise KeyError(f'invalid path arg type {path}')
    sequence = list(map(lambda x: x.code, chain.residues()))
    ca_xyz = list(map(get_CA_xyz, chain.residues()))
    ca_xyz = np.array(ca_xyz, dtype=np.float32)
    if get_pdb_ss:
        secondary = list(map(get_ss_label, chain.residues()))
        return ca_xyz, sequence, secondary
    else:
        return ca_xyz, sequence
    
    

def calculate_box_size(path: str, factor: int = 1.5) -> Tuple[float]:
    '''
    calculate box size for given structure
    factor defines box width with formula 
     max_CA_dist * factor = box width
    params:
        path (str) pdb structure
        factor (int)
    return
        box_size (tuple[float])
    '''
    xyz, _ = parse_xyz(path)
    dist = cdist(xyz, xyz)
    max_dist = dist.max()/10 # scale to nanometers
    box_wall = max_dist*factor
    box_wall = np.clip(box_wall, a_min=5, a_max=np.nan).astype(np.float16)
    box_size = (box_wall, box_wall, box_wall)
    return box_size