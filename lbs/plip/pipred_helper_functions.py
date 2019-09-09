###!!!###

### EVERYTHING HERE NEEDS REVISION IF STILL USEFUL ###

###!!!###


### functions

## for each seq/ss pair in dict

# check if pi-helix present - yes: pass, no: continue

# get all pi-helix regions

from Bio import SeqIO, PDB
from plip.modules.preparation import PDBComplex

def get_all_pi_helix_regions(pdb_data):

    ''' find all pi-helical regions in sequence '''
    ''' input: two-element list; first element is sequence (string) and second element is corresponding '''
    ''' secondary structure assignment (string); input example: ['MASLLFLPHTTPAPT', 'CCCCCCCCHHHHHHH'] '''
    ''' output: list of pi-helical sequence regions, each element of list is sequence of one region (string) '''
    ''' output example: ['YEAAVRKG', 'FIVEIIL'] '''

    pi_helix_indices  = []
    pi_helix_regions  = []
    regions_sequences = []

    # find all indices where pi-helix is present
    for ss in enumerate(pdb_data[1]):
        if ss[1] == 'I':
            pi_helix_indices.append(ss[0])

    # find coninuous regions, create list of list of continuous pi-helix indices
    idx_pos = pi_helix_indices[0]
    region  = [idx_pos]
    for idx in pi_helix_indices[1:]:
        if idx_pos + 1 == idx:
            idx_pos += 1
            region.append(idx)
        else:
            pi_helix_regions.append(region)
            region  = [idx]
            idx_pos = idx
    pi_helix_regions.append(region)

    # convert pi-helix indices to pi-helix sequences
    for region in pi_helix_regions:
        regions_sequences.append(''.join([ pdb_data[0][idx] for idx in region ]))

    return regions_sequences

def get_sequence_biopython(pdbpath, pdb_id):

    ''' get sequence from pdb file using Biopython '''
    ''' input: path to pdb file, pdb id with chain, example input: '/home/pdb/1ztm.pdb', '1ztm_A' '''
    ''' output: sequence (string) of desired chain '''

    sequence = ''
    handle   = open(pdbpath, "rU")

    for record in SeqIO.parse(handle, "pdb-atom"):
        # go through sequence records in pdb and select only that with desired chain id
        if record.id.split(':')[-1] == pdb_id.split('_')[-1]:
            sequence = record.seq

    handle.close()

    return sequence

'''
def get_first_residue_id(pdbname, pdbpath, pdb_id):

    ###FIXME ### - czy to jeszcze uzywane?

    '''  '''

    p          = PDB.PDBParser()
    structure  = p.get_structure(pdbname, pdbpath)
    model      = structure[0]
    chain_data = []

    for chain in model:
        if chain.id == pdb_id.split('_')[-1]:
            chain_data = chain

    for res in chain_data:
        print res
'''

def get_first_residue_id_dssp(pdbname, pdbpath, pdb_id):

    ''' get id of first residue in pdb file (pdb numbering) '''
    ''' use dssp data where pdb ids are keys of dssp data dictionary '''
    ''' input: pdbname, path to pdb file, pdb id with chain, example input: 'pdb1ztm', '/home/pdb/1ztm.pdb', '1ztm_A' '''
    ''' output: integer id of first residue in pdb, example output: 44 '''

    p          = PDB.PDBParser()
    structure  = p.get_structure(pdbname, pdbpath)
    model      = structure[0]
    dssp       = PDB.DSSP(model, pdbpath)
    chain_id   = pdb_id.split('_')[-1] # fetch chain id from pdb_id
    chain_data = []

    for res in list(dssp.keys()):
        # consider only residues of desired chain
        if res[0] == chain_id:
            chain_data.append(res)

    return chain_data[0][1][1]

def plip_analysis(pdbpath, pdb_id):

    ''' run PLIP analysis for pdb file '''

    mol = PDBComplex()
    mol.load_pdb(pdbpath)
    chain_id = pdb_id.split('_')[-1]
    ligand_interactions_dict = {}

    # invastigate only proteins with following ligands
    ligand_whitelist = set(['NAD','FAD','NDP','NAP','NAI','ADP','AMP','ATP','FMN','GTP','ADN','SAH','GDP',
                        '5GP','ADX','LSS','1ZZ','COA','ACO','01A','HMG','SAI','GRA','ST9','MCA','SAM',
                        'HEM','CLA','BCB','BCL','HEA','HEC','HEB','HAS'])



    # plip analysis
    ligand_list_from_mol = [ x.hetid for x in mol.ligands ]
    print('=> ' + str(len(ligand_list_from_mol)) + ' ligands to process')

    # if any ligand of interest present
#     if not ligand_whitelist.intersection(set(ligand_list_from_mol)):
#         print('=> No ligands of interest found')
#         #print(ligand_list_from_mol)
#         return ligand_interactions_dict

    mol.analyze()
    full_data_mol[pdb_id] = mol

    for bsid in [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]:
        # for every ligand encountered in pdb
        print bsid
        #print '|',

        interactions         = mol.interaction_sets[bsid] # Contains all interaction data
        print interactions.interacting_res
        interacting_residues = [ int(res[:-1]) for res in interactions.interacting_res if res[-1] == chain_id ]

        if interacting_residues:
            ligand_interactions_dict[bsid] = interacting_residues
    print '=> ' + str(len(ligand_interactions_dict)) + ' ligands interacting with chain ' + chain_id + ' found'
    print

    return ligand_interactions_dict

def plip_analysis_targets(pdbpath, pdb_id, pi_helix_pdb_ids):

    mol = PDBComplex()
    mol.load_pdb(pdbpath)
    chain_id = pdb_id.split('_')[-1]
    ligand_interactions = []

    # plip analysis
    ligand_list_from_mol = [ x.hetid for x in mol.ligands ]
    print('=> ' + str(len(ligand_list_from_mol)) + ' ligands to process')

    mol.analyze()
    full_data_mol[pdb_id] = mol

    print(target_pdb_plip_data[pdb_id])

    for bsid in [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]:
        # for every ligand encountered in pdb

        print bsid

        if bsid.split(':')[0] not in target_pdb_plip_data[pdb_id].keys():
            print 'BSID not in PLIP data', bsid
            continue

        if (not target_pdb_plip_data[pdb_id][bsid.split(':')[0]] or bsid.split(':')[1] != chain_id):
            continue
        print 'passed filters:', bsid
        #print '|',

        interactions         = mol.interaction_sets[bsid] # Contains all interaction data

        pi_helix_pdb_ids_all = [ res for reg in pi_helix_pdb_ids for res in reg ]
        #print pi_helix_pdb_ids
        #print pi_helix_pdb_ids_all
        for inter in interactions.all_itypes:
            if inter.resnr in pi_helix_pdb_ids_all:
                inter_type = str(type(inter)).strip('<>\'').split('.')[-1]
                print inter_type, inter.resnr, inter.restype, inter.restype_l
                ligand_interactions.append(inter)


    return ligand_interactions

def find_pi_helix_in_pdb(pi_helix_reg_sequences, pdb_sequence, first_res_pdb_id, skipped_count):

    #TODO: describe skipped_count and/or rewrite this functionality

    ''' find pdb residue ids corresponding to pi-helical regions detected in sequences '''
    ''' this has to be done when sequence do not corresponds perfectly to pdb '''
    ''' input: list of pi-helix sequences, sequence extracted from pdb, id of first residue in pdb '''
    ''' example input: [IHAVQAL, DLLAAMA], 'SLSTEAIHAVQALKRLTAADRSPPAATAAASAALGRLLRADLLAAMAELQR', 4 '''
    ''' output: list of pi-helix regions, each region is list of residue ids '''
    ''' example outpur: [[10, 11, 12, 13, 14, 15, 16], [44, 45, 46, 47, 48, 49, 50]]'''

    pi_helix_pdb_resi = []

    for region in pi_helix_reg_sequences:
        # if there is more than one occurance of region sequence in pdb_sequence then region is skipped
        if pdb_sequence.count(region) > 1:
            print 'Ambiguous region, skipping'
            skipped_count += 1
            continue
        reg_start_idx = pdb_sequence.find(region)
        pi_helix_pdb_resi.append([ resi + first_res_pdb_id for resi in range(reg_start_idx,reg_start_idx+len(region))])

    return pi_helix_pdb_resi

def check_if_interactions_with_pi_helix(lig_inter_dict, pi_helix_pdb_ids, data_pdb):

    regions_number = len(pi_helix_pdb_ids)

    for ligand in lig_inter_dict.items():
        #print ligand
        ligand_id = ligand[0].split(':')[0]
        #print ligand_id
        for pi_helix_region in pi_helix_pdb_ids:
            if set(ligand[1]) & set(pi_helix_region):
                # found interacting pi-helix region
                # update interactions dictionary
                print 'interaction with ', ligand_id, ' found'
                if ligand_id in data_pdb:
                    data_pdb[ligand_id] = True
                else:
                    data_pdb[ligand_id] = True
            else:
                if ligand_id not in data_pdb:
                    data_pdb[ligand_id] = False
