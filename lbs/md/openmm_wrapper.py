import os
import json
from typing import NamedTuple, Tuple, Union
from tempfile import NamedTemporaryFile

import pandas as pd
import numpy as np
from openmm.app import ForceField, Modeller, Topology
from openmm.app.simulation import Simulation
from openmm.openmm import LangevinMiddleIntegrator
from openmm.app import PME, HBonds, StateDataReporter, PDBReporter
from openmm.app.pdbfile import PDBFile
from openmm.unit import nanometer, picosecond, kelvin, Quantity
from openmm.vec3 import Vec3
try:
    from pdbfixer import PDBFixer
except ImportError as e:
    print(e)
    
from .tools import remove_hetatom


class Params:
    '''
    boxSize: Tuple[float] = (5.0, 5.0, 5.0) # nanometer
    temperature: float = 300 #kelvins
    frictionCoeff: float = 1 # 1/picosecond
    stepSize: float = 0.002 # picosecond
    numSteps: int = 10000 
    '''
    boxSize: Tuple[float] = (6.0, 6.0, 6.0) # nanometer
    temperature: float = 300 #kelvins
    frictionCoeff: float = 1 # 1/picosecond
    stepSize: float = 0.002 # picosecond
    numSteps: int = 50000
    saveStep: int = 1000 # how often openmm write snapshot
    # list of params to backup
    simattr = [
        'boxSize',
        'temperature',
        'frictionCoeff',
        'stepSize',
        'numSteps',
        'saveStep'
        ]
    def __repr__(self):
        out = f'''
        boxSize: {self.boxSize}
        temperature: {self.temperature}
        frictionCoeff: {self.frictionCoeff}
        stepSize: {self.stepSize}
        numSteps: {self.numSteps}
        saveStep: {self.saveStep}
        '''
        return out

    @property
    def simulationTimeNano(self):
        '''
        simulation time in pico seconds
        '''
        return round(self.numSteps*self.stepSize/1000, 4)

    def to_json(self) -> dict:
        '''
        return simulation params as json
        '''
        paramsjson = dict()
        for attr_name in self.simattr:
            attr = getattr(self, attr_name)
            # cast to valid json types
            if isinstance(attr, (np.float16, np.float32, np.float64, np.float128)):
                attr = float(attr)
            elif isinstance(attr, (np.int8, np.int16, np.int32, np.int64)):
                attr = int(attr)
            elif isinstance(attr, tuple):
                attr = tuple([float(val) for val in attr])
            paramsjson[attr_name] = attr
        paramsjson["simulationTimeNano"] = int(self.simulationTimeNano)
        return paramsjson
        
        
class OpenMM:
    '''
    main class for MD simulations via openMM package
    OpenMM FAQ:
    https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan
    '''
    def __init__(self, params: Params, device: str):
        self.params = params
        self.device = device

    def prepare_pdb(self, path_pdb: Union[str, PDBFile], rm_residue=[]):
        '''
        prepare pdb file: add hydrogens, add solvent
        return:
            pdb
        '''
        from_fixer = False
        if isinstance(path_pdb, str):
            assert os.path.isfile(path_pdb)
            pdb_tmp = PDBFile(path_pdb)
        else:
            from_fixer = True
            pdb_tmp = path_pdb
        # feed it into Modeller and add missing atoms
        # quasi centering
        
        positions_arr = self._asarray(pdb_tmp.positions._value)
        positions_arr = positions_arr - np.mean(positions_arr, axis=0)
        shift = self._asarray(self.params.boxSize)/2
        positions_arr += shift
        modeller = Modeller(pdb_tmp.topology, positions_arr)
        if len(rm_residue) != 0:
            modeller.delete(toDelete=rm_residue)
        if not from_fixer:
            _ = modeller.addHydrogens(self.forcefield)
        modeller.addSolvent(self.forcefield, boxSize=Vec3(*self.params.boxSize)*nanometer)
        return modeller
    
    def prepare_components(self):
        '''
        initialize frocefield and integrator objects
        '''
        temperature = float(self.params.temperature)*kelvin
        frictionCoeff = float(self.params.frictionCoeff)/picosecond
        stepSize = float(self.params.stepSize)*picosecond
        self.forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        self.integrator = LangevinMiddleIntegrator(temperature,
                                                   frictionCoeff,
                                                   stepSize)
        
    def run(self, pdb: Topology, path_output: str, remove_hetatm: bool = True) -> pd.DataFrame:
        '''
        run simulation
        Params:
            pdb (Topology)
            path_output (str): file to save output simulation
            remove_hetatm (bool)
        Returns:
            df (pd.DataFrame)
        '''
        assert isinstance(path_output, str)
        assert path_output.endswith('.pdb'), "simulation `path_output` must end with .pdb extension"
        tf = NamedTemporaryFile(delete=False)
        _file = tf.name
        tf.close()
        # save used params
        paramsdict = self.params.to_json()
        paramsfile = path_output.replace('.pdb', '_params.json')
        with open(paramsfile, 'wt') as f:
            json.dump(paramsdict, f)
        # initialize sumulation environemnt
        system = self.forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer,
                                 constraints=HBonds)
        simulation = Simulation(pdb.topology, system, self.integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        with open(_file, 'wt') as tmp_file:
            simulation.reporters.append(PDBReporter(path_output, self.params.saveStep))
            simulation.reporters.append(StateDataReporter(tmp_file,
                                                          self.params.saveStep,
                                                          step=True,
                                                          potentialEnergy=True,
                                                          totalEnergy=True,
                                                          density=True,
                                                          temperature=True))
            simulation.step(self.params.numSteps)

        df = pd.read_csv(_file)
        os.remove(_file)
        # remove water and other stuff
        if remove_hetatm:
            remove_hetatom(path_output, path_output)
        return df
    
    @staticmethod
    def _asarray(iterable, unit = nanometer):
        '''
        convert iterable to quantity array
        '''
        arr = np.array(iterable)
        arr = Quantity(value=arr, unit=unit)
        return arr

    def _save_params(self):
        pass
        
        
class Patcher:
    '''
    Main class for fixing common pdb issues such as:
    * remove HETATOMS
    * add missing hydrogens
    * remove/handle custom residues
    '''
    residues = (
    'ALA', 'ARG', 'ASH', 'ASN',
    'ASP', 'CYM', 'CYS', 'CYX',
    'GLH', 'GLN', 'GLU', 'GLY',
    'HID', 'HIE', 'HIP', 'HYP',
    'HIS',
    'ILE', 'LEU', 'LYN', 'LYS',
    'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
    )
    def __init__(self, forcefield):
        '''
        params:
            forcefield
        '''
        self.ff = forcefield
        self.template_db = self.ff._templates
    
    def read_and_repair(self, path_pdb: str):
        '''
        params:
            path_pdb (str) path to structrue
        return:
            pdb (PDBFixer object)
            invalid_residues (list[residues]) residues to remove
        '''
        assert os.path.isfile(path_pdb)
        fixer = PDBFixer(filename=path_pdb)
        #fixer.removeHeterogens(keepWater=False)
        #fixer.addMissingHydrogens()
        #fixer.findNonstandardResidues()
        #fixer.replaceNonstandardResidues()

        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)


        invalid_residues = self._check_residues(fixer.topology)
        
        return fixer, invalid_residues
        
    def _check_residues(self, topology):
        '''
        iterates over structure for non-residue objects
        '''
        invalid_residues = list()
        for residue in topology.residues():
            if residue.name not in self.residues:
                invalid_residues.append(residue)
        return invalid_residues 
    
    def get_template_dict(self, name):
        '''
        return forcefield `name` templateData object
        '''
        if name in self.template_db:
            template = self.template_db[name].atoms
            map_name_type = {atom.name : atom.type for atom in template}
        else:
            map_name_type = None
        return map_name_type
    
    def set_template_atom_type(self, template, map_name_type):
        for atom in template.atoms:
            atom.type = map_name_type[atom.name]
        return template
    
    
    
def find_residue_atoms(topology: Topology, chain: str, residue: str):
    '''
    pick residue from openmm topology object
    '''
    chains_available = {chain.id: chain for chain in topology.chains()}
    if chain in chains_available:
        data = chains_available[chain]
    else:
        raise KeyError(f'chain: {chain} not in {chains_available.keys()}')
    residues_available = {res.id: res for res in data.residues()}
    if residue in residues_available:
        data = residues_available[residue]
    else:
        raise KeyError(f'residue: {residue} not in {residues_available.keys()}')
    return data
