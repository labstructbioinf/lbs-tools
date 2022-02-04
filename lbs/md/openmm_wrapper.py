import os
from typing import NamedTuple, Tuple
from tempfile import TemporaryFile

import pandas as pd
import numpy as np
from openmm.app import ForceField, Modeller, Topology
from openmm.app.simulation import Simulation
from openmm.openmm import LangevinMiddleIntegrator
from openmm.app import PME, HBonds, StateDataReporter, PDBReporter
from openmm.app.pdbfile import PDBFile
from openmm.unit import nanometer, picosecond, kelvin, Quantity
from openmm.vec3 import Vec3


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
    numSteps: int = 10000
        
        
class OpenMM:
    '''
    OpenMM FAQ
    https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan
    '''
    def __init__(self, params: Params):
        self.params = params

    def prepare_pdb(self, path_pdb: str):
        '''
        prepare pdb file: add hydrogens, add solvent
        return:
            pdb
        '''
        assert os.path.isfile(path_pdb)
        # load pdb file
        pdb_tmp = PDBFile(path_pdb)
        # feed it into Modeller and add missing atoms
        # quasi centering
        positions_arr = pdb_tmp.getPositions(True)
        positions_arr = positions_arr - np.mean(positions_arr, axis=0)
        shift = np.array(self.params.boxSize)
        shift = self._as_array(shift)
        positions_arr += shift
        modeller = Modeller(pdb_tmp.topology, positions_arr)
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
        
    def run(self, pdb: Topology, path_output: str) -> pd.DataFrame:
        '''
        run simulation
        params:
            pdb (Topology)
            path_output (str)
        return:
            df (pd.DataFrame)
        '''
        _file = 'tmp.csv'
        system = self.forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer,
                                 constraints=HBonds)
        simulation = Simulation(pdb.topology, system, self.integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        with open(_file, 'wt') as tmp_file:
            simulation.reporters.append(PDBReporter(path_output, 1000))
            simulation.reporters.append(StateDataReporter(tmp_file, 1000, step=True,
                    potentialEnergy=True, density=True, temperature=True))
            simulation.step(self.params.numSteps)
        df = pd.read_csv(_file)
        os.remove(_file)
        return df
    
    @staticmethod
    def _as_array(iterable, unit = nanometer):
        '''
        convert iterable to quantity array
        '''
        arr = np.array(iterable)
        arr = Quantity(value=arr, unit=unit)
        return arr
        
        
    
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