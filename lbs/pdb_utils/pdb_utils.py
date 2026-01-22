import os
from typing import Union

from biopandas.pdb import PandasPdb
import numpy as np


class PDBUtils:
    def __init__(self, pdb_file: Union[os.PathLike, str]) -> None:
        self.pdb_file = pdb_file
        assert os.path.exists(pdb_file), f"{pdb_file} does not exist."

    def calculate_distance_between_residues(
        self,
        residue1: int,
        residue2: int,
        ligand=False,
        residue1_chain: str = "A",
        residue2_chain: str = "B",
    ) -> float:
        """
        Calculate the distance between two residues.

        Parameters:
            residue1 (str): The first residue.
            residue2 (str): The second residue or ligand

        Returns:
            float: The distance between the two residues.
        """
        pdb = PandasPdb().read_pdb(self.pdb_file)
        residue1 = pdb.df["ATOM"][(pdb.df["ATOM"]["residue_number"] == residue1) & (pdb.df["ATOM"]["chain_id"] == residue1_chain)]
        if ligand:
            residue2 = pdb.df["HETATM"][pdb.df["HETATM"]["residue_name"] == residue2]
        else:
            residue2 = pdb.df["ATOM"][(pdb.df["ATOM"]["residue_number"] == residue2) & (pdb.df["ATOM"]["chain_id"] == residue2_chain)]
        distance = np.linalg.norm(residue1[['x_coord', 'y_coord', 'z_coord']].values - residue2[['x_coord', 'y_coord', 'z_coord']].values)
        return round(distance,1)

    def calculate_distance_between_atoms(
        self, atom1: int, atom2: int, ligand=False
    ) -> float:
        """
        Calculate the distance between two atoms.

        Parameters:
            atom1 (str): The first atom.
            atom2 (str): The second atom of protein or ligand.

        Returns:
            float: The distance between the two atoms.
        """
        pdb = PandasPdb().read_pdb(self.pdb_file)
        atom1 = pdb.df["ATOM"][pdb.df["ATOM"]["atom_number"] == atom1]
        if ligand:
            atom2 = pdb.df["HETATM"][pdb.df["HETATM"]["atom_number"] == atom2]
        else:
            atom2 = pdb.df["ATOM"][pdb.df["ATOM"]["atom_number"] == atom2]
        distance = np.linalg.norm(atom1[['x_coord', 'y_coord', 'z_coord']].values - atom2[['x_coord', 'y_coord', 'z_coord']].values)
        return round(distance,1)

    def truncate_structure_to_backbone(self):
        """
        Truncate the structure to only include the backbone atoms.

        Returns:
            PandasPdb: The truncated structure.
        """
        pdb = PandasPdb().read_pdb(self.pdb_file)
        pdb.df["ATOM"] = pdb.df["ATOM"][
            pdb.df["ATOM"]["atom_name"].isin(["N", "CA", "C", "O"])
        ]
        return pdb

    def truncate_structure_by_residue_indices(
        self, end: int, start: int = 1, chain: str = "A"
    ):
        """
        Truncate the structure to only include the selected residues.

        Parameters:
            start (int): The start residue index.
            end (int): The end residue index.
            chain (str): The chain ID.

        Returns:
            PandasPdb: The truncated structure.
        """
        pdb = PandasPdb().read_pdb(self.pdb_file)
        pdb.df["ATOM"] = pdb.df["ATOM"][
            (pdb.df["ATOM"]["residue_number"] >= start)
            & (pdb.df["ATOM"]["residue_number"] <= end)
            & (pdb.df["ATOM"]["chain_id"] == chain)
        ]
        return pdb
