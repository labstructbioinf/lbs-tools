import subprocess
import os
from typing import Union

import os
import subprocess
from typing import Union


class RunUSAlign:
    def __init__(self, usalign_binary: Union[os.PathLike, str]) -> None:
        """
        Initialize the RunUSAlign object.

        Parameters:
            usalign_binary (Union[os.PathLike, str]): The path to the USAlign binary executable.

        Returns:
            None
        """
        self.usalign_binary = usalign_binary

    def run_usalign(
        self,
        reference_pdb: Union[os.PathLike, str],
        target_pdb: Union[os.PathLike, str],
        output_filename: str,
        mm=1,
        ter=0,
        outfmt=2,
        pymol=False,
        matrix=False,
        pdb_rank=None,
    ) -> int:
        """
        Run the USAlign program with the given parameters.

        Parameters:
            reference_pdb (Union[os.PathLike, str]): The path to the reference PDB file.
            target_pdb (Union[os.PathLike, str]): The path to the target PDB file.
            output_filename (str): The name of the output file to save the USAlign results.
            mm (int): The mm option for USAlign (default: 1).
            ter (int): The ter option for USAlign (default: 0).
            outfmt (int): The output format for USAlign results (default: 2).
            pymol (bool): Whether to generate a PyMOL session file (default: False).
            matrix (bool): Whether to generate a rotation matrix file (default: False).
            pdb_rank (Union[None, int]): The rank of the PDB file (default: None).

        Returns:
            int: 0 if the analysis ran successfully, 1 otherwise.
        """
        print(f"Processing {reference_pdb} with {target_pdb}...")

        target_pdb_parent_directory  = os.path.basename(target_pdb)

        assert isinstance(reference_pdb, str) and isinstance(
            target_pdb, str
        ), "Invalid path or file does not exist"

        if pymol == False:
            command = (
                f"{self.usalign_binary} {target_pdb} {reference_pdb} -mm {mm} -ter {ter} -outfmt {outfmt} -o "
                f"{os.path.join(target_pdb_parent_directory, f'usalign_{pdb_rank}.txt')} "
                f"-m {os.path.join(target_pdb_parent_directory, f'matrix_{pdb_rank}.dat')} "
                f">> {output_filename}"
            )

        if matrix == False:
            command = (
                f"{self.usalign_binary} {target_pdb} {reference_pdb} -mm {mm} -ter {ter} -outfmt {outfmt} -o "
                f"{os.path.join(target_pdb_parent_directory, f'usalign_{pdb_rank}.txt')}  "
                f"{os.path.join(target_pdb_parent_directory, f'usalign_pymol_{pdb_rank}.pse')} "
                f">> {output_filename}"
            )

        if pymol == False and matrix == False:
            command = (
                f"{self.usalign_binary} {target_pdb} {reference_pdb} -mm {mm} -ter {ter} -outfmt {outfmt} -o "
                f"{os.path.join(target_pdb_parent_directory, f'usalign_{pdb_rank}.dat')} "
                f">> {output_filename}"
            )

        if pymol == True and matrix == False:
            command = (
                f"{self.usalign_binary} {target_pdb} {reference_pdb} -mm {mm} -ter {ter} -outfmt {outfmt} -o "
                f"{os.path.join(target_pdb_parent_directory, f'usalign_pymol_{pdb_rank}.pse')} "
                f"-m {os.path.join(target_pdb_parent_directory, f'matrix_{pdb_rank}.dat')} "
                f">> {output_filename}"
            )

        result = subprocess.run(command, shell=True)
        print(command)
        if result.returncode == 0:
            print("Analysis ran successfully")
            return 0
        else:
            with open("usalign_errors.txt", "a+") as fout:
                print("Analysis encountered an error")
                fout.write(f"{reference_pdb} FAILED WITH {target_pdb}\n")
            return 1
        
    def generate_slurm_job_input(self,
        reference_pdb: Union[os.PathLike, str],
        target_pdb_paths: Union[os.PathLike, str, list],
        input_filelist: str,
        input_batch_file: str,
        output_filename: str,
        mm=1,
        ter=0,
        outfmt=2,
        nodes_excluded=[0,1,2,3],
        memory=8,
        )  -> str:

        assert isinstance(reference_pdb, str), "Invalid path or file does not exist"

        command = ( f"for pdb in $(cat {input_filelist}.txt); do " +
                f"{self.usalign_binary} $pdb {reference_pdb} -mm {mm} -ter {ter} -outfmt {outfmt} -o >> {output_filename}"
            )

        if os.path.exists(input_filelist):
            with open(input_batch_file,'w+') as batch_in:
                batch_in.write(
                    f"#!/bin/bash \n\n#SBATCH -p cpu\n#SBATCH -n 1\n#SBATCH -N 1\n#SBATCH --exclude=edi0{nodes_excluded}\n#SBATCH --mem={memory}GB\n#SBATCH -J {output_filename}\n\n"
                )
                batch_in.write(command)
            return f"{input_batch_file} generated"
        
        else:
            with open(input_batch_file,'w+') as batch_in:
                batch_in.write(
                    f"#!/bin/bash \n\n#SBATCH -p cpu\n#SBATCH -n 1\n#SBATCH -N 1\n#SBATCH --exclude=edi0{nodes_excluded}\n#SBATCH --mem={memory}GB\n#SBATCH -J {output_filename}\n\n"
                )
                batch_in.write(command)
            with open(input_filelist+'txt', 'w') as fin:
                for pdb in target_pdb_paths:
                    fin.write(pdb + '\n')


        return f"{input_batch_file} and {input_filelist}.txt generated"

# RunUSAlign("/home/nfs/rmadaj/bins/usalign/USalign").generate_slurm_job_input('test.pdb', ['test1.pdb', 'test2.pdb'], 'test', 'test.sh', 'test.txt')