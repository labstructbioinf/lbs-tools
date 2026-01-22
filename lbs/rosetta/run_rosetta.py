# if for some reason one cannot use  pyrosetta...

import os
import subprocess
from typing import List


class RunRosetta:
    def __init__(self, rosetta_bin_directory: str) -> None:
        self.bin_dir = rosetta_bin_directory
        assert os.path.isdir(
            self.bin_dir
        ), f"Rosetta bin directory {self.bin_dir} does not exist"

    def run_fast_relax(
        self, input_file: str, out_dir: str, extra_flags: str = "", silent_file: bool = False,
    ) -> None:
        """
        Run minimisation on a input file, output in the same directory where the input_file is
        """
        os.makedirs(out_dir, exist_ok=True)
        outpath = out_dir
        assert os.path.isfile(input_file), f"File {input_file} does not exist"
        
        if not silent_file:
            command = f"{self.bin_dir}/minimize.linuxgccrelease -s {input_file} -out:path:all {outpath} -out:file:scorefile 'score_minim.sc' -ignore_zero_occupancy false \
    -relax:constrain_relax_to_start_coords -relax:default_repeats 5 -relax:ramp_constraints false -in:file:fullatom -ex1 -ex2 -nstruct 1 {extra_flags}"
            
            subprocess.run(command, shell=True, check=True)
        else:
            command = f"{self.bin_dir}/minimize.linuxgccrelease -in:file:silent {input_file} -out:path:all {outpath} -out:file:scorefile 'score_minim.sc' -ignore_zero_occupancy false \
    -relax:constrain_relax_to_start_coords -relax:default_repeats 5 -relax:ramp_constraints false -in:file:fullatom -ex1 -ex2 -nstruct 1 {extra_flags}"
            
            subprocess.run(command, shell=True, check=True)


    def run_interface_analysis(
        self,
        input_file: str,
        out_file: str,
        interchain_interface: list,
        extra_flags: str = "",
        silent_file: bool = False,
    ) -> None:
        """
        Run interface analysis on a pdb file
        """
        assert os.path.isfile(input_file), f"File {input_file} does not exist"


        if not silent_file:
            command = f"{self.bin_dir}/InterfaceAnalyzer.linuxgccrelease -s {input_file} -interface {'_'.join(interchain_interface)} -pack_separated {extra_flags} -out:file:score_only -out:file:scorefile {out_file}"
            subprocess.run(command, shell=True, check=True)
        else:
            command = f"{self.bin_dir}/InterfaceAnalyzer.linuxgccrelease -in:file:silent {input_file} -interface {'_'.join(interchain_interface)} -pack_separated {extra_flags} -out:file:score_only -out:file:scorefile {out_file}"
            subprocess.run(command, shell=True, check=True)

    def run_residue_energy_breakdown(
        self, input_file: str, out_file: str, extra_flags: str, silent_file: bool = False,
    ) -> None:
        """
        Run decomposition on a input_file
        """
        assert os.path.isfile(input_file), f"File {input_file} does not exist"
        if not silent_file:
            command = f"{self.bin_dir}/residue_energy_breakdown.linuxgccrelease -s {input_file} -out:file:silent {out_file} {extra_flags}"
            subprocess.run(command, shell=True, check=True)
        else:
            command = f"{self.bin_dir}/residue_energy_breakdown.linuxgccrelease -in:file:silent {input_file} -out:file:silent {out_file} {extra_flags}"
            subprocess.run(command, shell=True, check=True)


    def extract_pdbs_by_ids(self, input_silent_file, tags: List[str], out_dir: str) -> None:
        """
        Extract pdbs from a silent file by tags
        """
        assert os.path.isfile(input_silent_file), f"File {input_silent_file} does not exist"
        os.makedirs(out_dir, exist_ok=True)
        command = f"{self.bin_dir}/extract_pdbs.linuxgccrelease -in:file:silent {input_silent_file} -tags {' '.join(tags)} -out:prefix {out_dir}/"
        subprocess.run(command, shell=True, check=True)