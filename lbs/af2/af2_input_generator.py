import os
from typing import List, Tuple, Union
import math
import re

class AF2InputGenerator:
    def __init__(
        self, input_fasta: Union[str, List], output_dir: Union[os.PathLike, str]
    ) -> None:
        """
        Initialize the AF2_input_generator object.

        Parameters:
            input_fasta (Union[str, List]): Either the path to the input FASTA,  A3M file or a list of FASTA sequences.
            output_dir (Union[os.PathLike, str]): The directory where the AF2 input files will be written.

        Returns:
            None
        """


        self.input_fasta = input_fasta
        self.output_dir = output_dir
        if isinstance(self.input_fasta, str):
            self.input_fasta = input_fasta
        elif isinstance(self.input_fasta, list):
            self.input_fasta = ":".join(self.input_fasta)

    def check_sequence(self) -> None:
        """
        Check if the input sequence is a valid protein sequence.

        Parameters:
            None

        Returns:
            None

        Raises:
            ValueError: If the input sequence contains invalid characters or is empty.
        """
        sequence = self.input_fasta.replace(" ", "")
        valid_chars = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY:]+$")
        if not valid_chars.match(sequence):
            raise ValueError("Invalid protein sequence. It contains invalid characters.")
        if not sequence:
            raise ValueError("Invalid protein sequence. It is empty.")

    def write_af2_inputs(self, oligo_state=1) -> str:
        """
        Write AF2 input files.

        Parameters:
            oligo_state (int): The number of oligo states.

        Returns:
            None
        """
        os.makedirs(self.output_dir, exist_ok=True)

        if self.input_fasta.endswith(".a3m"):
            # copy  file to output dir
            os.system(f"cp {self.input_fasta} {self.output_dir}")

            return "a3m detected, copying file to output dir"

        else:
            self.check_sequence()

            output_name = (
                self.output_dir.split("/")[-1]
                if "/" in self.output_dir
                else self.output_dir
            )
            if oligo_state > 1:
                self.input_fasta = ":".join([self.input_fasta] * oligo_state)

            with open(os.path.join(self.output_dir, "af2input_fasta.csv"), "w") as fout:
                fout.write("id,sequence\n")
                fout.write(f"{output_name},{self.input_fasta}")

                return "fasta detected, writing af2_input.fasta"

    def write_af2_exec(
        self,
        n_cores=7,
        amber=True,
        num_models=5,
        num_recycles=5,
        nodes_excluded=[0,1,2,3],
        memory=16,
        use_dropout=False,
        num_seeds=1,
        max_seqs=False,
        max_extra_seqs=False
    ) -> None:
        """
        Write the AF2 execution script.

        Parameters:
            n_cores (int): The number of CPU cores to use.
            amber (bool): Whether to use Amber relaxation (True) or not (False).
            num_models (int): The number of models to generate.
            num_recycles (int): The number of recycle steps.
            nodes_excluded (List[int]): List of node numbers to exclude.
            memory (int): The memory (in GB) required for the execution.
            use_dropout (bool): Whether to use dropout for models (True) or not (False).
            num_seeds (int): The number of random seeds to use.
            max_seqs (int): The maximum number of sequences (should be a power of 2).
            max_extra_seqs (int): The maximum number of extra sequences.


        Returns:
            None
        """
        if amber:
            amber = " --amber"
        else:
            amber = ""
        if use_dropout:
            use_dropout = " --use-dropout"
        else:
            use_dropout = ""
        if max_seqs:
            if not (math.log2(max_seqs).is_integer()):
                raise ValueError("max_seq must be a power of 2.")
            
            max_seqs = f" --max-seqs {max_seqs}"
        else:
            max_seqs = ""
        if max_extra_seqs:
            if not (math.log2(max_seqs).is_integer()):
                raise ValueError("max_extra_seqs must be a power of 2.")
            max_extra_seqs = f" --max-extra-seq {max_extra_seqs}"
        else:
            max_extra_seqs = ""

        # remove whitespace from nodes_excluded




        nodes_excluded = str(nodes_excluded).replace(" ", "")



        with open(os.path.join(self.output_dir, "af2_run.sh"), "w") as af2exec:
            af2exec.write(
                f"#!/bin/bash \n\n#SBATCH -p gpu\n#SBATCH -n {n_cores}\n#SBATCH -N 1\n#SBATCH --gres=gpu:1\n#SBATCH --exclude=edi0{nodes_excluded}\n#SBATCH --mem={memory}GB\n#SBATCH -J {self.output_dir}\nsource /opt/miniconda3/bin/activate cf_1.5\n"
            )  # necessary imports for AF2 to run
            if self.input_fasta.endswith(".a3m"):
                af2exec.write(
                    f"colabfold_batch *.a3m . --num-models {num_models} --num-recycle {num_recycles} {amber} {use_dropout} --num-seeds {num_seeds} {max_seqs} {max_extra_seqs}"
                )
            else:
                af2exec.write(
                    f"colabfold_batch af2input_fasta.csv . --num-models {num_models} --num-recycle {num_recycles} {amber} {use_dropout} --num-seeds {num_seeds} {max_seqs} {max_extra_seqs}"
                )
