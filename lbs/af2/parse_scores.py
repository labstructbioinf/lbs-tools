from array import array
import json
import os
import re
import numpy as np
from statistics import mean
import pandas as pd
from Bio import SeqIO
import warnings

warnings.filterwarnings("ignore")


class AF2ScoreParser:

    """
    Class parsing json files generated by AF2
    Required files: *scores.json
    Optional: pdbfiles
    af2_dir - directory with json and (preferably) pdb files after AF2 run
    domain length - length of each assymetric unit. Can be explicitly put in class instance
    oligo_state - number of assymetric units in structure. Can be explicitly put in class instance

    sample usage:
    instance = AF2ScoreParser('./test/').obtain_data()

    Returns dataframe
    """

    def __init__(self, af2_dir: str, domain_length=None, oligo_state=None) -> None:
        self.af2_dir = [os.path.join(af2_dir, x) for x in os.listdir(af2_dir)]
        # print(self.af2_dir)
        self.domain_length = domain_length
        self.oligo_state = oligo_state

    """Extracts length of assymetric unit"""

    def get_domain_length(self):
        pdb_list = [x for x in self.af2_dir if x.endswith(".pdb")]
        # print(pdb_list)
        if len(pdb_list) != 0:
            with open(pdb_list[0], "r") as pdb_file:
                for record in SeqIO.parse(pdb_file, "pdb-atom"):
                    assymetric_unit = len(record.seq)
                return assymetric_unit
        else:
            print(
                "No PDB files found, nor domain length and oligo state called, aborting\nPlease copy PDB files OR set domain_length and oligo_state accordingly"
            )

    """Extracts number of assymetric units"""

    def get_domain_number(self):
        pdb_list = [x for x in self.af2_dir if x.endswith(".pdb")]
        if len(pdb_list) != 0:
            with open(pdb_list[0], "r") as pdb_file:
                counter = 0
                for record in SeqIO.parse(pdb_file, "pdb-atom"):
                    counter += 1
                return counter
        else:
            print(
                "No PDB files found, nor domain length and oligo state called, aborting\nPlease copy PDB files OR set domain_length and oligo_state accordingly"
            )

    """Basing on json content extract plddt, ptm and pae array"""

    def parse_file_data(self, filename):
        with open(filename) as file:
            file = file.read()

        parsed = json.loads(file)
        plddt = round(mean(parsed["plddt"]), 2)
        ptm = parsed["ptm"]
        pae = np.asarray(parsed["pae"], dtype="float")
        return plddt, ptm, pae

    """Basing on filename extracts number or assymetric units or tries to extract it from the
    its content, along with rank and model number"""

    def get_file_info(self, filename):

        clean_path = filename.split("/")[-1]
        # print(re.findall(r"(.*)_(?:relaxed|unrelaxed)", clean_path)[0])
        _id = re.findall(r"(.*)_(?:relaxed|unrelaxed)", clean_path)[0]
        if not self.oligo_state:
            if "mer" in _id:
                mer = int(re.findall(r"(\d+)(?:mer)", clean_path)[0])
            elif "mer" not in _id:
                mer = self.get_domain_number()
            else:
                print("Assuming olgimer state equal 1 (monomer)")
                mer = 1
        else:
            mer = self.oligo_state
        rank = int(re.findall(r"(?:rank)_(\d+)", clean_path)[0])
        model = int(re.findall(r"(?:model)_(\d+)", clean_path)[0])

        return _id, rank, model, mer

    """Function calculating pae numbers and yielding list with results"""

    def gather_data(
        self,
        pae_arr: array,
        mer: int,
        rank: int,
        model: int,
        ptm: float,
        _id: str,
        plddt: float,
    ):
        if mer != 0:
            if mer > 1:
                assert int(pae_arr.shape[0] / self.domain_length) == int(
                    mer
                ), "Oligomerisation state provided does not match pae array shape"
            else:
                print("Monomeric structure being analyzed...")
            total_intra_pae = []
            total_inter_pae = []

            for pos, k in enumerate(np.hsplit(pae_arr, mer)):

                col = np.asarray(np.vsplit(k, mer))

                intra_pae = col[pos]
                intra_pae_score = intra_pae.mean()

                inter_cols = [
                    np.roll(np.arange(mer), -1)[pos],
                    np.roll(np.arange(mer), 1)[pos],
                ]

                inter_pae = col[inter_cols]
                inter_pae_score = inter_pae.mean()

                total_intra_pae.append(intra_pae_score)
                total_inter_pae.append(inter_pae_score)

            total_intra_pae = mean(total_intra_pae)
            total_inter_pae = mean(total_inter_pae)

            return [
                _id,
                mer,
                total_intra_pae,
                total_inter_pae,
                ptm,
                int(rank),
                int(model),
                plddt,
            ]
        else:
            print("***WARNING***, No oligomerisation found or specified!")

    """Data collections and yieding final dataframe"""

    def obtain_data(self):
        colnames = [
            "id",
            "oligomerisation_state",
            "intra_pae",
            "inter_pae",
            "ptm",
            "model",
            "model_rank",
            "plddt",
        ]
        dataframe = pd.DataFrame()

        if not self.domain_length:
            self.domain_length = self.get_domain_length()

        for file in self.af2_dir:
            if file.endswith("scores.json"):
                _id, rank, model, mer = self.get_file_info(file)
                plddt, ptm, pae_array = self.parse_file_data(file)

                df_to_append = pd.DataFrame(
                    self.gather_data(
                        _id=_id,
                        mer=mer,
                        rank=rank,
                        model=model,
                        ptm=ptm,
                        pae_arr=pae_array,
                        plddt=plddt,
                    )
                )

                dataframe = pd.concat(
                    [dataframe, df_to_append.T], axis=0, ignore_index=True
                )

        dataframe.columns = colnames
        # dataframe = dataframe.reset_index()

        return dataframe