import json
import os
from statistics import mean
import re

import numpy as np
import pandas as pd
from typing import Union, List

class AF2Analysis:


    def __init__(self, json_output: Union[str, List]) -> None:
        """
        Initialize the AF2_analysis object.

        Parameters:
            json_output (Union[str, List]): Either the path to the AF2 JSON output file or a list of JSON files.

        Returns:
            None
        """
        self.json_output = json_output
        if isinstance(self.json_output, str):
            self.json_output = [json_output]
        elif isinstance(self.json_output, list):
            self.json_output = json_output

    
    def read_json(self) -> pd.DataFrame:
        """
        Read the AF2 JSON output file.

        Parameters:
            None

        Returns:
            pd.DataFrame: The DataFrame containing the AF2 JSON output data.
        """
        df = pd.DataFrame()
        for _json in self.json_output:

            with open(_json) as f:
                _file = f.read()  
                parsed = json.loads(_file)
            
            plddt = mean(parsed['plddt'])
            ptm = float(parsed['ptm'])
            pae = np.asarray(parsed['pae'], dtype='float')
            pae = (pae + pae.T) / 2
            pdb_file = _json.replace('.json', '.pdb').replace('scores','unrelaxed')
            rank = re.findall('_rank_(\d{3})_',_json)
            pae = np.asarray(parsed['pae'], dtype='float')
            pae = (pae + pae.T) / 2

            j_df = pd.DataFrame({
                'json': [_json],  
                'plddt': [plddt],  
                'ptm': [ptm],
                'pae': [pae],  
                'pdb_file': [pdb_file], 
                'rank': [rank][0] 
            })
            df = pd.concat([df,j_df], axis=0, ignore_index=True)

        return df   
