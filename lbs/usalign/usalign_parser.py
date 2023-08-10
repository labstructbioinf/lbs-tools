import os
from typing import Union

import pandas as pd

import os
import pandas as pd
from typing import Union


class USalign_parser:
    def __init__(self, user_input: Union[str, pd.DataFrame]) -> None:
        """
        Initialize the USalign_parser object.

        Parameters:
            user_input (Union[str, pd.DataFrame]): Either the path to the USAlign output file (CSV format)
                                                  or a pandas DataFrame containing USAlign results.

        Returns:
            None
        """
        if isinstance(user_input, str):
            assert os.path.exists(user_input), "Invalid path or file does not exist"
            # load csv, drop all lines starting with #



            self.df = pd.read_csv(
                user_input,
                sep="\s+",
                comment="#",
                skip_blank_lines=True,
                names=[
                    "target",
                    "template",
                    "tm1",
                    "tm2",
                    "rmsd",
                    "id1",
                    "id2",
                    "idali",
                    "docked_seqlength",
                    "template_seqlength",
                    "aligned_length",
                ],
            ).reset_index(drop=True)
            assert self.df.shape[0] > 0, "Empty dataframe"
        elif isinstance(user_input, pd.DataFrame):
            self.df = user_input

    def read_usalign_output(self):
        """
        Process the USAlign output DataFrame and extract relevant information.

        Returns:
            pd.DataFrame: The processed DataFrame with additional columns.
        """
        self.df["target_path"] = self.df["target"].apply(lambda x: x.split(':')[0])
        self.df["target"] = (
            self.df["target"].str.split("/").str[-1].str.split(".").str[0]
        )

        return self.df

    def get_top_n_by_selected_column(
        self, column: str, n: int = 10, ascending: bool = False
    ):
        """
        Get the top or bottom 'n' rows based on the values in the specified 'column'.

        Parameters:
            column (str): The column name based on which to select the top or bottom 'n' rows.
            n (int): The number of rows to select (default: 10).
            ascending (bool): If True, select the top 'n' rows; otherwise, select the bottom 'n' rows (default: False).

        Returns:
            pd.DataFrame: A DataFrame containing the top or bottom 'n' rows based on the specified 'column'.
        """
        return self.df.sort_values(by=column, ascending=ascending).head(n)
    
    def add_column(self, column: str, value: Union[str, int, float]):
        """
        Add a column to the DataFrame with the specified value.

        Parameters:
            column (str): The name of the column to add.
            value (str): The value to add to the column.

        Returns:
            None
        """
        self.df[column] = value