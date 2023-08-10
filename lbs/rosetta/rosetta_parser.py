import os
import typing

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


import os
import pandas as pd
import typing
import seaborn as sns
import matplotlib.pyplot as plt


class RosettaParser:
    """
    Class to parse Rosetta output files.

    Parameters:
        user_input (Union[os.PathLike, pd.DataFrame]): Either the path to the Rosetta output file (CSV format)
                                                      or a pandas DataFrame containing Rosetta results.

    Returns:
        None
    """

    def __init__(self, user_input: typing.Union[os.PathLike, pd.DataFrame]):
        if isinstance(user_input, pd.DataFrame):
            self.df = user_input

        elif isinstance(user_input, str):
            assert os.path.exists(user_input), "Invalid path or file does not exist"

            with open(user_input, "r") as f:
                if f.readline().split()[0] == "SEQUENCE:":
                    self.df = pd.read_csv(user_input, sep="\s+", skiprows=1)
                else:
                    self.df = pd.read_csv(user_input, sep="\s+")
            assert self.df.shape[0] > 0, "Empty dataframe"

    def read_dataframe(self):
        """
        Get the parsed DataFrame.

        Returns:
            pd.DataFrame: The parsed DataFrame containing the Rosetta results.
        """
        return self.df

    def get_top_n_by_selected_column(self, column: str, n: int = 10):
        """
        Get the top 'n' rows based on the values in the specified 'column'.

        Parameters:
            column (str): The column name based on which to select the top 'n' rows.
            n (int): The number of rows to select (default: 10).

        Returns:
            pd.DataFrame: A DataFrame containing the top 'n' rows based on the specified 'column'.
        """
        return self.df.sort_values(by=column, ascending=False).head(n)

    def add_column_to_dataframe(
        self, column: str, values: typing.Union[typing.List, str, float, int]
    ):
        """
        Add a new column to the DataFrame.

        Parameters:
            column (str): The name of the new column.
            values (Union[List, str, float, int]): The values to be added to the new column.

        Returns:
            None
        """
        self.df[column] = values

    def change_strings_in_selected_column_by_user_input(
        self, column: str, original: typing.Union[str, int, float], replacement: str
    ):
        """
        Replace specific strings in the selected column with the given replacement.

        Parameters:
            column (str): The column name where replacement will be performed.
            original (Union[str, int, float]): The string or value to be replaced.
            replacement (str): The string to replace the original values.

        Returns:
            None
        """
        if original == "all":
            self.df[column] = replacement
        else:
            self.df[column] = self.df[column].str.replace(str(original), replacement)

    def plot_top_n_by_selected_column(self, column: str, n: int = 10):
        """
        Plot a bar chart of the top 'n' rows based on the values in the specified 'column'.

        Parameters:
            column (str): The column name based on which to select the top 'n' rows.
            n (int): The number of rows to select (default: 10).

        Returns:
            None
        """
        top_n = self.get_top_n_by_selected_column(column, n)
        sns.barplot(x=top_n["description"], y=top_n[column])
        plt.xticks(rotation=90)
        plt.show()

    def plot_histogram(self, column: str):
        """
        Plot a histogram for the values in the specified column.

        Parameters:
            column (str): The column name for which to plot the histogram.

        Returns:
            None
        """
        sns.histplot(self.df[column])
        plt.show()
