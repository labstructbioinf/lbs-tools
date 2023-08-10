import numpy as np
import sys
import glob

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

sys.path.append("tools/amber")

import tools.amber.process_amber as ac



class ProcessMultipleAmberOutputs:
    """
    Class to process multiple Amber simulation outputs and visualize the results.

    Parameters:
        input_csvs (list): List of paths to the CSV files containing Amber simulation data.
        sim_length (int): The length of the simulation in nanoseconds (ns).
        number_of_reps (int): The number of simulation replicates.
        trajectory_dump_frequency (int, optional): The frequency at which trajectory frames were saved (default is 5000).
        timestep (float, optional): The timestep used in the simulation in femtoseconds (fs) (default is 0.002).
        trajectory_names (list, optional): List of names for each trajectory (default is None).

    Raises:
        AssertionError: If any of the input CSV files do not exist, have fewer than 2 columns,
                       or the number of rows is not divisible by the number of replicates,
                       or if the number of rows in any CSV is not consistent with the simulation parameters.

    Returns:
        None
    """

    def __init__(
        self,
        input_csvs: list,
        sim_length: int,
        number_of_reps: int,
        trajectory_dump_frequency: int = 5000,
        timestep: float = 0.2,
        trajectory_names: list = None,
    ) -> None:
        self.input_csvs = input_csvs
        self.sim_length = sim_length
        self.number_of_reps = number_of_reps
        self.trajectory_dump_frequency = trajectory_dump_frequency
        self.timestep = timestep
        self.df = pd.DataFrame()

        self.df = pd.DataFrame()
        if trajectory_names is None:
            self.trajectory_names = [csv.split("/")[-1] for csv in self.input_csvs]
            for csv, name in zip(self.input_csvs, self.trajectory_names):
                single_df = ac.ProcessAmber(
                    csv,
                    self.sim_length,
                    self.number_of_reps,
                    self.trajectory_dump_frequency,
                    self.timestep,
                ).read_dataframe()
                single_df["structure"] = name
                self.df = pd.concat([self.df, single_df], axis=0, ignore_index=True)
        else:
            self.trajectory_names = trajectory_names
            assert len(self.trajectory_names) == len(
                self.input_csvs
            ), f"Number of trajectory names does not match number of input csvs ({len(self.trajectory_names)} != {len(self.input_csvs)})"
            for csv, name in zip(self.input_csvs, self.trajectory_names):
                single_df = ac.ProcessAmber(
                    csv,
                    self.sim_length,
                    self.number_of_reps,
                    self.trajectory_dump_frequency,
                    self.timestep,
                ).read_dataframe()
                single_df["structure"] = name
                self.df = pd.concat([self.df, single_df], axis=0, ignore_index=True)
        self.x_axis_points_number = int(
            (self.sim_length / timestep)
            * 1000
            / trajectory_dump_frequency
            * self.number_of_reps
            + 10
        )

    def read_dataframe(self) -> pd.DataFrame:
        """
        Get the DataFrame containing Amber simulation data.

        Returns:
            pd.DataFrame: The DataFrame containing Amber simulation data.
        """
        return self.df

    def melt_dataframe(self) -> pd.DataFrame:
        """
        Melt the DataFrame so that it uses "#Frame" and "structure" as id vars and the rest as value vars.

        Returns:
            pd.DataFrame: The melted DataFrame.
        """
        return pd.melt(
            self.df,
            id_vars=["#Frame", "structure"],
            var_name="name",
            value_name="value",
        )

    def analyze_feature(self, feature: str) -> pd.DataFrame:
        """
        Analyze the specified feature from the DataFrame.

        Parameters:
            feature (str): The name of the feature to analyze.

        Returns:
            pd.DataFrame: A DataFrame containing the mean and standard deviation of the specified feature for each trajectory.
        """
        df = (
            self.melt_dataframe().loc[self.melt_dataframe()["name"] == feature]
            .groupby("structure")['value']
            .agg({"mean", "std"})
            .round(2)
        )
        return df

    def vertical_lines(self, color):
        """
        Add vertical lines to the plot.

        Parameters:
            color: The color of the vertical lines.

        Returns:
            None
        """
        line_position = [
            x
            for x in range(
                0,
                self.x_axis_points_number,
                int(self.x_axis_points_number / self.number_of_reps),
            )
        ]
        for x in line_position:
            plt.axvline(x, color="w", linewidth=3)
            plt.axvline(x, color="black", linestyle=":")

    def plot_feature(
        self,
        feature: str,
        custom_ticks: list = None,
        vertical_lines: bool = False,
        png: bool = False,
        col_wrap: int = 2,
        height: int = 4,
        aspect: float = 1.5,
    ):
        """
        Plot the specified feature for each trajectory.

        Parameters:
            feature (str): The name of the feature to plot.
            custom_ticks (list, optional): Custom labels for x-axis ticks (default is None).
            vertical_lines (bool, optional): Whether to add vertical lines to the plot (default is False).

        Returns:
            str: A message indicating the figure has been saved.
        """
        if custom_ticks is None:
            custom_ticks = [str(x) for x in range(1, self.number_of_reps + 1)]

        feature_df = self.melt_dataframe().loc[self.melt_dataframe()["name"] == feature]
        g = sns.FacetGrid(
            feature_df,
            col="structure",
            col_wrap=col_wrap,
            height=height,
            aspect=aspect,
        )
        g.map(sns.lineplot, "#Frame", "value")
        sns.set(font_scale=1.5, style="whitegrid")

        if vertical_lines:
            g.map(self.vertical_lines)

        g.set(
            xlabel="Simulation",
            ylabel=feature,
            xticks=np.arange(
                (self.x_axis_points_number / self.number_of_reps) / 2,
                self.x_axis_points_number
                + (self.x_axis_points_number / self.number_of_reps / 2),
                self.x_axis_points_number / self.number_of_reps,
            ),
            xticklabels=custom_ticks,
        )
        if png:
            g.savefig(f"{feature}.png")
            return f"Figure {feature}.png saved to the current directory"
        else:
            return "Done!"


t=ProcessMultipleAmberOutputs(    input_csvs = glob.glob('./files/amber_csvs/*.csv'),
    sim_length=30,
    number_of_reps=10,
    trajectory_dump_frequency=5000,
    timestep=0.002)