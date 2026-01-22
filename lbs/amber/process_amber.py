import os
import pandas as pd


class ProcessAmber:
    """
    Class to process Amber simulation data from a CSV file.

    Parameters:
        input_csv (str): The path to the CSV file containing Amber simulation data.
        sim_length (int): The length of the simulation in nanoseconds (ns).
        number_of_reps (int): The number of simulation replicates.
        trajectory_dump_frequency (int, optional): The frequency at which trajectory frames were saved (default is 5000).
        timestep (float, optional): The timestep used in the simulation in femtoseconds (fs) (default is 0.002).

    Raises:
        AssertionError: If the input CSV file does not exist, has fewer than 2 columns,
                       or the number of rows is not divisible by the number of replicates,
                       or if the number of rows in the CSV is not consistent with the simulation parameters.

    Returns:
        None
    """

    def __init__(
        self,
        input_csv: str,
        sim_length: int,
        number_of_reps: int,
        trajectory_dump_frequency: int = 5000,
        timestep: float = 0.002,
    ) -> None:
        self.sim_length = sim_length
        self.number_of_reps = number_of_reps
        self.df = pd.read_csv(input_csv, sep="\s+")
        # assert input csv exists
        assert os.path.exists(input_csv), "Input csv does not exist"
        # assert that number of columns in csv is higher than 2
        assert self.df.shape[1] >= 2, "Number of columns in csv is less than 2"
        # assert that number of rows in csv divided by number of reps is an integer
        assert (
            self.df.shape[0] % self.number_of_reps == 0
        ), "Number of rows in csv is not divisible by number of reps, check sim length"
        # assert that sim_length /  timestep * 1000 *  number_of_reps + number_of_reps is equal to number of rows in csv
        frame_numer_expected = (
            self.sim_length / timestep
        ) * 1000 / trajectory_dump_frequency * self.number_of_reps + 10
        assert (
            frame_numer_expected == self.df.shape[0]
        ), f"Number of rows in csv is not equal to sim_length / timestep * 1000 * number_of_reps + number_of_reps  in {input_csv}. {frame_numer_expected} != {self.df.shape[0]}"

    def read_dataframe(self) -> pd.DataFrame:
        """
        Get the DataFrame containing Amber simulation data.

        Returns:
            pd.DataFrame: The DataFrame containing Amber simulation data.
        """
        return self.df

    def obtain_statistics_from_df(self) -> pd.DataFrame:
        """
        Obtain statistics from the DataFrame.

        Returns:
            pd.DataFrame: A DataFrame containing the mean and standard deviation of each column in the original DataFrame.
        """
        # Drop #Frame column
        self.df = self.df.drop("#Frame", axis=1)
        # Obtain mean and standard deviation from the DataFrame by pivoting the columns
        self.df_mean = self.df.mean(axis=0)
        self.df_std = self.df.std(axis=0)
        # Combine both DataFrames into one
        self.df_mean = pd.concat([self.df_mean, self.df_std], axis=1)
        self.df_mean.columns = ["mean", "std"]

        return self.df_mean
