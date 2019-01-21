import pandas as pd


class Plugin(object):

    def __init__(self):
        pass

    @staticmethod
    def add_col(df, data_dict, added_col_name):
        data_dict_fixed = {key: [value] for key, value in data_dict.items()}
        data_df = pd.DataFrame.from_dict(data_dict_fixed, orient='index', columns=[added_col_name])
        df = pd.concat((df, data_df), axis=1, sort=True)
        return df
