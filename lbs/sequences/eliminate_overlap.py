import pandas as pd

class eliminate_overlap:
    def __init__(self, df, min_overlap=25):
        self.df = df
        self.min_overlap = min_overlap
        self.filtered_df = self._process()

    def _overlap(self, range1, range2):
        start1, end1 = range1
        start2, end2 = range2
        overlap_length = min(end1, end2) - max(start1, start2) + 1
        return overlap_length >= self.min_overlap

    def _group_overlapping_ranges(self, df):
        groups = []
        current_group = []
        df_sorted = df.sort_values(by='start').reset_index(drop=True)

        for idx, row in df_sorted.iterrows():
            if not current_group:
                current_group.append(idx)
            else:
                last_row = df_sorted.loc[current_group[-1]]
                if self._overlap((row['start'], row['stop']), (last_row['start'], last_row['stop'])):
                    current_group.append(idx)
                else:
                    groups.append(current_group)
                    current_group = [idx]

        if current_group:
            groups.append(current_group)

        return groups, df_sorted

    def _select_best_from_groups(self, groups, df_sorted):
        filtered_rows = []
        for group in groups:
            group_df = df_sorted.loc[group]
            best_row = group_df.loc[group_df['eval'].idxmin()]
            filtered_rows.append(best_row)
        return pd.DataFrame(filtered_rows)

    def _process(self):
        results = []
        for pid, group in self.df.groupby("target_acc"):
            groups, sorted_group = self._group_overlapping_ranges(group)
            filtered = self._select_best_from_groups(groups, sorted_group)
            results.append(filtered)
        return pd.concat(results, ignore_index=True)

    def __call__(self):
        return self.filtered_df

