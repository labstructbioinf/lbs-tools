from array import array
import json
import os
import re
import numpy as np
from statistics import mean
import pandas as pd
import warnings

warnings.filterwarnings("ignore")


class AF2ScoreParser:
    """
    Parse AF2 / ColabFold outputs for many jobs in a common root directory.

    Expected directory layout (compatible with af2run):
        root_dir/
            job1/
                job1.csv
                job1.done.txt
                job1_scores_rank_001_alphafold2_multimer_v3_model_2_seed_048.json
                job1_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_048.pdb
                ...
            job2/
                job2.csv
                job2.done.txt
                ...

    Usage
    -----
    parser = AF2ScoreParser("/path/to/root_dir")
    df = parser.obtain_data()

    The returned DataFrame has one row per model across all finished jobs.
    """

    def __init__(self, root_dir: str) -> None:
        self.root_dir = root_dir

    # --- helpers for one job -------------------------------------------------

    @staticmethod
    def _get_chain_lengths(job_dir: str, job_name: str):
        """
        Infer per-chain sequence lengths from the colabfold CSV in a job directory.

        Expects CSV of the form:
            id,sequence
            job_name,SEQ1:SEQ2:...

        Returns
        -------
        list[int]
            List of chain lengths (one entry per chain).
        """
        csv_path = os.path.join(job_dir, f"{job_name}.csv")
        if not os.path.isfile(csv_path):
            # fallback: any CSV in this dir
            csv_list = [os.path.join(job_dir, x) for x in os.listdir(job_dir) if x.endswith(".csv")]
            if not csv_list:
                raise RuntimeError(f"No CSV file found in job directory: {job_dir}")
            csv_path = csv_list[0]

        with open(csv_path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            raise RuntimeError(f"CSV file {csv_path} does not contain data lines.")

        data = lines[1].split(",")
        if len(data) < 2:
            raise RuntimeError(f"CSV file {csv_path} does not contain a 'sequence' field.")

        seq_field = data[1].strip()
        if not seq_field:
            raise RuntimeError(f"Empty sequence field in CSV file: {csv_path}")

        seqs = seq_field.split(":")
        chain_lengths = [len(s) for s in seqs]
        return chain_lengths

    @staticmethod
    def _parse_scores_json(filename: str):
        """
        Extract pLDDT (full vector + mean), pTM, ipTM, ranking_confidence
        and PAE array from a scores.json-like file.
        """
        with open(filename) as file:
            parsed = json.load(file)

        plddt_vec = parsed["plddt"]
        plddt_mean = round(mean(plddt_vec), 2)
        ptm = parsed.get("ptm", None)
        iptm = parsed.get("iptm", None)

        # ColabFold / AF2 multimer ranking confidence:
        # ranking_confidence = 0.2 * pTM + 0.8 * ipTM
        ranking_confidence = None
        if ptm is not None and iptm is not None:
            ranking_confidence = 0.2 * ptm + 0.8 * iptm

        pae = np.asarray(parsed["pae"], dtype="float")
        return plddt_vec, plddt_mean, ptm, iptm, ranking_confidence, pae

    @staticmethod
    def _get_model_id_rank_model_seed(filename: str):
        """
        Extract:
          - id    (job/base name),
          - rank  (1..N),
          - model (1..5),
          - seed  (integer).

        Expected pattern, e.g.:
          BC_BC_natural_scores_rank_001_alphafold2_multimer_v3_model_2_seed_048.json
        """
        clean_path = os.path.basename(filename)
        m = re.match(
            r"^(.*)_scores_rank_(\d+)_.*_model_(\d+)_seed_(\d+)\.json$", clean_path
        )
        if not m:
            raise ValueError(f"Unrecognized scores filename format: {clean_path}")
        _id = m.group(1)
        rank = int(m.group(2))
        model = int(m.group(3))
        seed = int(m.group(4))
        return _id, rank, model, seed

    @staticmethod
    def _get_pdb_path_for_scores(job_dir: str, scores_fname: str) -> str:
        """
        Given a scores JSON filename, try to find the corresponding PDB file.

        Typical mapping:
            *_scores_rank_XXX_..._model_Y_seed_Z.json
        ->  *_unrelaxed_rank_XXX_..._model_Y_seed_Z.pdb
            or *_relaxed_rank_XXX_..._model_Y_seed_Z.pdb

        Returns
        -------
        str
            Full path to PDB if found, otherwise None.
        """
        base = os.path.basename(scores_fname)

        # try unrelaxed
        pdb_unrelaxed = base.replace("_scores_", "_unrelaxed_")
        pdb_unrelaxed = os.path.splitext(pdb_unrelaxed)[0] + ".pdb"
        pdb_unrelaxed_path = os.path.join(job_dir, pdb_unrelaxed)
        if os.path.isfile(pdb_unrelaxed_path):
            return pdb_unrelaxed_path

        # try relaxed
        pdb_relaxed = base.replace("_scores_", "_relaxed_")
        pdb_relaxed = os.path.splitext(pdb_relaxed)[0] + ".pdb"
        pdb_relaxed_path = os.path.join(job_dir, pdb_relaxed)
        if os.path.isfile(pdb_relaxed_path):
            return pdb_relaxed_path

        # nothing found
        return None

    @staticmethod
    def _gather_data_for_model(
        pae_arr,
        chain_lengths,
        rank: int,
        model: int,
        seed: int,
        ptm: float,
        iptm: float,
        ranking_confidence: float,
        _id: str,
        plddt_vec,
        plddt_mean: float,
        job_name: str,
        pdb_path: str,
        model_id: str,
    ):
        """
        Compute intra- and inter-chain PAE scores for a single model and
        slice pLDDT per chain.

        Definitions
        -----------
        - intra_pae:
            list of length 'mer' (number of chains),
            intra_pae[i] = mean PAE within chain i (diagonal block).
        - inter_pae:
            mer == 1 (monomer): []
            mer == 2 (dimer)  : [ mean PAE over interface A<->B ]
            mer >= 3          : list of length 'mer', where inter_pae[i] is the
                                mean PAE over all interfaces of chain i with all
                                other chains (average over pairs (i,j), j!=i).

        pLDDT:
            - plddt_per_chain: list of lists, one list per chain (raw pLDDT).

        Returns
        -------
        list
            Single result row with aggregate metrics for this model.
        """
        mer = len(chain_lengths)
        if mer <= 0:
            print(f"***WARNING***, no chains found for job {job_name}!")
            return None

        total_len = pae_arr.shape[0]
        if total_len != sum(chain_lengths):
            raise ValueError(
                f"PAE matrix size ({total_len}) does not match sum of chain lengths "
                f"({sum(chain_lengths)}) in job {job_name}."
            )
        if len(plddt_vec) != total_len:
            raise ValueError(
                f"pLDDT length ({len(plddt_vec)}) does not match PAE size ({total_len}) in job {job_name}."
            )

        # chain boundaries: [0, L1, L1+L2, ...]
        boundaries = np.cumsum([0] + chain_lengths)

        # pLDDT per chain
        plddt_per_chain = []
        for i in range(mer):
            s, e = boundaries[i], boundaries[i + 1]
            plddt_per_chain.append(plddt_vec[s:e])

        # build block matrix of PAE submatrices [mer x mer]
        blocks = []
        for i in range(mer):
            row_blocks = []
            r_start, r_end = boundaries[i], boundaries[i + 1]
            for j in range(mer):
                c_start, c_end = boundaries[j], boundaries[j + 1]
                row_blocks.append(pae_arr[r_start:r_end, c_start:c_end])
            blocks.append(row_blocks)

        # intra per chain
        intra_list = []
        for i in range(mer):
            intra_block = blocks[i][i]
            intra_list.append(float(intra_block.mean()))

        # inter: najpierw pary (i<j)
        pair_means = np.full((mer, mer), np.nan)
        for i in range(mer):
            for j in range(i + 1, mer):
                pair_vals = np.concatenate(
                    [blocks[i][j].ravel(), blocks[j][i].ravel()]
                )
                m = float(pair_vals.mean())
                pair_means[i, j] = m
                pair_means[j, i] = m

        # inter list wg specyfikacji
        if mer == 1:
            inter_list = []
        elif mer == 2:
            # jedna para A-B
            inter_list = [float(pair_means[0, 1])]
        else:
            # dla każdego łańcucha średnia po wszystkich interfejsach z innymi
            inter_list = []
            for i in range(mer):
                vals = pair_means[i, :]
                vals = vals[~np.isnan(vals)]
                vals = vals[vals != 0]
                if len(vals) == 0:
                    inter_list.append(float("nan"))
                else:
                    inter_list.append(float(vals.mean()))

        intra_mean = float(mean(intra_list)) if len(intra_list) > 0 else float("nan")
        if mer == 1:
            inter_mean = float("nan")
        elif mer == 2:
            inter_mean = float(pair_means[0, 1])
        else:
            inter_mean = float(mean(inter_list))

        return [
            job_name,
            pdb_path,
            model_id,
            _id,
            mer,
            intra_list,
            inter_list,
            intra_mean,
            inter_mean,
            ptm,
            iptm,
            ranking_confidence,
            int(model),
            int(rank),
            int(seed),
            plddt_mean,
            plddt_per_chain,
            chain_lengths,   # NEW
            pae_arr,         # NEW
        ]

    # --- main API ------------------------------------------------------------

    def obtain_data(
        self,
        force_foldseek: bool = False,
        tmscore_threshold: float = 0.95,
        limit_job_names: list = [],
    ) -> pd.DataFrame:
        """
        Scan the root directory for finished jobs and collect metrics
        for all models in those jobs.
    
        A "finished" job is defined as a subdirectory whose name is `job_name`
        and which contains a file `{job_name}.done.txt`.
        """
        base_cols = [
            "job_name",
            "pdb_path",
            "model_id",
            "id",
            "oligomerisation_state",
            "intra_pae",
            "inter_pae",
            "intra_pae_mean",
            "inter_pae_mean",
            "ptm",
            "iptm",
            "ranking_confidence",
            "model",
            "model_rank",
            "seed",
            "plddt_mean",
            "plddt_per_chain",
            "chain_lengths",   # NEW
            "pae",             # NEW: pełna macierz PAE jako np.ndarray
        ]

        # global dataframe, na starcie z samymi bazowymi kolumnami
        dataframe = pd.DataFrame(columns=base_cols)

        # parse all jobs or only selected ones
        job_names = os.listdir(self.root_dir)
        if limit_job_names != []:
            job_names = list(set(job_names) & set(limit_job_names))
        
        # iterate over job subdirectories in root_dir
        for job_name in job_names:
    
            job_dir = os.path.join(self.root_dir, job_name)
            if not os.path.isdir(job_dir):
                continue
    
            done_flag = os.path.join(job_dir, f"{job_name}.done.txt")
            if not os.path.isfile(done_flag):
                continue
    
            print(f'parsing scores for {job_name}')
    
            chain_lengths = self._get_chain_lengths(job_dir, job_name)
    
            # ---------- PRE DF (run-local) ----------------------------------
    
            pre_rows = []
    
            for fname in os.listdir(job_dir):
                if not (fname.endswith(".json") and "scores" in fname):
                    continue
    
                scores_path = os.path.join(job_dir, fname)
    
                _id, rank, model, seed = self._get_model_id_rank_model_seed(scores_path)
    
                (
                    plddt_vec,
                    plddt_mean,
                    ptm,
                    iptm,
                    ranking_confidence,
                    pae_array,
                ) = self._parse_scores_json(scores_path)
    
                pdb_path = self._get_pdb_path_for_scores(job_dir, fname)
                model_id = os.path.splitext(os.path.basename(pdb_path))[0]
    
                row = self._gather_data_for_model(
                    pae_arr=pae_array,
                    chain_lengths=chain_lengths,
                    rank=rank,
                    model=model,
                    seed=seed,
                    ptm=ptm,
                    iptm=iptm,
                    ranking_confidence=ranking_confidence,
                    _id=_id,
                    plddt_vec=plddt_vec,
                    plddt_mean=plddt_mean,
                    job_name=job_name,
                    pdb_path=pdb_path,
                    model_id=model_id
                )
    
                if row is not None:
                    pre_rows.append(row)
    
            # jeśli nic nie ma, przechodzimy dalej
            if not pre_rows:
                continue
    
            pre_df = pd.DataFrame(pre_rows, columns=base_cols)
    
            # ---------- CLUSTERING -------------------------------------------
    
            print(f'clustering models for {job_name}')
    
            pdb_paths = pre_df["pdb_path"].dropna().tolist()
            if len(pdb_paths) == 0:
                warnings.warn(f"No PDBs for job {job_name}; skipping clustering.")
                dataframe = pd.concat([dataframe, pre_df], ignore_index=True)
                continue
    
            fs_work = os.path.join(job_dir, "_foldseek_cluster")
            sel_dir = os.path.join(fs_work, "selected")
            tmp_dir = "/tmp"
            out_prefix = os.path.join(fs_work, "af2multimer")
            clu_file = out_prefix + "_cluster.tsv"
    
            # czy uruchomić foldseeka?
            run_foldseek = force_foldseek or (not os.path.isfile(clu_file))
    
            if run_foldseek:
                os.makedirs(sel_dir, exist_ok=True)
                os.makedirs(tmp_dir, exist_ok=True)
    
                # symlink only PDBs from THIS job
                for p in pdb_paths:
                    if p and os.path.isfile(p):
                        dst = os.path.join(sel_dir, os.path.basename(p))
                        if not os.path.exists(dst):
                            os.symlink(os.path.abspath(p), dst)
    
                foldseek_bin = os.path.expanduser("~/apps/foldseek/bin/foldseek")
    
                cmd = [
                    foldseek_bin,
                    "easy-multimercluster",
                    sel_dir,
                    out_prefix,
                    tmp_dir,
                    "-c", "0.8",
                    "--cov-mode", "0",
                    "--tmscore-threshold", str(tmscore_threshold),
                    "--threads", "10",
                    " > /dev/null",
                ]
    
                print(" ".join(cmd))
    
                ret = os.system(" ".join(cmd))
                if ret != 0:
                    raise RuntimeError(f"Foldseek failed for job {job_name}")
    
            else:
                print(f"using existing Foldseek clusters for {job_name}: {clu_file}")
    
            if not os.path.isfile(clu_file):
                raise RuntimeError(f"Foldseek output missing: {clu_file}")
    
            # parse clusters
            clu_map = {}
            with open(clu_file) as f:
                for l in f:
                    rep, mem = l.strip().split("\t")
                    clu_map[mem] = rep
    
            # encode cluster names as unique numbers
            reps = sorted(set(clu_map.values()))
            rep_to_id = {rep: i for i, rep in enumerate(reps, 1)}
    
            # annotate PRE DF
            pre_df["cluster_rep"] = pre_df["pdb_path"].apply(
                lambda p: clu_map.get(os.path.splitext(os.path.basename(p))[0], None) if p else None
            )
            pre_df["cluster_id"] = pre_df["cluster_rep"].apply(
                lambda r: rep_to_id.get(r, None)
            )
    
            # ---------- append to GLOBAL DATAFRAME ----------------------------
    
            dataframe = pd.concat([dataframe, pre_df], ignore_index=True)
    
            # debug
            #break

        assert dataframe.model_id.is_unique

        # mark clusters' representatives
        if len(dataframe)>0:
            dataframe['is_rep'] = dataframe['model_id'] == dataframe['cluster_rep']
        
        return dataframe


