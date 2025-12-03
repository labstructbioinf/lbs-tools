import os
import re
  
class af2run:
    """
    Helper class to prepare input CSV files and a Slurm array script
    for running colabfold_batch on a cluster.

    sequences_dict:
        dict of the form:
        {
            "job1": ["SEQ_A"],                    # monomer
            "job2": ["SEQ_A", "SEQ_B"],          # dimer / PPI
            "job3": ["SEQ1", "SEQ2", "SEQ3"],    # 3-chain complex
            ...
        }

    jobs_dir:
        directory where subdirectories for each job and CSV files are created.

    name:
        Slurm job name (used in -J).
    """

    def is_valid_slurm_job_name(self, name: str) -> bool:
        return bool(re.fullmatch(r"[A-Za-z0-9\-_]+", name)) and len(name) <= 128
    
    def __init__(self, sequences_dict, jobs_dir, name, colabfold_config=None):
        self.sequences_dict = sequences_dict
        self.jobs_dir = jobs_dir
        self.name = name
        assert self.is_valid_slurm_job_name(self.name)
        
        # Optional configuration for colabfold_batch CLI flags
        self.colabfold_config = colabfold_config or {}
    
    def _build_colabfold_parameters(self, extra_parameters: str = "") -> str:
        """
        Convert self.colabfold_config into a CLI parameter string
        for colabfold_batch.

        Supported keys:
          - stop_at_score -> --stop-at-score
          - num_recycle -> --num-recycle
          - recycle_early_stop_tolerance -> --recycle-early-stop-tolerance
          - num_ensemble -> --num-ensemble
          - num_seeds -> --num-seeds
          - random_seed -> --random-seed
          - num_models -> --num-models
          - recompile_padding -> --recompile-padding
          - model_order -> --model-order
          - host_url -> --host-url
          - data -> --data
          - msa_mode -> --msa-mode
          - model_type -> --model-type
          - amber (bool) -> --amber
          - num_relax -> --num-relax
          - templates (bool) -> --templates
          - custom_template_path -> --custom-template-path
          - rank -> --rank
          - pair_mode -> --pair-mode
          - sort_queries_by -> --sort-queries-by
          - save_single_representations (bool) -> --save-single-representations
          - save_pair_representations (bool) -> --save-pair-representations
          - use_dropout (bool) -> --use-dropout
          - max_seq -> --max-seq
          - max_extra_seq -> --max-extra-seq
          - max_msa -> --max-msa  (special: expects 'max_seq:max_extra_seq')
          - disable_cluster_profile (bool) -> --disable-cluster-profile
          - zip_results (bool) -> --zip
          - use_gpu_relax (bool) -> --use-gpu-relax
          - save_all (bool) -> --save-all
          - save_recycles (bool) -> --save-recycles
          - overwrite_existing_results (bool) -> --overwrite-existing-results
          - disable_unified_memory (bool) -> --disable-unified-memory

        Note on max_msa:
          - If 'max_msa' is provided, it should be either:
              * a string "512:1024", or
              * a 2-element list/tuple [512, 1024] / (512, 1024)
            and in that case --max-msa is emitted and max_seq / max_extra_seq
            from the config are ignored.
          - If 'max_msa' is not provided, max_seq and max_extra_seq are used
            as separate flags (if present).
        """
        flag_map = {
            "stop_at_score": "--stop-at-score",
            "num_recycle": "--num-recycle",
            "recycle_early_stop_tolerance": "--recycle-early-stop-tolerance",
            "num_ensemble": "--num-ensemble",
            "num_seeds": "--num-seeds",
            "random_seed": "--random-seed",
            "num_models": "--num-models",
            "recompile_padding": "--recompile-padding",
            "model_order": "--model-order",
            "host_url": "--host-url",
            "data": "--data",
            "msa_mode": "--msa-mode",
            "model_type": "--model-type",
            "amber": "--amber",
            "num_relax": "--num-relax",
            "templates": "--templates",
            "custom_template_path": "--custom-template-path",
            "rank": "--rank",
            "pair_mode": "--pair-mode",
            "sort_queries_by": "--sort-queries-by",
            "save_single_representations": "--save-single-representations",
            "save_pair_representations": "--save-pair-representations",
            "use_dropout": "--use-dropout",
            "max_seq": "--max-seq",
            "max_extra_seq": "--max-extra-seq",
            # "max_msa" handled specially below
            "disable_cluster_profile": "--disable-cluster-profile",
            "zip_results": "--zip",
            "use_gpu_relax": "--use-gpu-relax",
            "save_all": "--save-all",
            "save_recycles": "--save-recycles",
            "overwrite_existing_results": "--overwrite-existing-results",
            "disable_unified_memory": "--disable-unified-memory",
        }

        parts = []

        # Special handling for max_msa
        use_max_msa = False
        if "max_msa" in self.colabfold_config and self.colabfold_config["max_msa"] is not None:
            v = self.colabfold_config["max_msa"]
            if isinstance(v, str):
                if ":" not in v:
                    raise ValueError("max_msa must be in the form 'max_seq:max_extra_seq', e.g. '512:1024'.")
                max_msa_str = v
            elif isinstance(v, (list, tuple)) and len(v) == 2:
                max_msa_str = f"{int(v[0])}:{int(v[1])}"
            else:
                raise ValueError("max_msa must be a string 'max_seq:max_extra_seq' or a 2-element list/tuple.")
            parts.append(f"--max-msa {max_msa_str}")
            use_max_msa = True

        for key, flag in flag_map.items():
            # skip max_seq / max_extra_seq if max_msa is explicitly used
            if use_max_msa and key in ("max_seq", "max_extra_seq"):
                continue

            if key not in self.colabfold_config:
                continue
            value = self.colabfold_config[key]
            if value is None:
                continue

            if isinstance(value, bool):
                if value:
                    parts.append(flag)
            else:
                parts.append(f"{flag} {value}")

        extra_parameters = extra_parameters.strip()
        if extra_parameters:
            parts.append(extra_parameters)

        return " ".join(parts)
        
    def prepare_input(self):
        """
        For each job, create a directory and a CSV file:

        id,sequence
        job_name,SEQ1:SEQ2:...

        where SEQ1, SEQ2, ... are chains in the complex (1 for monomer, 2+ for PPI).
        """
        os.makedirs(self.jobs_dir, exist_ok=True)
        self.model_paths = {}
        
        for job_name, seq_list in self.sequences_dict.items():
            model_path = os.path.join(self.jobs_dir, job_name)
            os.makedirs(model_path, exist_ok=True)
            self.model_paths[job_name] = model_path
            
            sequence_field = ":".join(seq_list)
            csv_path = os.path.join(model_path, job_name + '.csv')
            txt = f"id,sequence\n{job_name},{sequence_field}\n"
            with open(csv_path, 'wt') as f:
                f.write(txt)
                
    def prepare_array_job(self, jobs_in_parallel: int, parameters: str = ""):
        """
        Create a Slurm array script that, for each job, runs:

        colabfold_batch {parameters} job_name.csv .
        """
        parameters = self._build_colabfold_parameters(parameters)

        todo = []
        for job_name, model_path in self.model_paths.items():
            done_flag = os.path.join(model_path, f'{job_name}.done.txt')
            if not os.path.isfile(done_flag):
                todo.append(job_name) 

        run_script_path = os.path.join(self.jobs_dir, 'run_array.sh')

        if not todo:
            print("No jobs to process.")
            with open(run_script_path, 'w') as run_script:
                run_script.write('# All requested jobs done.\n')
            return

        with open(run_script_path, 'w') as run_script:
            run_script.write(f"""#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 4
#SBATCH --gres=gpu:1
#SBATCH --exclude=edi00
#SBATCH --mem=16GB
#SBATCH -J {self.name}
#SBATCH --array=1-{len(todo)}%{jobs_in_parallel}

directories=(
""")
            for t in todo:
                run_script.write(f'  "{t}"\n')

            run_script.write(f""")

sleep $((RANDOM % 60 + 1))
cd {self.jobs_dir}/${{directories[$SLURM_ARRAY_TASK_ID-1]}}
source /opt/miniconda3/bin/activate cf_1.5
export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64:$LD_LIBRARY_PATH
colabfold_batch {parameters} ${{directories[$SLURM_ARRAY_TASK_ID-1]}}.csv .
""")
