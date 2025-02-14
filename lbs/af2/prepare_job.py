import os
import re
  
class af2run:
    
    def is_valid_slurm_job_name(self, name):
        return bool(re.fullmatch(r"[A-Za-z0-9\-_]+", name)) and len(name) <= 128
    
    def __init__(self, sequences_dict, jobs_dir, rep, name):
        self.sequences_dict = sequences_dict
        self.jobs_dir = jobs_dir
        self.rep = rep
        self.name = name
        assert self.is_valid_slurm_job_name(self.name)
        
    def prepare_input(self):
        os.makedirs(self.jobs_dir, exist_ok=True)
        
        self.model_paths = {}
        
        for i, seq in self.sequences_dict.items():
            model_path = os.path.join(self.jobs_dir, i)
            os.makedirs(model_path, exist_ok=True)
            
            self.model_paths[i] = model_path
            
            seq = [seq] * self.rep
            afin_path = os.path.join(model_path, i+'.csv')
            txt=f"id,sequence\n{i},{':'.join(seq)}\n"
            with open(afin_path, 'wt') as f:
                f.write(txt)   
                
    def prepare_array_job(self, jobs_in_parallel):
        todo = []

        for i, model_path in self.model_paths.items():
            done_flag = os.path.join(model_path, f'{i}.done.txt')
            if not os.path.isfile(done_flag):
                todo.append(i) 

        if not todo:
            print("No jobs to process.")
            return

        run_script_path = os.path.join(self.jobs_dir, 'run_array.sh')

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

cd {self.jobs_dir}/${{directories[$SLURM_ARRAY_TASK_ID-1]}}
source /opt/miniconda3/bin/activate cf_1.5
export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64:$LD_LIBRARY_PATH
colabfold_batch ${{directories[$SLURM_ARRAY_TASK_ID-1]}}.csv .
""")


