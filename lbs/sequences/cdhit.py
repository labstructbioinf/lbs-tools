import os
import tempfile
import subprocess
from Bio import SeqIO

class cdhit:
    def __init__(self, sequences, cdhitbin="", cpus=1, debug=False):
        if not os.path.isfile(cdhitbin):
            raise ValueError(f"CD-HIT binary not found: {cdhitbin}")
        self.cdhitbin = cdhitbin
        self.sequences = sequences
        self.cpus = cpus
        self.debug = debug

    def run(self, identity=0.9, coverage=0.0, maxlendiff=0.0):
        if not (0.4 <= identity <= 1.0):
            raise ValueError("identity must be between 0.4 and 1.0")

        # determine word size and minimal length required
        if identity > 0.7:
            word = 5
            min_len = 11
        elif identity > 0.6:
            word = 4
            min_len = 9
        elif identity > 0.5:
            word = 3
            min_len = 7
        else:
            word = 2
            min_len = 5

        shortest_seq_len = min(len(s.seq) for s in self.sequences)
        if shortest_seq_len < min_len:
            raise ValueError(
                f"Cannot cluster: identity={identity} requires sequences ≥{min_len} aa, "
                f"but the shortest is {shortest_seq_len} aa. Try lowering identity or filtering short sequences."
            )

        if shortest_seq_len < 11:
            print(
                f"[CD-HIT NOTICE] Using short sequences (shortest={shortest_seq_len}). "
                f"Word length set to {word}, enforcing min seq len ≥ {min_len}"
            )

        seen_ids = set()

        with tempfile.NamedTemporaryFile(mode='wt', delete=False, suffix=".fasta") as inf:
            for s in self.sequences:
                ident = s.id
                if '...' in ident or '>' in ident:
                    raise ValueError(f"Invalid sequence ID: {ident}")
                if ident in seen_ids:
                    raise ValueError(f"Duplicate sequence ID: {ident}")
                seen_ids.add(ident)
                inf.write(s.format('fasta'))
            input_path = inf.name

        maxdesclen = max(len(s.id) for s in self.sequences)
        output_path = tempfile.mktemp(suffix=".fasta")
        clstr_file = output_path + ".clstr"

        if self.debug:
            print(f"[CD-HIT DEBUG] input: {input_path}")
            print(f"[CD-HIT DEBUG] output: {output_path}")
            print(f"[CD-HIT DEBUG] cluster file: {clstr_file}")

        cmd = [
            self.cdhitbin,
            "-i", input_path,
            "-o", output_path,
            "-d", str(maxdesclen + 1),
            "-T", str(self.cpus),
            "-c", str(identity),
            "-n", str(word),
            "-A", str(coverage),
            "-s", str(maxlendiff),
            "-l", str(min_len)
        ]

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            raise RuntimeError(f"CD-HIT failed. Check input: {input_path} and output: {output_path}")

        # parse clustering result
        clusters = {}
        current_cluster = None
        with open(clstr_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">Cluster"):
                    current_cluster = int(line.split()[-1])
                    clusters[current_cluster] = []
                elif line:
                    start = line.find('>') + 1
                    end = line.find('...')
                    seq_id = line[start:end]
                    clusters[current_cluster].append(seq_id)

        # cleanup (only if not in debug mode)
        if not self.debug:
            os.remove(input_path)
            os.remove(output_path)
            os.remove(clstr_file)

        return clusters
