# EXERCISE 1

import joblib as jl
import subprocess as sbp

GENOME_PATH = "/path/to/human.fsa"
SPLITS_DIR = "/path/to/human_splits"  # Files are here, named as split_1.fsa, split_2.fsa, etc.
OUT_DIR = "/path/to/revcomp_splits"
N_SEQUENCES = 24


def sys_call(cmd):  # Send command to system
    job = sbp.run(cmd.split(), stdout=sbp.PIPE, universal_newlines=True)
    result = job.stdout.split("\n")
    return result


# Generate commands for each split
jobs = [f"python3 worker.py {SPLITS_DIR}/split_{i + 1}.fsa {OUT_DIR}/revcomp_split_{i + 1}.fsa"
        for i in range(N_SEQUENCES)]

result = jl.Parallel(n_jobs=8)(jl.delayed(sys_call)(job) for job in jobs)  # Run each job in list
