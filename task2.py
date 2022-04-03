# Exercise 2

import joblib as jl
import re

# INDEXING #

INFILE = "/mnt/c/Users/dani_/Downloads/HPCweek4_in/human.fsa"
CHUNKSIZE = 2 ** 25

infile = open(INFILE, "rb")

position_before_chunk = 0  # Save position before reading the chunk
chunk = infile.read(CHUNKSIZE)
header_starts = []
while True:
    first = False
    matches = re.finditer(b">", chunk)  # Find init positions of fasta headers
    positions = [pos.start() for pos in matches]  # Convert iterator into list of positions
    if len(positions) > 0:
        # Save found positions as the position in the current chunk plus all previously read chunks
        header_starts.extend([pos + position_before_chunk for pos in positions])
    position_before_chunk += len(chunk)  # Add current chunk length to cumulative position
    if len(chunk) < CHUNKSIZE:
        break  # Break out of the while loop if finished reading in entire file
    chunk = infile.read(CHUNKSIZE)

# Sequence ends will be the previous position to the header starts (to account for the newline)
seq_ends = [pos - 1 for pos in header_starts][1:]  # Delete first element (will always be -1, before the first header)
seq_ends.append(position_before_chunk - 1)  # Add ending position (-1 to account for newline)

# Header ends and sequence starts are a bit trickier. I read a 300 char chunk starting from the header start and find
# the position of the first newline. The header end will be the position of the newline in the chunk plus the position
# from which I started reading the chunk. The sequence start will be that plus 1 (skipping the newline).
header_ends = []
seq_starts = []
for h_start in header_starts:
    infile.seek(h_start)
    chunk = infile.read(300)
    nl_pos = chunk.find(b"\n")
    header_ends.append(h_start + nl_pos)
    seq_starts.append(h_start + nl_pos + 1)

# Save indices to a list of tuples
indices = [(header_starts[i], header_ends[i], seq_starts[i], seq_ends[i]) for i in range(len(seq_starts))]

# PARALLELISING #

N_SEQUENCES = 24
OUT_DIR = "/path/to/revcomp_splits"  # Splits will be here


# Worker function to call with joblib.Parallel(). Almost identical to the worker program.
def worker_function(indices, n_seq, outfile):
    with open(INFILE, "rb") as genome:
        genome.seek(indices[2])
        seq = genome.read(indices[3] - indices[2])

    translation_table = bytes.maketrans(b"atcg", b"tagc")
    out_seq = seq.translate(translation_table)  # Complementary
    out_seq = out_seq[::-1]  # Reverse
    header = b">seq" + bytes(f"{n_seq}", "utf-8") + \
             b" A: " + bytes(str(out_seq.count(b"a")), "utf-8") + \
             b" C: " + bytes(str(out_seq.count(b"c")), "utf-8") + \
             b" T: " + bytes(str(out_seq.count(b"t")), "utf-8") + \
             b" G: " + bytes(str(out_seq.count(b"g")), "utf-8") + \
             b" N: " + bytes(str(out_seq.count(b"n")), "utf-8")  # Format header
    with open(outfile, "wb") as out:
        out.write(header)
        out.write(b"\n")
        out.write(out_seq)


# Generate list of tuples with the arguments for each call of the worker function
args = [(seq_indices, i, f"{OUT_DIR}/revcomp_split_{i + 1}.fsa") for i, seq_indices in enumerate(indices)]

# Spawn processes
result = jl.Parallel(n_jobs=8)(jl.delayed(worker_function)(args_set[0], args_set[1], args_set[2]) for args_set in args)
print(*result, sep="\n")  # Print messages (will only do so with errors)
