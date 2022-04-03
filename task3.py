# Exercise 3

import joblib as jl
import re

# INDEXING #

GENOME_COPY = "/path/to/human.fsa"
CHUNK_SIZE = 2 ** 25

infile = open(GENOME_COPY, "rb")

# Indexing is identical to the previous one
position_before_chunk = 0
chunk = infile.read(CHUNK_SIZE)
header_starts = []
while True:
    first = False
    matches = re.finditer(b">", chunk)
    positions = [pos.start() for pos in matches]
    if len(positions) > 0:
        header_starts.extend([pos + position_before_chunk for pos in positions])
    position_before_chunk += len(chunk)
    if len(chunk) < CHUNK_SIZE:
        break
    chunk = infile.read(CHUNK_SIZE)

seq_ends = [pos - 1 for pos in header_starts][1:]
seq_ends.append(position_before_chunk - 1)

header_ends = []
seq_starts = []
for h_start in header_starts:
    infile.seek(h_start)
    chunk = infile.read(300)
    nl_pos = chunk.find(b"\n")
    header_ends.append(h_start + nl_pos)
    seq_starts.append(h_start + nl_pos + 1)

indices = [[header_starts[i], header_ends[i], seq_starts[i], seq_ends[i]] for i in range(len(seq_starts))]

print(*indices, sep="\n")

# PARALLELISING #

N_SEQUENCES = 24


# Worker function to be calles with joblib
def worker_function(indices, n_seq):
    genome_file = open(GENOME_COPY, "r+b")  # Open file in read+write binary mode
    genome_file.seek(indices[2])  # Go to sequence position
    out_seq = genome_file.read(indices[3] - indices[2] - 1)  # Read chunk of same size as sequence, -1 to remove newline

    translation_table = bytes.maketrans(b"atcg", b"tagc")
    out_seq = out_seq.translate(translation_table)  # Complementary
    out_seq = out_seq[::-1]  # Reverse
    new_header = b">s" + bytes(f"{n_seq + 1}", "utf-8") + \
                 b" A: " + bytes(str(out_seq.count(b"a")), "utf-8") + \
                 b" C: " + bytes(str(out_seq.count(b"c")), "utf-8") + \
                 b" T: " + bytes(str(out_seq.count(b"t")), "utf-8") + \
                 b" G: " + bytes(str(out_seq.count(b"g")), "utf-8") + \
                 b" N: " + bytes(str(out_seq.count(b"n")), "utf-8") + b"\n"  # Format header

    header_size = indices[2] - indices[0]  # Calculate size of header

    # Ensure new header is exactly the same size, to not cause the rest of the file to dislodge
    if len(new_header) > header_size:
        new_header = new_header[:header_size - 1] + b"\n"
    elif len(new_header) < header_size:
        new_header = new_header[:-1] + b" " * (header_size - len(new_header)) + b"\n"

    genome_file.seek(indices[0])  # Go to header position
    genome_file.write(new_header) # Replace header
    genome_file.seek(indices[2])  # Go to sequence position
    genome_file.write(out_seq)    # Replace sequence
    genome_file.close()


# Generate list of tuples with the arguments for each call of the worker function
args = [(seq_indices, i) for i, seq_indices in enumerate(indices)]

# Spawn processes
result = jl.Parallel(n_jobs=8)(jl.delayed(worker_function)(args_set[0], args_set[1]) for args_set in args)
print(*result, sep="\n")  # Print any messages

