import sys
import time

# Receive paths from args
INFILE_PATH = sys.argv[1]
OUTFILE_PATH = sys.argv[2]

out_seq = []

infile = open(INFILE_PATH, "rb")
for line in infile:
    if line[0].to_bytes(1, "little") == b">":
        out_seq_header = line  # Save header line
    else:
        out_seq.append(line)  # Add sequence lines to list

out_seq = b"".join(out_seq)  # Join list

# Use translation table for complementary
translationTable = bytes.maketrans(b"atcg", b"tagc")
out_seq = out_seq.translate(translationTable)

out_seq = out_seq[::-1]  # Reverse

# Format header with sequence counts
out_seq_header = b">seq" + bytes(f"{sys.argv[1]}", "utf-8") + \
                 b" A: " + bytes(str(out_seq.count(b"a")), "utf-8") + \
                 b" C: " + bytes(str(out_seq.count(b"c")), "utf-8") + \
                 b" T: " + bytes(str(out_seq.count(b"t")), "utf-8") + \
                 b" G: " + bytes(str(out_seq.count(b"g")), "utf-8") + \
                 b" N: " + bytes(str(out_seq.count(b"n")), "utf-8")

infile.close()

# Write out file
outfile = open(OUTFILE_PATH, "wb")
outfile.write(out_seq_header)
outfile.write(b"\n")
outfile.write(out_seq)
outfile.close()
