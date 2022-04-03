# Worker script for exercise 1

import sys
import os

# Receive paths from args
INFILE_PATH = sys.argv[1]
OUTFILE_PATH = sys.argv[2]

infile = open(INFILE_PATH, "rb")
size = os.path.getsize(INFILE_PATH)
out_seq = infile.read(size)
out_seq = out_seq[out_seq.find(b"\n"):]

# Use translation table for complementary
translationTable = bytes.maketrans(b"atcg", b"tagc")
out_seq = out_seq.translate(translationTable)

out_seq = out_seq[::-1]  # Reverse

# Format header with sequence counts
out_seq_header = b">seq" + bytes(f"{sys.argv[3]}", "utf-8") + \
                 b" A: " + bytes(str(out_seq.count(b"a")), "utf-8") + \
                 b" C: " + bytes(str(out_seq.count(b"c")), "utf-8") + \
                 b" T: " + bytes(str(out_seq.count(b"t")), "utf-8") + \
                 b" G: " + bytes(str(out_seq.count(b"g")), "utf-8") + \
                 b" N: " + bytes(str(out_seq.count(b"n")), "utf-8") + b"\n"

infile.close()

# Write out file
outfile = open(OUTFILE_PATH, "wb")
outfile.write(out_seq_header)
outfile.write(out_seq)
outfile.close()
