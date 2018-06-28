import model_global_local

import numpy

# This is part of a series of scripts which models sequencing reads and aims
# to determine how frequently randomly positioned reads are identified as 
# co-located reads, by an algorithm equivalent to "suprDUPr".

# Computes duplicates under an alternative distribution of reads, in
# which they are sampled from a collection of "sub-libraries". These
# sub-libraries could represent e.g. transcripts in RNA-seq or amplified
# genes in amplicons.

print("-- Duplicate count toy model --")

print("")
print("Coordinate ranges:")
print("  x: ", model_global_local.x_lim_tot)
print("  y: ", model_global_local.y_lim_tot)
print("Number of tiles:", model_global_local.n_tiles)
print("")

print("Distance thresholds for local duplicates:")
print("  x: ", model_global_local.x_lim_dist)
print("  y: ", model_global_local.x_lim_dist)
print("")

# Number of possible substrings inside a sub-library. Note that results
# depend on the number of reads simulated.

models = [
        ([5000, 5000, 5000, 5000], [0.005, 0.004, 0.003, 0.002]),
        ([5000, 5000, 5000, 5000], [0.01, 0.01, 0.008, 0.008]),
        ([5000, 5000, 5000, 5000], [0.1, 0.05, 0.008, 0.008]),
        ([5000], [0.2]),
        ]

global_local = map(model_global_local.analyse_with_sublibraries, models)
total_reads = model_global_local.total_reads
print("[Sublibrary_Lengths] [Sublibrary_Fractions] global_dup local_dups")
for o, (g, l) in zip(models, global_local):
    print("["+", ".join("%3.1f" % oi for oi in o[0]) + "]", o[1], g / total_reads, l / total_reads)

