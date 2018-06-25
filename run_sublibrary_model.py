import model_global_local

import numpy
import multiprocessing

pool = multiprocessing.Pool()

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

# Number of possible substrings inside a sub-library. In general we want total_reads
# to just decide the precision, so the simulated coverage should stay constant. The 
# number of substrings should be expressed as a function of total_reads. As an 
# arbitrary selection, let's look for a 5000 bp transcript in a human transcriptome 
# sequenced at equivalent of 30X coverage.

slf = (1 / (3e9 * 30)) * model_global_local.total_reads

models = [
        ([5000*slf, 5000*slf, 5000*slf, 5000*slf], [0.005, 0.004, 0.003, 0.002]),
        ([5000*slf, 5000*slf, 5000*slf, 5000*slf], [0.01, 0.01, 0.008, 0.008]),
        ([5000*slf, 5000*slf, 5000*slf, 5000*slf], [0.1, 0.05, 0.008, 0.008]),
        ]

global_local = pool.map(model_global_local.analyse_with_sublibraries, models)
total_reads = model_global_local.total_reads
for o, (g, l) in zip(models, global_local):
    print("["+", ".join("%3.1f" % oi for oi in o[0]) + "]", o[1], g / total_reads, l / total_reads)

