import model_global_local

import numpy
import multiprocessing

# This exists primarily to run the code in model_global_local.py
# It's a separate file so we can use multiprocessing -- with a single file there 
# would be some weird bugs.

pool = multiprocessing.Pool()

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

n_origs = [5, 10, 20, 50] + list(range(100, 8000, 100)) + list(range(8000, 20000, 500))
global_local = pool.map(model_global_local.analyse, n_origs)
total_reads = model_global_local.total_reads
for o, (g, l) in zip(n_origs, global_local):
    print(o / total_reads, g / total_reads, l / total_reads)


