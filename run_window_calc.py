import model_global_local

import numpy
import multiprocessing

pool = multiprocessing.Pool()

# Determine the ratio of reads within search window, to the total number of
# reads. This is approximately the ratio of areas. It is also affected by reads
# on the edge of the tile, which don't have a full search area.

NTEST = 10000
print("Computing the search window ratio using", NTEST, "samples")
# Using a parallel pool here just doesn't work, no matter what
num_in = pool.map(model_global_local.eval_in_range, (tuple(c) for c in model_global_local.coords[0:NTEST]))
print("In range percentage:", sum(num_in) * 100.0/ (len(model_global_local.coords) * NTEST), "%")

