import model_global_local

import numpy
import multiprocessing

pool = multiprocessing.Pool()

# Computes duplicates under an alternative distribution of reads, in
# which they are sampled from a collection of "sub-libraries". These
# sub-libraries could represent e.g. transcripts in RNA-seq or amplified
# genes in amplicons.

num_in = pool.map(model_global_local.eval_in_range, (tuple(c) for c in model_global_local.coords[0:NTEST]))
print("In range percentage:", sum(num_in) * 100.0/ (len(model_global_local.coords) * NTEST), "%")

