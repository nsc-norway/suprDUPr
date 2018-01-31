# Toy model for determining the relationship between local and global duplicates

import numpy

print("-- Duplicate count toy model --")

# Area parameters
x_lim_tot = (1100, 33000)
y_lim_tot = (1000, 50000)

# Max distance along the axis (so, the window size is about twice the *_lim_dist)
x_lim_dist = (2500, 2500)
y_lim_dist = (2500, 2500)
n_tiles = 112


# Simulate PCR duplicates -- all reads generated from a limited space of original sequences
# nseq_orig sequences amplified up to nseq_reads. The sequences are just represented as 
# integers.
nseq_orig = int(1e6)
nseq_reads = int(10e6)


# Generate nseq_reads by sampling randomly from the space of [0,nseq_orig].
reads = numpy.random.random_integers(1, nseq_orig, (nseq_reads))
xs = numpy.random.random_integers(*x_lim_tot, (nseq_reads))
ys = numpy.random.random_integers(*x_lim_tot, (nseq_reads))
tiles = numpy.random.random_integers(1, n_tiles, (nseq_reads))


# ** Global duplication ratio: count number of duplicates for each unique sequence, regardless of position**
# First get unique sequences (numbers) and the count of each. Only the count is needed for the global duplication ratio.
unique_seqs, indexes, counts = numpy.unique(reads, return_inverse=True, return_counts=True)

num_global_duplicates = (counts - 1).sum()

print("Global duplication rate: {0:.1f} %".format(num_global_duplicates * 100.0 / nseq_reads))

# ** Local duplication ratio **
# Count the number of reads within a threshold distance, and in the same tile

local_duplicate_counter = 0

for i in range(len(unique_seqs)):
    indexes = numpy.argwhere(indexes == i)
    xs, ys, tiles = [], [], []
    for j, cand in enumerate(indexes):
        pass


print("Local duplication rate: {0:.1f} %".format(local_duplicate_counter * 100.0 / nseq_reads))

