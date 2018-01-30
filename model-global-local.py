# Model for determining the relationship between local and global duplicates

import random

# Area parameters
x_lim_tot = (1100, 33000)
y_lim_tot = (1000, 50000)

x_win_tot = (5000, 5000)
y_win_tot = (5000, 5000)


# Simulate PCR duplicates -- all reads generated from a limited space of original sequences
# nseq_orig sequences amplified up to nseq_reads.

nseq_orig = 1e6
nseq_reads = 10e6


# Generate reads at random positions. Number of replicates given by the binomial distribution.


# X

