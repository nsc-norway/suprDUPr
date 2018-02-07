# Toy model for determining the relationship between local and global duplicates

import numpy
import multiprocessing

print("-- Duplicate count toy model --")

# Area parameters
x_lim_tot = (1100, 33000)
y_lim_tot = (1000, 50000)
n_tiles = 112

print("")
print("Coordinate ranges:")
print("  x: ", x_lim_tot)
print("  y: ", y_lim_tot)
print("Number of tiles:", n_tiles)
print("")

# Max distance along the axis (so, the window size is about twice the *_lim_dist)
x_lim_dist = 2500
y_lim_dist = 2500
print("Distance thresholds for local duplicates:")
print("  x: ", x_lim_dist)
print("  y: ", x_lim_dist)
print("")

# Generate nseq_reads by sampling randomly from the space of [0,nseq_orig].
def generate_reads(nseq_orig, nseq_reads):
    # Simulate that we include first nseq_orig reads, representing the original library
    # that undergoes PCR
    #reads = numpy.arange(min(nseq_reads, nseq_orig)) 
    # Then add on duplicates if we require more reads
    reads = numpy.random.random_integers(1, nseq_orig, nseq_reads)
    xs = numpy.random.random_integers(*x_lim_tot, (nseq_reads))
    ys = numpy.random.random_integers(*x_lim_tot, (nseq_reads))
    tiles = numpy.random.random_integers(1, n_tiles, (nseq_reads))
    return reads, xs, ys, tiles


# ** Function to compute the duplication ratios, global and local **
def get_duplicate_counts(reads, xs, ys, tiles, x_lim_dist, y_lim_dist):

    # ** Global duplication ratio: count number of duplicates for each unique sequence, regardless
    # of position**
    
    # First get unique sequences (numbers) and the count of each. Only the count is needed for the
    # global duplication ratio.
    
    # We also request return_inverse, which returns an array which is as long as reads, which 
    # contains an index into unique_seqs (array of uniques) for each element in reads. The 
    # variable indexes is used in the next section (local duplicates).
    unique_seqs, indexes, counts = numpy.unique(reads, return_inverse=True, return_counts=True)

    num_global_duplicates = (counts - 1).sum()


    # ** Local duplication ratio **
    # Count the number of reads within a threshold distance, and in the same tile

    local_duplicate_counter = 0

    # Looping over all unique simulated sequences. We only need their index into the
    # unique_seqs array and their count.
    for i, count in ((i, count) for i, count in enumerate(counts) if count != 1):
        # The positions of all duplicates of this read in the arrays reads, xs, ys, tiles
        dup_indexes = numpy.argwhere(indexes == i)

        # Create a matrix with columns tile, y, x, rows for each of the duplicates
        # in this group of sequences
        dup_coords = numpy.concatenate(
            (tiles[dup_indexes],xs[dup_indexes],ys[dup_indexes]),
            axis=1
            )

        for j in range(1, count):
            if numpy.any(
                (dup_coords[j,0] == dup_coords[0:j,0]) &
                (numpy.fabs(dup_coords[j,1] - dup_coords[0:j,1]) < y_lim_dist) &
                (numpy.fabs(dup_coords[j,2] - dup_coords[0:j,2]) < x_lim_dist)
                ):
                local_duplicate_counter += 1

    return num_global_duplicates, local_duplicate_counter

total_reads = 10000
def analyse(library_size):
    reads, xs, ys, tiles = generate_reads(library_size, total_reads)
    return get_duplicate_counts(reads, xs, ys, tiles, x_lim_dist, y_lim_dist)

pool = multiprocessing.Pool()
n_origs = [1, 10, 20, 50] + list(range(100, 8000, 100)) + list(range(8000, 20000, 500))
global_local = pool.map(analyse, n_origs)
for o, (g, l) in zip(n_origs, global_local):
    print(o, g / total_reads, l / total_reads)


