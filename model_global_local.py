# Toy model for determining the relationship between local and global duplicates

import numpy

# Area parameters
x_lim_tot = (1100, 33000)
y_lim_tot = (1000, 50000)
n_tiles = 112

# Max distance along the axis (so, the window size is about twice the *_lim_dist)
x_lim_dist = 2500
y_lim_dist = 2500


# Generate nseq_reads by sampling randomly from the space of [0,nseq_orig].
def generate_reads(nseq_orig, nseq_reads):
    # Simulate that we include first nseq_orig reads, representing the original library
    # that undergoes PCR. 
    #reads = numpy.arange(min(nseq_reads, nseq_orig)) 
    # Then add on duplicates if we require more reads
    reads = numpy.random.randint(0, nseq_orig, nseq_reads)
    xs = numpy.random.randint(*x_lim_tot, (nseq_reads))
    ys = numpy.random.randint(*y_lim_tot, (nseq_reads))
    tiles = numpy.random.randint(0, n_tiles, (nseq_reads))
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
        dup_coords = numpy.stack(
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


total_reads = 10000000

# Code for analysis of window size
reads, xs, ys, tiles = generate_reads(10, total_reads) 
coords = numpy.stack((tiles, ys, xs), axis=1)

def eval_in_range(ctyx):
    tile, y, x = ctyx
    return numpy.sum(
                (coords[:,0] == tile) &
                (numpy.fabs(coords[:,1] - y) < y_lim_dist) &
                (numpy.fabs(coords[:,2] - x) < x_lim_dist)
                )

# Code for main analysis
def analyse(library_size):
    reads, xs, ys, tiles = generate_reads(library_size, total_reads)
    return get_duplicate_counts(reads, xs, ys, tiles, x_lim_dist, y_lim_dist)

# Analysis with sub-libraries comprising a 
def analyse_with_sublibraries(param):
    sub_library_sizes, sub_fractions_of_reads = param
    remaining = 1.0 - sum(sub_fractions_of_reads)
    if remaining > 0.00001: # Generate moderately random reads for remaining
        sub_library_sizes.append(total_reads * remaining)
        sub_fractions_of_reads.append(remaining)

    reads, xs, ys, tiles = [numpy.empty(total_reads) for _ in range(4)]
    start_read = 0
    for sub_library_size, sub_fraction_of_reads in\
                    zip(sub_library_sizes, sub_fractions_of_reads):
        num_reads = int(total_reads*sub_fraction_of_reads)
        end_read = start_read + num_reads
        print("sub_library_size", sub_library_size)
        if sub_library_size < 1:
            sub_library_size = 1
        reads[start_read:end_read],\
           xs[start_read:end_read],\
           ys[start_read:end_read],\
        tiles[start_read:end_read] = generate_reads(sub_library_size, num_reads)
        start_read += num_reads

    return get_duplicate_counts(reads, xs, ys, tiles, x_lim_dist, y_lim_dist)


