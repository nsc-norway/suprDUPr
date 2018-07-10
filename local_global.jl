# Rewrite in Julia language (requires GroupSlices)
using GroupSlices

# Area parameters
x_lim_tot = (1100, 33000)
y_lim_tot = (1000, 50000)
n_tiles = 112

# Max distance along the axis (so, the window size is about twice the *_lim_dist)
x_lim_dist = 2500
y_lim_dist = 2500

# Total reads determines the precision of the simulation. Many parameters are 
# based on total_reads, and will give incorrect results if this is changed.
total_reads = Int64(300e6)

function generate_reads(nseq_orig, nseq_reads)
    # Simulated reads: There are nseq_orig different possibilities, and nseq_reads
    reads = rand(1:nseq_orig, nseq_reads)
    xs = rand(x_lim_tot[1]:x_lim_tot[2], nseq_reads)
    ys = rand(y_lim_tot[1]:y_lim_tot[2], nseq_reads)
    tiles = rand(1:n_tiles, nseq_reads)
    return hcat(reads, xs, ys, tiles)
end

function get_duplicate_counts(reads_xs_ys_tiles, x_lim_dist, y_lim_dist)
    println("Calling groupinds...")
    groupindes = groupinds(reads_xs_ys_tiles[:, 1])
    println("End for groupindes call...")

    num_global_duplicates = length(reads_xs_ys_tiles) - length(groupindes)
    num_local_duplicates = 0

    lastreport = Dates.now()
    for (i, dup_indexes) = enumerate(groupindes)
        if length(dup_indexes) > 1
            dup_coords = reads_xs_ys_tiles[dup_indexes, [4, 2, 3]]
            for j = 2:size(dup_coords, 1)
                if any(
                    (dup_coords[j,1] .== dup_coords[1:j,1]) .&
                    (abs.(dup_coords[j,2] - dup_coords[1:j,2]) .< y_lim_dist) .&
                    (abs.(dup_coords[j,3] - dup_coords[1:j,3]) .< x_lim_dist)
                    )
                    num_local_duplicates += 1
                end
            end
            now = Dates.now()
            if now - lastreport > Dates.Millisecond(120000)
                lastreport = now
                println("Completed ", i, " of ", length(groupindes), " groupsize=", length(dup_indexes),
                            " local_dups=", num_local_duplicates)
            end
        end
    end

    return num_global_duplicates, num_local_duplicates
end

function analyse_with_sublibraries(sub_library_sizes, sub_fractions_of_reads)
    remaining = 1.0 - sum(sub_fractions_of_reads)
    if remaining > 0.00001 # Generate moderately random reads for remaining
        push!(sub_library_sizes, total_reads * remaining)
        push!(sub_fractions_of_reads, remaining)
    end

    reads_xs_ys_tiles = Array{Int64,2}(total_reads,4)
    start_read = 1
    for (sub_library_size, sub_fraction_of_reads) = zip(sub_library_sizes, sub_fractions_of_reads)
        num_reads = Int64(floor(total_reads*sub_fraction_of_reads))
        end_read = Int64(floor(start_read + num_reads - 1))
        if sub_library_size < 1
            sub_library_size = 1
        end

        reads_xs_ys_tiles[start_read:end_read, :] = generate_reads(sub_library_size, num_reads)
        start_read += num_reads
    end

    return get_duplicate_counts(reads_xs_ys_tiles, x_lim_dist, y_lim_dist)
end

(g, l) = analyse_with_sublibraries([5000], [0.2])

println("[Sublibrary_Lengths] [Sublibrary_Fractions] global_dup local_dups")
println("5000 / 0.2] ", g / total_reads, " ", l / total_reads)


