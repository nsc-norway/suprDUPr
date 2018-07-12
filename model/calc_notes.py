# Notes about calculation of contribution to duplicate estimate from 
# randomly positioned pairs of duplicates.

# This calculation focuses on a small sub-population of reads taken 
# from substrings of a short sequence. The motivation is from RNA-seq
# and amplicons, which may produce some overrepresented sequences, e.g.
# very common transcripts.

# See run_model.py etc. for a model of PCR duplication.

import scipy.stats

# 350 million is typical for a HiSeq 4000
Nreads = 350e6
# Fraction of the reads which are from the sub-string (transcript, etc.) in 
# question, and its length in nucleotides.
f = 0.2
l = 5000

# Total area of flow cell, and search area, both expressed as number of pixels.
# Area approximately as for one lane of HiSeq 4000.
ntiles = 112
A = ntiles * 32000 * 48000
# Default search window in suprDUPr is +/- 2500 pixels
a = 5001*5001

print("Substring fraction f:                    {0:5.1f} %".format(f*100))
mu = Nreads*(f/l)*(a/A)
print("Area ratio:                            {0:7.1e}".format(a/A))
print("Mean number of reads in range:           {0:5.1f}".format(mu))

Pdup = (1-scipy.stats.poisson.pmf(0, mu))
print("Probability to find a duplicate nearby:  {0:5.1f} %".format(Pdup*100))
print("Contribution to duplication probability: {0:5.1f} %".format(f*Pdup*100))

