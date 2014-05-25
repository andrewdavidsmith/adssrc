adssrc
======

Some of the code I use in my own data analysis. Some of this might make it into smithlabcode repositories eventually.

collapsebed
-----------

This program takes a set of genomic intervals, and identifies each maximal
contiguous sub-iterval that overlaps a fixed number of the input intervals.
Then the set of sub-intervals overlapping more than some specified number
of the original intervals are reported.

majormethstate
--------------

The proagram looks at each "epiread" in an epiread formatted file (of the
kind used in amrfinder) and determines if the number of methylated states in
the read is lower than 0.5. If so, then the states are inverted: all C changed
to T, and all T changed to C. Then the epiread is printed to a file or stdout.
If you don't know why this might be a useful thing to do, then you probably
don't need to do it.

roimethstat2
------------

In theory this program will do the identical computation that the roimethstat
program from methpipe does... The difference is that when you have only few
regions of interest, rather than loading the entire file of single-site methylation
levels, the ones that are required are found by doing a binary search on disk
for the proper location. Hopefully this will make it into methpipe at some point.

smoothmeth
----------

In methpipe the hmr program uses a 2-state HMM to identify hypomethylated regions (HMRs)
using a beta-binomial model for emissions. This program takes the same approach for
smoothing methylation levels. Each of the 2 states has an average methylation level,
and the posterior probabilities for state occupancy at any CpG site are extracted from
the HMM. So if p is the posterior for being in the high methylation state, and h is the
methylation level for the "high" state and l is the methylation level for the "low"
state, then the smoothed level is ph + (1 - p)*l. This program might eventually replace
how we build tracks for single-CpG methylation levels.

tsscpgplot
----------

I use this program to make meta-gene plots of methylation levels around transcription
start sites (TSS) but also around other fixed landmarks in the genome. The program
takes a BED format file, and processes it in a stranded way, but currently it only
works if the BED interval is a single site, so the interval has size 1.

methcounts2
-----------

This program is a totally rewritten version of methcounts, which currently always
prints information about all cytosines, with stranded information, and also data
on whether or not it seems like a particular CpG site has been lost by deamination
which otherwise would appear like an unmethylated site. This program should at
some point get worked into methpipe.

