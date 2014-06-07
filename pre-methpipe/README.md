pre-methpipe
============

This directory is my own staging ground for improvements to specific
programs in methpipe. I hope members of my lab (or anyone else) will
improve these, or at least test them, and help me get better
functionality into methpipe.

roimethstat2
------------

In theory this program will do the identical computation that the
roimethstat program from methpipe does... The difference is that when
you have only few regions of interest, rather than loading the entire
file of single-site methylation levels, the ones that are required are
found by doing a binary search on disk for the proper
location. Hopefully this will make it into methpipe at some point.

methcounts2
-----------

This program is a totally rewritten version of methcounts, which
currently always prints information about all cytosines, with stranded
information, and also data on whether or not it seems like a
particular CpG site has been lost by deamination which otherwise would
appear like an unmethylated site. This program should at some point
get worked into methpipe.


allelicmeth2
------------

This one has a strange history. Our first method for identifying
allele-specific methylation computed an allelic score at each CpG site
(corresponding to a pair of CpGs) and then identified peaks in that
score. It was replaced when Fang developed the amrfinder program. But
this approach is still good for visualization, and so I'm trying to
revive it by re-writing a cleaner version.

bsrate2
-------

The only difference between this code and the regular bsrate is that I
want to add some statistics showing the distribution across
molecules. For example, it would be nice to know if the unconverted
sites are distributed uniformly at random among reads, or if there are
some really bad reads that have an unusual amount of unconverted
cytosines.
