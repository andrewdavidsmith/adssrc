adssrc
======

Some of the code I use in my own data analysis.

collapsebed
===========

This program takes a set of genomic intervals, and identifies each maximal
contiguous sub-iterval that overlaps a fixed number of the input intervals.
Then the set of sub-intervals overlapping more than some specified number
of the original intervals are reported.

majormethstate
==============

The proagram looks at each "epiread" in an epiread formatted file (of the
kind used in amrfinder) and determines if the number of methylated states in
the read is lower than 0.5. If so, then the states are inverted: all C changed
to T, and all T changed to C. Then the epiread is printed to a file or stdout.
If you don't know why this might be a useful thing to do, then you probably
don't need to do it.
