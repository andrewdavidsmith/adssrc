pst
===

This is my implementation of a probabilistic suffix tree. I'm not citing
any papers here because I can't remember who came up with this idea. Basically
a "PST" is a variable order Markov model. It is constructed by building a lexico
graphically ordered trie, like a suffix tree (uncompressed). Each node has a
probability distribution for the letter that would come after the sequence
obtained by matching (backwards) down the suffix tree to that node. All leafs
are at the same depth to start, but then a pruning happens, bottom-up, and
a leaf node is removed if its distribution is not different from the one associated
with its parent node.

I implemented the following:

* Construction/training the model based on a set of sequences supplied in a FASTA file.

* Simulation of a sequence from the trained model.

* Scoring a sequence, by log-likelihood, based on the model.
