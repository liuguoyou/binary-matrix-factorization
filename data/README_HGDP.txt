Given an input file with one gene per row, and one individual per column
With m rows and n columns (m genes, n individuals), the parser produces 4 
files:

stat.txt ...... four probability vectors, corresopnding to the observed frequencies of A,T,G and C in the pairs of nucleotides of each rows
mode.pbm ...... A pseudo-binary image in PBM format where each row encodes the most frequent of A,T,G and C using four binary pixels as follows
                  A --> 0 0 0 1
                  T --> 0 0 1 0
                  G --> 0 1 0 0
                  C --> 1 0 0 0
mask.pbm ...... A pseudo-binary image of size (2m x n) where every pair of rows corresponds to one row from the original file. 
                The matrix contains two ones, one on top of each other, if the corresponding row in the input file is a known datum
                and 0 if it is missing ("--" in the input file). Thus, this file indicates the locations where data is known with 1s and
                0s where it is not.
dist.pbm ...... Given the computed modes which are stored in mode.pbm, this file encodes the "genetic distance" between the nucleotides
                of a row and its mode. Again, this matrix is encoded as a pseudo-image of 2m rows and n columns, where each pair
                of rows encodes a single row from the original file as follows. Given the j-th pair of the input column,
                the corresponding pair of elements in the corresponding two consecutive rows is given as follows (above, below):
                 (0,0) if both nucleotides equal the mode (ex. AA and the mode is A)
                 (0,1) if one of the nucleotides differs from the mode (ex. AG or GA and the mode is A)
                 (1,0) if both nucleotides differ from the mode (ex. GG and the mode is A)
                 (0,0) is also assigned when the curresponding datum is missing.

