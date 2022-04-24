# CMSC 858D Spring 2022 HW 1

MG Hirsch

See instructions here: https://rob-p.github.io/CMSC858D_S22/assignments/02_suffix_arrays

Github repository: https://github.com/mghirsch42/CMSC858D-Sp22-HW2


## Part 1: Suffix Array Construction

### Scaling

### Most challenging part

Setting up the prefix table was the hardest part for me. At first, I tried to mainly copy the query algorithm from the slides, but I found that some modifications were needed to get the lowest and highest index. I had to back away from the code and write out the different cases that could arise and would need different handling: when the query was equal, larger, or smaller to the current suffix, and for each of these if the query was one above the low index, one or two below the high index, or neither. I had to repeat this process two or three times after coded it when I was trying to fix bugs.

## Part 2: Suffix Array Querying

### Scaling


### Most challenging part

I basically used the same code as when setting up the prefix table, and as mentioned above, I had some challenges with that. Aside from that, I think the hardest part of this was just trying different test reference texts and queries and debugging the different cases of results that could come out of the high and low index methods.


## Code

### Suffix array construction

I read in the FASTA file using code modified from https://rosettacode.org/wiki/FASTA_format#C.2B.2B. Then I clean the reference, turning it into all capital letters and replacing any characters that aren't A, C, G, or T into one of them uniformly at random. I used SDSL-Lite to construct the suffix array. The suffix array is stored in an integer vector.

### Prefix table construction

To build the prefix table, first, I get all k-mers from the reference string. I do not include k-mers that are not included in the reference, as those would need to search the entire suffix array so there is no reason to store their range. Then for each unique k-mer, I calculate the range of the suffix array by doing two slightly modified binary searches. These searches will return the indices of the suffix array that contain the lowest and highest matching suffixes up to k. I store this in a map where the keys are the k-mers and the values are tuples containing the low and high values of the suffix array range.

### Querying

To query the suffix array, I get an initial range from the prefix table by looking up the length k prefix of the query. Then I use the same methods as in the prefix table construction to find the high and low indices of matching suffixes, starting at the indices from the prefix table if any were found. If using simple acceleration, these methods will keep track of the longest common prefix during the binary search and only do character comparisons on characters after that prefix. If no matches are found, the method reports a count of 0. Otherwise, it reports the count as the high index - low index + 1. It also returns the values of the suffix array at those indices. These results are written to a text file.

### Saving and loading

I use the serialization library Cereal. Because I stored my prefix table as an unordered map from strings to tuples and Cereal does not support that (or at least I couldn't get it to work), I first turn the map into a string and then serialize the string. Then when reading, I read the serialized string into the string variable and use string methods to put the values into the map. I also serialize the reference and the prefix table k value. These did not need extra work to serialize. I do not save the length because I can calucate that from the reference after reading it in.

## Resources

To construct the suffix array, I use SDSL-lite (https://github.com/simongog/sdsl-lite).

To read the FASTA file, I use code modified from https://rosettacode.org/wiki/FASTA_format#C.2B.2B. 

To serialize the data, I use Cereal (https://uscilab.github.io/cereal/index.html).