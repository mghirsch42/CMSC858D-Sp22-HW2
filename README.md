# CMSC 858D Spring 2022 Homework 2

See full instructions here: https://rob-p.github.io/CMSC858D_S22/assignments/02_suffix_arrays

This assignment is written in C++.

## How to run

Compile with:
```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib runner.cpp -o runner suffix_array.cpp -lsdsl -ldivsufsort -ldivsufsort64
```

## Code

## Resources

This program uses sdsl-lite (https://github.com/simongog/sdsl-lite) to construct the suffix array.

To read the FASTA file, I use code modified from https://rosettacode.org/wiki/FASTA_format#C.2B.2B. 