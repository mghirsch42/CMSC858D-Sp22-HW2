# CMSC 858D Spring 2022 Homework 2

See full instructions here: https://rob-p.github.io/CMSC858D_S22/assignments/02_suffix_arrays

Github link: https://github.com/mghirsch42/CMSC858D-Sp22-HW2

This assignment is written in C++.

See a full description of the code, results, and assignment questions in the writeup.

## How to run

Compile with:
```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib -I cereal/include buildsa.cpp -o buildsa  suffix_array.cpp -lsdsl -ldivsufsort -ldivsufsort64 
```
for buildsa or 
```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib -I cereal/include querysa.cpp -o querysa  suffix_array.cpp -lsdsl -ldivsufsort -ldivsufsort64 
```
for querysa.

Run with:
```
./buildsa (optional: --preftab <k>) reference output (optional: --query_mode <q>) (optional: --time_fname <f>)
```
for buildsa or
```
./querysa index queries query_mode output (optional: --time_fname <f>)
```
for querysa.

## Resources

To construct the suffix array, I use SDSL-lite (https://github.com/simongog/sdsl-lite).

To read the FASTA file, I use code modified from https://rosettacode.org/wiki/FASTA_format#C.2B.2B. 

To serialize the data, I use Cereal (https://uscilab.github.io/cereal/index.html).