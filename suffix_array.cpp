#include <string>
#include <iostream>
#include <fstream>
#include <sdsl/construct_sa.hpp>
#include "suffix_array.hpp"


using namespace std;
using namespace sdsl;

SuffixArray::SuffixArray() {
    length = 0;
};


bool SuffixArray::buildsa(string& reference_fname, string& output_fname) {
    map<string, string> fasta = this->read_fasta(reference_fname);
    map<string, string>::iterator fasta_it = fasta.begin();
    reference = &(fasta_it->second);
    length = reference->length();
    sdsl::algorithm::calculate_sa(reinterpret_cast<const unsigned char*>((*reference).c_str()) , reference->length(), sa);
    // Need to add saving to a file but when I try to open a file I get a seg fault. 
    // if I comment out the calculate_sa call, I can open and write to a file fine.
    return true;
};

// Code modified from https://rosettacode.org/wiki/FASTA_format#C.2B.2B
map<string, string> SuffixArray::read_fasta(string& fname){
    ifstream input(fname);
    if (!input.good()) {
        cerr << "Error opening query file" << endl;
    }
    map<string, string> seqs;
    string line, name, content;
    while(getline(input, line).good()) {
        if(line.empty() || line[0] == '>') {
            if(!name.empty()) {
                if(!content.empty()) {
                    seqs[name] = content;
                    // cout << name << " " << content << endl;
                }
                name.clear();
            }
            if(!line.empty()) {
                name = line.substr(1);
            }
            content.clear();
        }
        else if(!name.empty()) {
            // can add error checking for invalid characters
            content += line;
        }
    }
    if(!name.empty()) {
        if(!content.empty()) {
            seqs[name] = content;
            // cout << name << " " << content << endl;
        }
    }
    input.close();
    return seqs;
};

void SuffixArray::print_sa() {
    if(length == 0) {
        cout << "Suffix array has not been constructed." << endl;
    }
    else {
        for (int i=0; i<length; i++) {
            cout << sa[i];
        }
        cout << endl;
    }
};

void SuffixArray::get_suffix(int index) {
    if (index >= length) {
        cerr << "Index " << index << " out of range for suffix array length " << length << endl;
    }
    string suffix = "";
    for (int i=sa[index]; i<length; i++) {
        suffix += (*reference)[i];
    }
    cout << suffix << endl;
    
};