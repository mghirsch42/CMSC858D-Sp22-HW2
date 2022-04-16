#ifndef INCLUDED_SUFFIXARRAY
#define INCLUDED_SUFFIXARRAY

#include <string>
#include <sdsl/construct_sa.hpp>

using namespace std;
using namespace sdsl;

class SuffixArray {
    public:
        // public members
        string *reference;
        int_vector<> sa;
        int length;

        SuffixArray();
        ~SuffixArray();

        bool buildsa(string& reference_fname, string& output_fname);
        bool querysa(string& index_fname, string& query_fname, string& query_mode, string& output_fname);        

        map<string,  string> read_fasta(string& reference_fname);
        void get_suffix(int index);
        void print_sa();
};

#endif