#ifndef INCLUDED_SUFFIXARRAY
#define INCLUDED_SUFFIXARRAY

#include <string>
#include <sdsl/construct_sa.hpp>

using namespace std;
using namespace sdsl;

class SuffixArray {
    public:
        string reference;
        int_vector<> sa;
        int length;
        int preftab_k;
        unordered_map<string, tuple<int, int>> preftab;
        string preftab_string;

        SuffixArray();
        ~SuffixArray();

        void buildsa(string& reference_fname, string& output_fname, int preftab_k=0);
        void buildsa(string& reference_fname, string& output_fname, string& query_mode, string& time_fname, int preftab_k=0);
        void querysa(string& index_fname, string& query_fname, string& query_mode, string& output_fname, string& time_fname);

        // Helper functions 
        map<string,  string> read_fasta(string& reference_fname);   
        void build_preftab(int k, bool query_accel);
        void all_kmers(int k, vector<string>& kmers);
        void clean_reference();
        string get_suffix(int index, int k);
        int low_index(string& query, int left_index=0, int right_index=-1, bool accel=false);
        int high_index(string& query, int left_index=0, int right_index=-1, bool accel=false);
        int lcp(string& a, string& b);

        void print_sa();
        void print_preftab();

        void preftab_to_string();
        void preftab_from_string();

        void save(string& output_fname);
        void load(string& input_fname);
        void save_readable(string& output_fname);

        template <typename S>
        void serialize(S& s);
};   

#endif