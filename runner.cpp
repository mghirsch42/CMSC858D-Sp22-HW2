#include <queue>
#include "suffix_array.hpp"


using namespace std;

int main() {
    string reference_fname = "data/ecoli_small.fa";
    string query_fname = "data/test_query.fa";
    string output_fname = "data/output.dat";
    SuffixArray* sa = new SuffixArray();
    sa->buildsa(reference_fname, output_fname, 3);
    sa->save(output_fname);
    // sa->preftab_from_string();
    SuffixArray* sa2 = new SuffixArray();
    sa2->load(output_fname);
    // sa2->print_sa();
    sa2->print_preftab();
    // sa->save_readable(output_fname);
    // string query_mode = "naive";
    string query_output = "data/query_out.txt";
    string query_mode = "simpaccel";
    sa2->querysa(reference_fname, query_fname, query_mode, query_output);
}