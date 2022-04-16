
#include "suffix_array.hpp"

using namespace std;

int main() {
    string reference_fname = "data/test_ref.fa";
    string query_fname = "data/query.fa";
    string output_fname = "data/output.dat";
    SuffixArray* sa = new SuffixArray();
    sa->buildsa(reference_fname, output_fname);
    sa->print_sa();
    sa->get_suffix(2);
}