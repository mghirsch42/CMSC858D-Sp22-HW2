#include "suffix_array.hpp"

using namespace std;

void print_help() {
    cout << "HELP:" << endl;
    cout << "This program four required arguments in the following order:" << endl;
    cout << "index : the path to the binary file containing the suffix array." << endl;
    cout << "queries : the pth to a FASTA file containing a set of query sequences." << endl;
    cout << "query_mode : the query mode to be used, must be either 'naive' or 'simpaccel'." << endl;
    cout << "output: the program will write a single binary output file with this name that contains the results of the queries." << endl; 
};

int main(int argc, char *argv[]) {
    if (argc < 5) {
        print_help();
        exit(0);
    }

    string index = argv[1];
    string queries = argv[2];
    string query_mode = argv[3];
    string output = argv[4];
    string time_fname = string();

    if (argc >= 7) {
        if (strcmp(argv[5], "--time_fname") == 0) {
            time_fname = argv[6];
        }
    }

    SuffixArray *sa = new SuffixArray();
    // sa->load(index);
    sa->querysa(index, queries, query_mode, output, time_fname);
    delete(sa);
};
