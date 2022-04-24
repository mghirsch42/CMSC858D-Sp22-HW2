#include "suffix_array.hpp"

using namespace std;

void print_help() {
    cout << "HELP:" << endl;
    cout << "This program takes one optional parameter and two required arguments in the following order:" << endl;
    cout << "--preftab <k> : optional parameters that specifies that a prefix table should be built with prefix length k." << endl;
    cout << "reference : the path to a FASTA format file containing the reference of which you will build the suffix array." << endl;
    cout << "output: the program will write a single binary output file with this name that contains the serialized data." << endl; 
};

int main(int argc, char *argv[]) {
    if (argc < 3) {
        print_help();
        exit(0);
    }

    int k = 0;
    string reference;
    string output;

    if (strcmp(argv[1], "--preftab") == 0) {
        try {
            k = stoi(argv[2]);
        }
        catch (const exception&) {
            cout << "--preftab must be followed by an integer." << endl;
            exit(0);
        }

        if (argc >= 5) {
            reference = argv[3];
            output = argv[4];
        }
        else {
            print_help();
            exit(0);
        }
    }
    else {
        if (argc >= 3) {
            reference = argv[1];
            output = argv[2];
        }
        else {
            print_help();
            exit(0);
        }
    }

    SuffixArray *sa = new SuffixArray();
    sa->buildsa(reference, output, k);
    delete(sa);
};
