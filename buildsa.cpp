#include "suffix_array.hpp"

using namespace std;

void print_help() {
    cout << "HELP:" << endl;
    cout << "This program takes three optional parameter and two required arguments in the following order:" << endl;
    cout << "--preftab <k> : optional parameter that specifies that a prefix table should be built with prefix length k." << endl;
    cout << "reference : the path to a FASTA format file containing the reference of which you will build the suffix array." << endl;
    cout << "output: the program will write a single binary output file with this name that contains the serialized data." << endl; 
    cout << "--query_mode <q> : optional parameter that specifies which query mode to use to build the prefix table." 
            "Either 'naive' or 'simpaccel'. Only used if --preftab parameter is also specified. Default is 'naive'." << endl;
    cout << "--time_fname <f> : optional parameter that specifies the file to append timing information to." << endl;
};

int main(int argc, char *argv[]) {
    if (argc < 3) {
        print_help();
        exit(0);
    }

    int k = 0;
    string reference;
    string output;
    string query_mode = "naive";
    string time_fname = string();

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
        
        if (argc >=7) {
            if (strcmp(argv[5], "--query_mode") == 0) {
                query_mode = argv[6];
            }
            else if (strcmp(argv[5], "--time_fname") == 0) {
                time_fname = argv[6];
            }
        }
        
        if (argc >=9) {
            if (strcmp(argv[7], "--query_mode") == 0) {
                query_mode = argv[8];
            }
            else if (strcmp(argv[7], "--time_fname") == 0) {
                time_fname = argv[8];
            }
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
        else if (argc >=5) {
            cout << "argc >= 5" << endl;
            if (strcmp(argv[3], "--query_mode") == 0) {
                query_mode = argv[4];
            }
            else if (strcmp(argv[3], "--time_fname") == 0) {
                time_fname = argv[4];
            }
        }
        else if (argc >=7) {
            if (strcmp(argv[5], "--query_mode") == 0) {
                query_mode = argv[6];
            }
            else if (strcmp(argv[5], "--time_fname") == 0) {
                time_fname = argv[6];
            }
        }
        else {
            print_help();
            exit(0);
        }
    }

    cout << "Reference: " << reference << ", Output: " << output << ", Preftab k: " << k << ", Query mode: " << query_mode << ", Time fname: " << time_fname << endl;
    SuffixArray *sa = new SuffixArray();
    sa->buildsa(reference, output, query_mode, time_fname, k);
    delete(sa);
};
