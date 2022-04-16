#ifndef INCLUDED_PREFIXTABLE
#define INCLUDED_PREFIXTABLE

#include <string>

using namespace std;

class PrefixTable {
    public:
        // public members
        PrefixTable();
        ~PrefixTable();

        bool buildpt(string& reference_fname, string& output_fname);
}

#endif