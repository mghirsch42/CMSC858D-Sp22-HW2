#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <chrono>
#include <sdsl/construct_sa.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/string.hpp>
#include "suffix_array.hpp"


using namespace std;
using namespace std::chrono;
using namespace sdsl;

SuffixArray::SuffixArray() {
    length = 0;
    preftab_k = 0;
};

SuffixArray::~SuffixArray() {

};

void SuffixArray::buildsa(string& reference_fname, string& output_fname, int preftab_k) {
    string query_mode = "naive";
    string time_fname = string();
    buildsa(reference_fname, output_fname, query_mode, time_fname, preftab_k);
}

void SuffixArray::buildsa(string& reference_fname, string& output_fname, string& query_mode, string& time_fname, int preftab_k) {
    cout << "Building the suffix array." << endl;

    // check that mode is valid
    if (query_mode.compare("naive") != 0 && query_mode.compare("simpaccel") != 0) {
        cout << "Query mode must be 'naive' or 'simpaccel'." << endl;
        return;
    }
    bool query_accel = (query_mode.compare("naive") == 0) ? false : true;
    if (query_accel) {
        cout << "Using simple acceleration for querying to build the prefix table." << endl;
    }

    high_resolution_clock::time_point start_time;
    high_resolution_clock::time_point end_time;
    duration<double, std::milli> diff_time;

    // Read the reference
    cout << "Reading reference from " << reference_fname << "." << endl;
    map<string, string> fasta = this->read_fasta(reference_fname);
    map<string, string>::iterator fasta_it = fasta.begin();
    reference = fasta_it->second;
    clean_reference();
    length = reference.length();

    cout << "Starting to build the suffix array." << endl;

    start_time = high_resolution_clock::now();

    // Create the suffix array
    sa = int_vector<>(length, 0, length);
    sdsl::algorithm::calculate_sa(reinterpret_cast<const unsigned char*>(reference.c_str()), reference.length(), sa);

    // Create the prefix table
    this->preftab_k = preftab_k;
    if (preftab_k > 0) {
        build_preftab(preftab_k, query_accel);
    }

    end_time = high_resolution_clock::now();
    diff_time = (end_time - start_time);

    cout << "Finished building the suffix array. Time: " << diff_time.count() << " ms" << endl;

    cout << "Saving suffix array to " << output_fname << "." << endl;
    save(output_fname);
    cout << "Suffix array saved." << endl;

    if (!time_fname.empty()) {
        cout << "Saving time information to " << time_fname << "." << endl;
        ofstream outfile;
        outfile.open(time_fname, std::ios_base::app);
        outfile << length << "," << preftab_k << "," << query_mode << "," << diff_time.count() << endl;
        cout << "Time information saved." << endl;
    }
};

void SuffixArray::querysa(string& index_fname, string& query_fname, string& query_mode, string& output_fname, string& time_fname) {
    cout << "Querying" << endl;

    // check that query mode is valid
    if (query_mode.compare("naive") != 0 && query_mode.compare("simpaccel") != 0) {
        cout << "Query mode must be 'naive' or 'simpaccel'." << endl;
        return;
    }

    bool query_accel = (query_mode.compare("naive") == 0) ? false : true;
    if (query_accel) {
        cout << "Using simple acceleration for querying." << endl;
    }
    // Load index
    cout << "Loading index from " << index_fname << "." << endl;
    load(index_fname);

    // Read queries
    cout << "Reading queries from " << query_fname << "." << endl;
    map<string, string> queries = read_fasta(query_fname);

    cout << "Running querying." << endl;
    high_resolution_clock::time_point start_time;
    high_resolution_clock::time_point end_time;
    duration<double, std::milli> diff_time;

    start_time = high_resolution_clock::now();

    string result = "";
    int q_len;
    for (map<string, string>::iterator it = queries.begin(); it!=queries.end(); it++) {
        result += it->first;

        string query = it->second;
        q_len = (int)query.length();

        int left_index = 0;
        int right_index = length;

        // If the prefix is in the prefix table, start there
        if (preftab_k > 0 && preftab.count(query.substr(0, preftab_k)) != 0) {
            left_index = get<0>(preftab[query.substr(0, preftab_k)]);
            right_index = get<1>(preftab[query.substr(0, preftab_k)]);
        }

        int count;
        if (left_index == -1 || right_index == -1) {
            count = 0;
        }
        else {
            if (right_index != left_index) {
                left_index = low_index(query, left_index, right_index, query_accel);
                if (left_index != -1) { // There's no reason to find the right index if the left index doesn't exist
                    right_index = high_index(query, left_index, right_index);
                }
            }
            else {
                if (get_suffix(left_index, query.length()) != query) {
                    left_index = right_index = -1;
                }
            }
        }

        if (left_index < 0 || right_index < 0) {
            count = 0;
        }
        else {
            count = right_index - left_index + 1;
        }

        if (count == 0) {
            result += "\t0\n";
        }
        else {
            result += "\t" + to_string(count);
            for (int i=0; i<count; i++) {
                int pos = sa[left_index + i];
                result += "\t" + to_string(pos);
            }
            result += "\n";
        }   
    }

    end_time = high_resolution_clock::now();
    diff_time = (end_time - start_time);

    cout << "Finished querying. Time: " << diff_time.count() << " ms for " << queries.size() << " queries." << endl;

    cout << "Saving query results to " << output_fname << "." << endl;
    ofstream outfile;
    outfile.open(output_fname);
    outfile << result;
    outfile.close();
    cout << "Finished saving query results." << endl;

    if (!time_fname.empty()) {
        cout << "Saving time information to " << time_fname << "." << endl;
        outfile.open(time_fname, std::ios_base::app);
        outfile << length << "," << preftab_k << "," << query_mode << "," << q_len << "," << queries.size() << "," << diff_time.count() << endl;
        cout << "Time information saved." << endl;
    }
};

map<string, string> SuffixArray::read_fasta(string& fname){
    // Code modified from https://rosettacode.org/wiki/FASTA_format#C.2B.2B
    ifstream input(fname);
    if (!input.good()) {
        cerr << "Error opening fasta file" << endl;
    }
    map<string, string> seqs;
    string line, name, content;
    while(getline(input, line).good()) {
        // cout << line << endl;
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

void SuffixArray::build_preftab(int k, bool query_accel) {
    vector<string> kmers;
    all_kmers(k, kmers);
    preftab = {};
    for (string s: kmers) {
        if (preftab.count(s) == 0) {
            int low = low_index(s, 0, -1, query_accel);
            int high = high_index(s, 0, -1, query_accel);
            preftab[s] = {low, high};
        }
    }
};

void SuffixArray::all_kmers(int k, vector<string>& kmers) {
    // cout << "all kmers" << endl;
    for (int i=0; i<length-k; i++) {
        kmers.push_back(reference.substr(i, k));
    }
};

void SuffixArray::clean_reference() {
    for(int i=0; i<length; i++) {
        char c = toupper((reference)[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            int r = rand() % 4;
            switch(r) {
                case 0:
                    c = 'A';
                    break;
                case 1:
                    c = 'C';
                    break;
                case 2:
                    c = 'G';
                    break;
                case 3:
                    c = 'T';
                    break;
                default:
                    c = 'A';
            }
        }
        (reference)[i] = c;
    }
};

string SuffixArray::get_suffix(int index, int k=-1) {
    if (index >= length) {
        return "";
    }
    string suffix = "";
    if (k<=0) {
        for (int i=sa[index]; i<length; i++) {
            suffix += (reference)[i];
        }
    }
    else {
        for (int i=sa[index]; i<sa[index]+k; i++) {
            suffix += (reference)[i];
        }
    }
    return suffix;
};

int SuffixArray::low_index(string& query, int left_index, int right_index, bool query_accel) {    
    int k = query.length();

    if (right_index == -1) {
        right_index = length;
    }

    // Variables used for simple acceleration. 
    // If using naive algorithm, these stay 0
    int l_lcp = 0;
    int r_lcp = 0;
    int c_lcp = 0;
    int m_lcp = 0;

    if (query_accel) {
        string l_suffix = get_suffix(left_index, -1);
        string r_suffix = get_suffix(right_index, -1);
        l_lcp = lcp(query, l_suffix);
        r_lcp = lcp(query, r_suffix);
    }
    if (get_suffix(left_index, -1).compare(c_lcp, k, query, c_lcp, k) == 0) {
        return left_index;
    }

    while (right_index > left_index) {
        int mid_index = (int) floor(left_index + (right_index - left_index) / 2);
        // cout << "Low: " << left_index << " High: " << right_index << " Mid: " << mid_index << endl;
        
        // Get current suffix to compare to
        string curr_suffix = get_suffix(mid_index, -1);

        if (query_accel) {
            m_lcp = lcp(query, curr_suffix);
            c_lcp = min(l_lcp, r_lcp);
            k = query.length() - c_lcp;
        }
        

        // Get comparison location and length; taking the min of the lengths of the strings
        // and the lcp values in case we have a suffix that is shorter than our query we don't
        // go out of range 
        int comp_loc = min(c_lcp, min((int)curr_suffix.length(), (int)query.length()));
        int comp_len = min(k, min((int)curr_suffix.length(), (int)query.length()));
        // cout << "comp loc " << comp_loc << " comp len " << comp_len << " query " << query.length() << " suffix " << curr_suffix.length() << endl;  
        // The query is equal to the current suffix
        if (curr_suffix.compare(comp_loc, comp_len, query, comp_loc, comp_len) == 0) {
            // cout << "if 1" << endl;
            if (mid_index == left_index + 1) {
                return mid_index;
            }
            right_index = mid_index;
            if (query_accel) {
                r_lcp = m_lcp;
            }
        } 
        // The query is larger than the current suffix
        else if (curr_suffix.compare(comp_loc, comp_len, query, comp_loc, comp_len) < 0) {
            // cout << "if 2" << endl;
            if (mid_index == right_index - 1) {
                // cout << "length get suffix right index " << get_suffix(right_index, -1).length() << endl;
                if (get_suffix(right_index, -1).compare(comp_loc, comp_len, query, comp_loc, comp_len) == 0) {
                    // cout << "return " << right_index << endl;
                    return right_index;
                }
                else {
                    // cout << "return -1" << endl;
                    return -1;
                }
            }
            // cout << "update" << endl;
            left_index = mid_index;
            if (query_accel) {
                l_lcp = m_lcp;
            }
        }
        // Query is smaller than the current suffix
        else if (curr_suffix.compare(comp_loc, comp_len, query, comp_loc, comp_len) > 0) {
            // cout << "if 3" << endl;
            if (mid_index == left_index + 1) {
                return -1;
            }
            right_index = mid_index;
            if (query_accel) {
                r_lcp = m_lcp;
            }
        }
    }
    return -1;
};

int SuffixArray::high_index(string& query, int left_index, int right_index, bool query_accel) {
    int k = query.length();

    if (right_index == -1) {
        right_index = length;
    }

    // Variables used for simple acceleration. 
    // If using naive algorithm, these stay 0
    int l_lcp = 0;
    int r_lcp = 0;
    int c_lcp = 0;
    int m_lcp = 0;

    if (query_accel) {
        string l_suffix = get_suffix(left_index, k);
        string r_suffix = get_suffix(right_index, k);
        l_lcp = lcp(query, l_suffix);
        r_lcp = lcp(query, r_suffix);
    }

    if (get_suffix(right_index, -1).compare(0, k, query, 0, k) == 0) {
        return right_index;
    }

    while (right_index > left_index) {
        int mid_index = (int) floor(left_index + (right_index - left_index) / 2);
        // cout << "Low: " << left_index << " High: " << right_index << " Mid: " << mid_index << endl;

        // Get current suffix to compare to
        string curr_suffix = get_suffix(mid_index, -1);

        if (query_accel) {
            m_lcp = lcp(query, curr_suffix);
            c_lcp = min(l_lcp, r_lcp);
            k = query.length() - c_lcp;
        }

        // Get comparison location and length; taking the min of the lengths of the strings
        // and the lcp values in case we have a suffix that is shorter than our query we don't
        // go out of range 
        int comp_loc = min(c_lcp, min((int)curr_suffix.length(), (int)query.length()));
        int comp_len = min(k, min((int)curr_suffix.length(), (int)query.length()));

        // The query is equal to the current suffix
        if (curr_suffix.compare(comp_loc, comp_len, query, comp_loc, comp_len) == 0) {
            if (mid_index == right_index - 1) {
                return mid_index;
            }
            left_index = mid_index;
            if (query_accel) {
                l_lcp = m_lcp;
            }
        }
        // The query is larger than the current suffix
        else if (curr_suffix.compare(comp_loc, comp_len, query, comp_loc, comp_len) < 0) { // query is larger
            if (mid_index == right_index - 1) {
                return -1;
            }
            left_index = mid_index;
            if (query_accel) {
                r_lcp = m_lcp;
            }
        }
        // Query is smaller than the current suffix
        else if (curr_suffix.compare(comp_loc, comp_len, query, comp_loc, comp_len) > 0) { // query is smaller
            if (mid_index == left_index + 1) {
                if(get_suffix(left_index, -1).compare(comp_loc, comp_len, query, comp_loc, comp_len) == 0) {
                    return left_index;
                }
                else {
                    return -1;
                }
            }            
            right_index = mid_index;
            if (query_accel) {
                r_lcp = m_lcp;
            }
        }
    }
    return -1;
};

int SuffixArray::lcp(string& a, string& b) {
    int lcp = 0;
    for (int i=0; i<min(a.length(), b.length()); i++) {
        if (a[i] == b[i]) {
            lcp++;
        }
        else {
            break;
        }
    }
    return lcp;
};

void SuffixArray::print_sa() {
    if(length == 0) {
        cout << "Suffix array has not been constructed." << endl;
    }
    else {
        for (int i=0; i<length; i++) {
            cout << i << ": " << sa[i] << ": " << get_suffix(i, -1) << endl;
        } 
        // cout << endl;
    }
};

void SuffixArray::print_preftab() {
    if(preftab_k == 0) {
        cout << "Prefix table has not been constructed." << endl;
    }
    else {
        for(unordered_map<string, tuple<int, int>>::iterator it=preftab.begin(); it!=preftab.end(); it++) {
            cout << it->first << ": (" << get<0>(it->second) << ", " << get<1>(it->second) << ")" << endl;
        }
    }
};

void SuffixArray::preftab_to_string() {
    preftab_string = "";
    for(unordered_map<string, tuple<int, int>>::iterator it=preftab.begin(); it!=preftab.end(); it++) {
        preftab_string += it->first + " " + to_string(get<0>(it->second)) + " " + to_string(get<1>(it->second)) + "\n";
    }
}

void SuffixArray::preftab_from_string() {
    stringstream ss (preftab_string);
    string line;
    preftab = {};

    while (getline(ss, line)) {
        // cout << line << endl;
        int pos_start = 0;
        int pos_end = line.find(" ", pos_start);
        string kmer = line.substr(pos_start, pos_end-pos_start);
        pos_start = pos_end + 1;
        pos_end = line.find(" ", pos_start);
        int start_index = stoi(line.substr(pos_start, pos_end-pos_start));
        pos_start = pos_end + 1;
        int end_index = stoi(line.substr(pos_start, line.length() - pos_start));
        preftab[kmer] = {start_index, end_index};
    }
}

void SuffixArray::save(string& output_fname) {

    preftab_to_string();

    ofstream outfile;
    outfile.open(output_fname, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(outfile);
    oarchive(reference, preftab_k, preftab_string);
    sa.serialize(outfile);
    outfile.close();
};

void SuffixArray::load(string& input_fname) {

    ifstream infile;
    infile.open(input_fname, std::ios::binary);
    cereal::BinaryInputArchive iarchive(infile);
    iarchive(reference, preftab_k, preftab_string);
    sa.load(infile);
    infile.close();
    length = reference.length();
    preftab_from_string();
};

void SuffixArray::save_readable(string& output_fname) {
    ofstream outfile;
    outfile.open(output_fname);

    // Save prefix table k
    outfile << "Prefix table k: " << preftab_k << endl;

    // If we have a prefix table, save that
    if (preftab_k > 0) {
        outfile << "Prefix table size: " << preftab.size() << endl;
        for (unordered_map<string, tuple<int, int>>::iterator it=preftab.begin(); it!=preftab.end(); it++) {
            outfile << it->first << ": ";
            outfile << "(" << get<0>(it->second);
            outfile << ", " << get<1>(it->second) << ")" << endl; 
        }
    }

    cout << endl;

    // Save the length
    outfile << "Reference length: " <<  length << endl;

    // Save the reference
    outfile << "Reference: " << endl << reference << endl;

    // Save the suffix array
    for (int i=0; i<length; i++) {
        outfile << i << ": " << sa[i] << ": " << get_suffix(i) << endl;
    } 
    outfile.close();
};

template <class Archive>
void SuffixArray::serialize(Archive &archive) {
    archive(reference, preftab_k, preftab_string);
};