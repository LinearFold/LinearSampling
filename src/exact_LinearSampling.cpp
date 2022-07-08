/*
 *LinearSampling.cpp*
 The main code for LinearSampling: Linear-Time Stochastic Sampling for RNA Secondary Structure.

 author: He Zhang
 edited by: 12/2019
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>

#include "exact_LinearSampling.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "exact_backtrace.cpp" // sampling
#define SPECIAL_HP

using namespace std;

unsigned long quickselect_partition(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper) {
    float pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);

    }
    return upper;
}

void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestC = new node_value_type[seq_length];

    bestH = new node_value_type*[seq_length];
    bestP = new node_value_type*[seq_length];
    bestM = new node_value_type*[seq_length];
    bestM2 = new node_value_type*[seq_length];
    bestMulti = new node_value_type*[seq_length];


    for (int i = 0; i < seq_length; i++){
        bestC[i] = VALUE_MIN;
        bestH[i] = new node_value_type[seq_length];
        bestP[i] = new node_value_type[seq_length];
        bestM[i] = new node_value_type[seq_length];
        bestM2[i] = new node_value_type[seq_length];
        bestMulti[i] = new node_value_type[seq_length];
    }

    for (int i = 0; i < seq_length; i++){
        for (int j = 0; j < seq_length; j++){
            bestH[i][j] = VALUE_MIN;
            bestP[i][j] = VALUE_MIN;
            bestM[i][j] = VALUE_MIN;
            bestM2[i][j] = VALUE_MIN;
            bestMulti[i][j] = VALUE_MIN;
        }
    }

    sortedP = new vector<int>[seq_length];
    nucs = new int[seq_length];
    scores.reserve(seq_length); 
    next_pair.resize(seq_length * NOTON, -1);
    prev_pair.resize(seq_length * NOTON, -1);
}

void BeamCKYParser::cleanup() {

    for (int i = 0; i< seq_length; i++){
        delete[] bestH[i];
    }
    for (int i = 0; i< seq_length; i++){
        delete[] bestP[i];
    }
    for (int i = 0; i< seq_length; i++){
        delete[] bestM[i];
    }
    for (int i = 0; i< seq_length; i++){
        delete[] bestM2[i];
    }
    for (int i = 0; i< seq_length; i++){
        delete[] bestMulti[i];
    }        

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;
    delete[] sortedP;

#ifndef non_saving
    for (int i = 0; i < 5; i++)
        delete[] samplestates[i];
    delete[] samplestates;
#endif
}

BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq) {
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            int next = -1;
            for (int j = seq_length-1; j >= 0; --j) { // going backward
                next_pair[nuci * seq_length + j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }

            int prev = -1;
            for (int j = 0; j <= seq_length-1; ++j) { // going forward
                prev_pair[nuci * seq_length + j] = prev;
                if (_allowed_pairs[nuci][nucs[j]]) prev = j;
            }
        }
    }

#ifdef SPECIAL_HP
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif

    if (!read_forest)
    {

    if(seq_length > 0) bestC[0] = 0.0;
    if(seq_length > 1) bestC[1] = 0.0;

    value_type newscore;
    // from left to right
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        node_value_type * beamstepH = bestH[j];
        node_value_type * beamstepMulti = bestMulti[j];
        node_value_type * beamstepP = bestP[j];
        node_value_type * beamstepM2 = bestM2[j];
        node_value_type * beamstepM = bestM[j];
        node_value_type & beamstepC = bestC[j];

        // beam of H
        {
            {
                int jnext = next_pair[nucj * seq_length + j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj * seq_length + jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                    int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                    if (jnext-j-1 == 4) // 6:tetra
                        tetra_hex_tri = if_tetraloops[j];
                    else if (jnext-j-1 == 6) // 8:hexa
                        tetra_hex_tri = if_hexaloops[j];
                    else if (jnext-j-1 == 3) // 5:tri
                        tetra_hex_tri = if_triloops[j];
#endif
                    newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                    Fast_LogPlusEquals(bestH[jnext][j], newscore/kT);
                }
            }

            {
                for (int i = 0; i < j; i++){
                    node_value_type & state = beamstepH[i];
                    if (state == VALUE_MIN) continue;

                    int nuci = nucs[i];
                    int jnext = next_pair[nuci * seq_length + j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-i-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[i];
                        else if (jnext-i-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[i];
                        else if (jnext-i-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[i];
#endif
                        newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][i], newscore/kT);
                    }

                    // 2. generate p(i, j)
                    {
                        Fast_LogPlusEquals(beamstepP[i], state);
                    }
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            for (int i = 0; i < j; i++){
                node_value_type & state = beamstepMulti[i];

                if (state == VALUE_MIN) continue;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci * seq_length + j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
                        newscore = - v_score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(bestMulti[jnext][i], state + newscore/kT);
                    }
                }

                // 2. generate P (i, j)
                {
                    value_type score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i], state + score_multi/kT);
                }
            }
        }

        // beam of P
        {
            for (int i = 0; i < j; i++){
                node_value_type & state = beamstepP[i];

                if (state == VALUE_MIN) continue;

                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                if (i >0 && j<seq_length-1) {
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp * seq_length + j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
                                int score_single = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                            nuci_1, nuci, nucj, nucj1);
                                 if (use_shape)
                                    score_single += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                Fast_LogPlusEquals(bestP[q][p], state + score_single/kT);
                            } else {
                                // single branch
                                int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(bestP[q][p], state + score_single/kT);
                            }
                            q = next_pair[nucp * seq_length + q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                    int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    Fast_LogPlusEquals(beamstepM[i], state + score_M1/kT);
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 ) {    
                    value_type M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = state + M1_score/kT;
                    for (int newi = 0; newi < k ; newi++){
                        node_value_type & m_state = bestM[k][newi];
                        if (m_state == VALUE_MIN) continue;
                        Fast_LogPlusEquals(beamstepM2[newi], m_state + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      node_value_type & prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                nucj, nucj1, seq_length);     
                        Fast_LogPlusEquals(beamstepC, prefix_C + state + score_external_paired/kT);  
                    } else {
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC, state + score_external_paired/kT);    
                    }
                }
            }
        }

        // beam of M2
        {
            for(int i = 0; i< j; i++) {
                node_value_type & state = beamstepM2[i];
                if (state == VALUE_MIN) continue;
                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp * seq_length + j];
                        if (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            newscore = - v_score_multi_unpaired(p+1, i-1) - v_score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(bestMulti[q][p], state + newscore/kT);  
                        }
                    }
                }

                // 2. M = M2
                {
                    Fast_LogPlusEquals(beamstepM[i], state);  
                }
            }
        }

        // beam of M
        {
            for(int i = 0; i< j; i++){
                node_value_type & state = beamstepM[i];

                if(state == VALUE_MIN) continue;

                if (j < seq_length-1) {
                    newscore = - v_score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(bestM[j+1][i], state + newscore/kT); 
                }
            }
        }

        // beam of C
        {
            if (j < seq_length-1) {
                value_type newscore;
                newscore = -v_score_external_unpaired(j+1, j+1);              
                Fast_LogPlusEquals(bestC[j+1], beamstepC + newscore/kT); 
            }
        }
    }  // end of for-loo j
    } // if there is no forest file
    // else  // read forest
      // load_forest();

    node_value_type & viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    if (is_verbose)
      printf("Partition Function Time: %.5f secs\n", parse_elapsed_time);   

    if(is_verbose) printf("Free Energy of Ensemble: %.2f kcal/mol\n", -kT * viterbi / 100.0);

    fflush(stdout);

    return {viterbi, 0, parse_elapsed_time};
}

// void BeamCKYParser::load_forest() {  
//   string line, type;
//   int i, j;
//   unordered_map<string, unordered_map<int, value_type> * > states;
//   states["P"] = bestP;
//   states["M"] = bestM;
//   states["M2"] = bestM2;
//   states["Multi"] = bestMulti;  

//   struct timeval forest_starttime, forest_endtime;

//   gettimeofday(&forest_starttime, NULL);
//   int num_nodes = 0;
//   while (true) {
//     if (!getline(cin, line)) break;
//     num_nodes ++;
//     istringstream iss(line);    
//     iss >> type;
//     if (type == "E") 
//       iss >> j >> bestC[j-1];
//     else  // P M M2 Multi
//       iss >> i >> j >> states[type][j-1][i-1];
//   }

//   gettimeofday(&forest_endtime, NULL);
//   double forest_elapsed_time = forest_endtime.tv_sec - forest_starttime.tv_sec + (forest_endtime.tv_usec-forest_starttime.tv_usec)/1000000.0;

//   printf("forest of %d nodes loaded in %.1f secs.\n", num_nodes, forest_elapsed_time);

// }

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             bool readforest,
                             string shape_file_path,
                             bool fasta)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      read_forest(readforest),
      is_fasta(fasta) {
    initialize();

    if (shape_file_path != "" ){
        use_shape = true;
        int position;
        string data;

        double temp_after_mb_shape;

        ifstream in(shape_file_path);

        if (!in.good()){
            cout<<"Reading SHAPE file error!"<<endl;
            assert(false);
        }

        // actually, we can combine the SHAPE_data and the energy_stack together
        while (!(in >> position >> data).fail()) {
            if (isdigit(int(data[0])) == 0){
                SHAPE_data.push_back(double((-1.000000)));
            }

            else {
                SHAPE_data.push_back(stod(data));
            }
        }

        for (int i = 0; i<SHAPE_data.size(); i++){
            temp_after_mb_shape = SHAPE_data[i] < 0 ? 0. : (m * log(SHAPE_data[i] + 1) + b);

            pseudo_energy_stack.push_back((int)roundf(temp_after_mb_shape * 100.));

            assert(pseudo_energy_stack.size() == i + 1 );
        }
    }
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}
// -------------------------------------------------------------

int main(int argc, char** argv){

    srand(std::chrono::system_clock::now().time_since_epoch().count()); // lhuang: not time(NULL)!

    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    int sample_number = 10;
    bool read_forest;
    bool fasta = false;
    string input_file;
    // SHAPE
    string shape_file_path = "";

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        sample_number = atoi(argv[4]);
        read_forest = atoi(argv[5]) == 1;
        fasta = atoi(argv[6]) == 1;
        input_file = argv[7];
        shape_file_path = argv[8];
    }

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    {
        ifstream infile(input_file);
        vector<string> input_lines;
        string line;
        if(input_file.size()){
            while (getline(infile, line))
                input_lines.push_back(line);
        }else{
            while (getline(cin, line))
                input_lines.push_back(line);
        }

        string rna_seq;
        vector<string> rna_seq_list, rna_name_list;
        if (fasta){
            for (string& seq : input_lines){
                if (seq.empty()) continue;
                else if (seq[0] == '>' or seq[0] == ';'){
                    rna_name_list.push_back(seq); // sequence name
                    if (!rna_seq.empty())
                        rna_seq_list.push_back(rna_seq);
                    rna_seq.clear();
                    continue;
                }else{
                    rtrim(seq);
                    rna_seq += seq;
                }
            }
            if (!rna_seq.empty())
                rna_seq_list.push_back(rna_seq);
        }else{
            for (string& seq : input_lines){
                if (seq.empty()) continue;
                if (!isalpha(seq[0])){
                    printf("Unrecognized sequence: %s\n", seq.c_str());
                    continue;
                }
                rna_seq_list.push_back(seq);
            }
        }

        for(int i = 0; i < rna_seq_list.size(); i++){
            if (rna_name_list.size() > i)
                printf("%s\n", rna_name_list[i].c_str());
            rna_seq = rna_seq_list[i];

            printf("%s\n", rna_seq.c_str());
            
            // convert to uppercase
            transform(rna_seq.begin(), rna_seq.end(), rna_seq.begin(), ::toupper);

            // convert T to U
            replace(rna_seq.begin(), rna_seq.end(), 'T', 'U');

            BeamCKYParser parser(beamsize, !sharpturn, is_verbose, read_forest, shape_file_path);

            parser.parse(rna_seq); // inside

            parser.sample(sample_number);
        }
    }

    if(is_verbose){
        gettimeofday(&total_endtime, NULL);
        double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
        printf("Total Time: %f secs\n", total_elapsed_time);
    }

    return 0;
}
