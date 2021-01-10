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

#include "LinearSampling.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "backtrace.cpp" // sampling
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

// in-place quick-select
float quickselect(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


float BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        float newalpha = (k >= 0 ? bestC[k].alpha : 0.0) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestC = new State[seq_length];

    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<int, State>[seq_length];
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];

    nucs = new int[seq_length];

    scores.reserve(seq_length); 
}

void BeamCKYParser::postprocess() {
    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
}


BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq) {
    struct timeval parse_starttime, parse_endtime;

    // number of states
    unsigned long nos_H = 0, nos_P = 0, nos_M2 = 0,
            nos_M = 0, nos_C = 0, nos_Multi = 0;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    int *next_pair = new int[NOTON * seq_length]{-1};
    int *prev_pair = new int[NOTON * seq_length]{-1};
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

    if(seq_length > 0) bestC[0].alpha = 0.0;
    if(seq_length > 1) bestC[1].alpha = 0.0;
    ++nos_C;

    value_type newscore;
    // from left to right
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);
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
                    Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
                    ++ nos_H;
                }
            }

            {
                for (auto &item : beamstepH) {
                    int i = item.first;
                    State &state = item.second;
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
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
                        ++nos_H;
                    }

                    // 2. generate p(i, j)
                    {
                        Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);
                        ++ nos_P;
                    }
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci * seq_length + j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
                        // 1. extend (i, j) to (i, jnext)
                        
                        newscore = - v_score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha + newscore/kT);
                        ++nos_Multi;
                    }
                }

                // 2. generate P (i, j)
                {
                    value_type score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + score_multi/kT);
                    ++ nos_P;
                }
            }
        }

        // beam of P
        {
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                if (i >0 && j<seq_length-1) {
                    // value_type precomputed;
                    // precomputed = 0;
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
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + score_single/kT);
                                ++nos_P;
                            } else {
                                // single branch
                                int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + score_single/kT);
                                ++nos_P;
                            }
                            q = next_pair[nucp * seq_length + q];
                        }
                    }
                }
                // printf(" helix / single at %d\n", j); fflush(stdout);

                // 2. M = P
                if(i > 0 && j < seq_length-1){
                    int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + score_M1/kT);
                    ++ nos_M;
                }
                // printf(" M = P at %d\n", j); fflush(stdout);

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
                    value_type M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = state.alpha + M1_score/kT;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                        ++nos_M2;
                    }
                }
                // printf(" M/M2 = M + P at %d\n", j); fflush(stdout);

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                nucj, nucj1, seq_length);     
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + score_external_paired/kT);  
                        ++ nos_C;
                    } else {
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + score_external_paired/kT);    
                        ++ nos_C;
                    }
                }
            }
        }
        // printf("P at %d\n", j); fflush(stdout);

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp * seq_length + j];
                        if (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            newscore = - v_score_multi_unpaired(p+1, i-1) - v_score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha + newscore/kT);  
                            ++ nos_Multi;
                        }
                    }
                }

                // 2. M = M2
                {
                    Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
                    ++ nos_M;
                }
            }
        }
        // printf("M2 at %d\n", j); fflush(stdout);

        // beam of M
        {
            if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
                    newscore = - v_score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha + newscore/kT); 
                    ++ nos_M;
                }
            }
        }
        // printf("M at %d\n", j); fflush(stdout);

        // beam of C
        {
            if (j < seq_length-1) {
                value_type newscore;
                newscore = -v_score_external_unpaired(j+1, j+1);              
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha + newscore/kT); 
                ++ nos_C;
            }
        }
    }  // end of for-loo j
    } // if there is no forest file
    else  // read forest
      load_forest();

    State& viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;

    if(is_verbose) printf("Free Energy of Ensemble: %.2f kcal/mol\n", -kT * viterbi.alpha / 100.0);

    fflush(stdout);

    sample(next_pair, prev_pair);

    postprocess();

    return {viterbi, nos_tot, parse_elapsed_time};
}

void BeamCKYParser::load_forest() {
  string line, type, ii, jj;
  int i, j;
  string alpha;
  unordered_map<int, State> *states;
  while (true) {
    if (!getline(cin, line)) break;
    istringstream iss(line);
    iss >> type;
    if (type == "E") {
      iss >> jj >> alpha;
      j = stoi(jj) - 1;
      bestC[j].alpha = stof(alpha);	
    }
    else { // P M M2 Multi
      iss >> ii >> jj >> alpha;
      i = stoi(ii) - 1;
      j = stoi(jj) - 1;
      if (type == "P")
	states = bestP;
      else if (type == "M")
	states = bestM;
      else if (type == "M2")
	states = bestM2;
      else if (type == "Multi")
	states = bestMulti; 
      	
      states[j][i].alpha = stof(alpha);
    }      
  }  
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             int samplenumber,
			     bool readforest)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      sample_number(samplenumber),
      read_forest(readforest) {
#ifdef lv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif
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

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        sample_number = atoi(argv[4]);
	read_forest = atoi(argv[5]) == 1;
    }

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    {
        for (string seq; getline(cin, seq);) {
            if (seq.length() == 0)
                continue;

            if (seq[0] == ';' || seq[0] == '>') {
                printf("%s\n", seq.c_str());
                continue;
            }

            if (!isalpha(seq[0])){
                printf("Unrecognized sequence: %s\n", seq.c_str());
                continue;
            }

            printf("%s\n", seq.c_str());
            
            // convert to uppercase
            transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

            // convert T to U
            replace(seq.begin(), seq.end(), 'T', 'U');

            BeamCKYParser parser(beamsize, !sharpturn, is_verbose, sample_number, read_forest);

            BeamCKYParser::DecoderResult result = parser.parse(seq);
        }
    }

    if(is_verbose){
        gettimeofday(&total_endtime, NULL);
        double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
        printf("Total Time: %f\n", total_elapsed_time);
    }

    return 0;
}
