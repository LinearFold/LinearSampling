/*
 *backtrace.cpp*
Backtracing code for LinearSampling: Linear-Time RNA Secondary Structure Sampling.

 author: He Zhang, Liang Zhang
 edited by: 10/2020
*/

#include <iostream>
#include <random>
#include <array>
#include <chrono>
#include <time.h>
#include <stdlib.h>
#include <cassert>

#include "LinearSampling.h"

using namespace std;


void BeamCKYParser::backtrack_beamC(int j, char* result){
    if(j == 0) return;
    auto &stateC = bestC[j];

    SampleState & samplestate = sampleC[j];
    if (!samplestate.visited) {

      samplestate.visited = true;
      
	vector <float> alphalist;

        float accu_alpha = 0.0;
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        // C = C + U
	float localZ = bestC[j].alpha;
        accu_alpha = bestC[j-1].alpha - localZ;
        update_if_better(samplestate, alphalist, accu_alpha, MANNER_C_eq_C_plus_U);
        // C = C + P
        for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm 
            int i = item.first;
            State& state = item.second;
            int nuci = nucs[i];
            int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

            int k = i - 1;
            if (k >= 0) {
                int nuck = nuci_1;
                int nuck1 = nuci;
                int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                     nucj, nucj1, seq_length);
                accu_alpha = bestC[k].alpha + bestP[j][i].alpha + score_external_paired/kT - localZ;

                update_if_better(samplestate, alphalist, accu_alpha, MANNER_C_eq_C_plus_P, k);
            } else {
                int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                         nucj, nucj1, seq_length);
                accu_alpha = bestP[j][0].alpha + score_external_paired/kT - localZ;
                update_if_better(samplestate, alphalist, accu_alpha, MANNER_C_eq_C_plus_P, -1);
            }
        }
        
        samplestate.distribution = discrete_distribution<> (alphalist.begin(), alphalist.end());
    
    }

    backtrack(0, j, result, samplestate);
}

void BeamCKYParser::backtrack_beamP(int i, int j, char* result){
    assert(i <= j - 4);
    auto &stateP = bestP[j][i];
 
    SampleState & samplestate = sampleP[j][i];

    if (!samplestate.visited) {
      samplestate.visited = true;

	vector <float> alphalist;

        float accu_alpha = 0.0;
        int newscore;

        int nuci = nucs[i];
        int nuci1 = (i+1) < seq_length ? nucs[i+1] : -1;
        int nucj = nucs[j];
        int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

	float localZ = bestP[j][i].alpha;

        {
            int p, q, nucq, nucq1, nucp, nucp_1;
            // helix or single_branch
            for (q = j - 1; q >= std::max(j - SINGLE_MAX_LEN, i+5); --q) { // no sharp turn
                nucq = nucs[q];
                nucq1 = nucs[q + 1];
                p = next_pair[nucq * seq_length + i]; 
                while (p != -1 && p <= q - 4 && ((p - i) + (j - q) - 2 <= SINGLE_MAX_LEN)) {
                    auto iterator = bestP[q].find (p);
                    if(iterator != bestP[q].end()) {
                        nucp = nucs[p];
                        nucp_1 = nucs[p - 1];

                        int score_single = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                            nucp_1, nucp, nucq, nucq1); // same for vienna
                        accu_alpha = bestP[q][p].alpha + score_single/kT - localZ;

                        if (((p - i) == 1) && ((j - q) == 1)){
                            update_if_better(samplestate, alphalist, accu_alpha, MANNER_HELIX);
                        }
                        else{
                            update_if_better(samplestate, alphalist, accu_alpha, MANNER_SINGLE, static_cast<char>(p - i), j - q);
                        }
                    }
                    p = next_pair[nucq * seq_length + p];
                }
            }
        }

        // hairpin
        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
        if (j-i-1 == 4) // 6:tetra
            tetra_hex_tri = if_tetraloops[i];
        else if (j-i-1 == 6) // 8:hexa
            tetra_hex_tri = if_hexaloops[i];
        else if (j-i-1 == 3) // 5:tri
            tetra_hex_tri = if_triloops[i];
#endif

        newscore = - v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
        accu_alpha = newscore/kT - localZ;
        update_if_better(samplestate, alphalist, accu_alpha, MANNER_HAIRPIN);
	
        auto iterator = bestMulti[j].find (i);
        if(iterator != bestMulti[j].end()) {
            int score_multi = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length);
            accu_alpha = bestMulti[j][i].alpha + score_multi/kT - localZ;
            update_if_better(samplestate, alphalist, accu_alpha, MANNER_P_eq_MULTI);
        }
        samplestate.distribution = discrete_distribution<> (alphalist.begin(), alphalist.end());

    }
    backtrack(i, j, result, samplestate);

}

void BeamCKYParser::backtrack_beamMulti(int i, int j, char* result){
    auto &stateMulti = bestMulti[j][i];

    SampleState & samplestate = sampleMulti[j][i];
    
    if (!samplestate.visited) {

      samplestate.visited = true;
      
	float localZ = bestMulti[j][i].alpha;
	vector <float> alphalist;
	float accu_alpha = 0.0;
        int nuci = nucs[i];
        int nuci1 = nucs[i+1];
        int jprev = prev_pair[nuci * seq_length + j];
        auto iterator = bestMulti[jprev].find (i);
        if (jprev > i+10 && iterator != bestMulti[jprev].end()) { // no sharp turn 
            accu_alpha = bestMulti[jprev][i].alpha - localZ;
            update_if_better(samplestate, alphalist, accu_alpha, MANNER_MULTI_JUMP, jprev);
        }

        for (int q = j-1; q >= i + 10; --q){
            for(auto& item : bestM2[q]){
                int p = item.first;
                if(p > i && (p - i) + (j - q) - 2 <= SINGLE_MAX_LEN){
		  accu_alpha = bestM2[q][p].alpha - localZ;

                    update_if_better(samplestate, alphalist, accu_alpha, MANNER_MULTI, static_cast<char>(p - i), j - q);
                }
            }
            if(_allowed_pairs[nuci][nucs[q]]) break;

        }
        samplestate.distribution = discrete_distribution<> (alphalist.begin(), alphalist.end());

    }
    backtrack(i, j, result, samplestate);

}

void BeamCKYParser::backtrack_beamM2(int i, int j, char* result){

    auto &stateM2 = bestM2[j][i];
    SampleState & samplestate = sampleM2[j][i];

    if (!samplestate.visited) {

      samplestate.visited = true;
      
	  float accu_alpha = 0.0;

	vector <float> alphalist;
        int nuci = nucs[i];
        int nuci1 = nucs[i+1];
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
	
	float localZ = bestM2[j][i].alpha;
	
        // M2 = M + P
        for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm 
            int k = item.first;
            if (k > i) {
                int m = k - 1;
                auto iterator = bestM[m].find(i);
                if(iterator != bestM[m].end()) {
                    int nuck = nucs[k];
                    int nuck_1 = (k-1>-1) ? nucs[k-1] : -1;
                    value_type M1_score = - v_score_M1(k, j, j, nuck_1, nuck, nucj, nucj1, seq_length);
                
                    accu_alpha = bestM[m][i].alpha + bestP[j][k].alpha + M1_score/kT - localZ;
                    update_if_better(samplestate, alphalist, accu_alpha, MANNER_M2_eq_M_plus_P, m);
                }
            }
        }

        samplestate.distribution = discrete_distribution<> (alphalist.begin(), alphalist.end());
    }

    backtrack(i, j, result, samplestate);

}

void BeamCKYParser::backtrack_beamM(int i, int j, char* result){
    auto &stateM = bestM[j][i];

    SampleState & samplestate = sampleM[j][i];

    if (!samplestate.visited) {

      samplestate.visited = true;
      
	vector <float> alphalist;
        float accu_alpha = 0.0;

        int nuci = nucs[i];
        int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

	float localZ = bestM[j][i].alpha;

        // M = M + U
        auto iterator = bestM[j-1].find(i);
        if(j > i+1 && iterator != bestM[j-1].end()) {
            accu_alpha = bestM[j-1][i].alpha - localZ;
            update_if_better(samplestate, alphalist, accu_alpha, MANNER_M_eq_M_plus_U);
        }

        iterator = bestP[j].find(i);
        if(iterator != bestP[j].end()) {
            int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
            accu_alpha = bestP[j][i].alpha + M1_score/kT - localZ;
            update_if_better(samplestate, alphalist, accu_alpha, MANNER_M_eq_P);
        }

        iterator = bestM2[j].find(i);
        if(iterator != bestM2[j].end()) {
	  accu_alpha = bestM2[j][i].alpha - localZ;
            update_if_better(samplestate, alphalist, accu_alpha, MANNER_M_eq_M2);
        }
        samplestate.distribution = discrete_distribution<> (alphalist.begin(), alphalist.end());


    }

    backtrack(i, j, result, samplestate);

}

void BeamCKYParser::backtrack(int i, int j, char* result, SampleState& samplestate){

  BackPointer backpointer = samplestate.tracelist.at(samplestate.distribution(generator)); // lhuang: buggy: vector index out of range for empty tracelist
    int p, q, k;
    switch(backpointer.manner) {
        case MANNER_H:
            // this state should not be traced
            break;
        case MANNER_HAIRPIN:
            {
                result[i] = '(';
                result[j] = ')';
            }
            break;
        case MANNER_SINGLE:
            {
                result[i] = '(';
                result[j] = ')';
                p = i + backpointer.trace.paddings.l1;
                q = j - backpointer.trace.paddings.l2;
                backtrack_beamP(p, q, result);
            }
            break;
        case MANNER_HELIX:
            {
                result[i] = '(';
                result[j] = ')';
                p = i + 1;
                q = j - 1;
                backtrack_beamP(p, q, result);
            }
            break;
        case MANNER_MULTI: 
            p = i + backpointer.trace.paddings.l1;
            q = j - backpointer.trace.paddings.l2;
            backtrack_beamM2(p, q, result);
            break;
        case MANNER_MULTI_JUMP:
            k = backpointer.trace.split;
            backtrack_beamMulti(i, k, result);
            break;
        case MANNER_P_eq_MULTI:
            result[i] = '(';
            result[j] = ')';
            backtrack_beamMulti(i, j, result);
            break;
        case MANNER_M2_eq_M_plus_P:
            k = backpointer.trace.split;
            backtrack_beamM(i, k, result);
            backtrack_beamP(k+1, j, result);
            break;
        case MANNER_M_eq_M2:
            backtrack_beamM2(i, j, result);
            break;
        case MANNER_M_eq_M_plus_U:
            backtrack_beamM(i, j-1, result);
            break;
        case MANNER_M_eq_P:
            backtrack_beamP(i, j, result);
            break;
        case MANNER_C_eq_C_plus_U:
            k = j - 1;
            if (k != -1) backtrack_beamC(k, result);
            break;
        case MANNER_C_eq_C_plus_P:
            {
                k = backpointer.trace.split;
                if (k != -1) {
                    backtrack_beamC(k, result);
                    backtrack_beamP(k+1, j, result);

                }
                else {
                    backtrack_beamP(i, j, result);
                }
            }
            break;
        default:  // MANNER_NONE or other cases
            printf("wrong manner at %d, %d: manner %d\n", i, j, backpointer.manner); fflush(stdout);
            assert(false);
    }
    return;
}

void BeamCKYParser::sample(int sample_number){

  // sampleC = new SampleState[seq_length];

    sampleH = new unordered_map<int, SampleState>[seq_length];
    sampleP = new unordered_map<int, SampleState>[seq_length];
    sampleM = new unordered_map<int, SampleState>[seq_length];
    sampleM2 = new unordered_map<int, SampleState>[seq_length];
    sampleMulti = new unordered_map<int, SampleState>[seq_length];

  char result[seq_length+1];

    generator = default_random_engine(rand());

    struct timeval parse_starttime, parse_endtime;
    gettimeofday(&parse_starttime, NULL);
    for(int i = 0; i < sample_number; i++){
        memset(result, '.', seq_length);
        result[seq_length] = 0;

	try {
	  backtrack_beamC(seq_length-1, result);
	  printf("%s\n", string(result).c_str());
	  }
	catch (const out_of_range & err) {
//	  if (is_verbose)
//	    printf("empty tracelist\n");
	  i--; // NB: hit vector index out of range exception
	}
    }   
    if(is_verbose){
        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        printf("Sequence_length: %d Sample Number: %d Sample Time: %f secs\n", seq_length, sample_number, parse_elapsed_time);
    }
    cleanup();    
    return;
}

