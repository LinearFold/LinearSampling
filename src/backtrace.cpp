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

default_random_engine generators[1];

void BeamCKYParser::backtrack_beamC(int *next_pair, int *prev_pair, int j, char* result){
    if(j == 0) return;
    auto &stateC = bestC[j];

    if (!stateC.visited){
        stateC.visited = true;

        float accu_alpha = 0.0;
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        // C = C + U
        accu_alpha = bestC[j-1].alpha - bestC[j].alpha;
        update_if_better(stateC, accu_alpha, MANNER_C_eq_C_plus_U);
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
                accu_alpha = bestC[k].alpha + bestP[j][i].alpha + score_external_paired/kT - bestC[j].alpha;

                update_if_better(stateC, accu_alpha, MANNER_C_eq_C_plus_P, k);
            } else {
                int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                         nucj, nucj1, seq_length);
                accu_alpha = bestP[j][0].alpha + score_external_paired/kT - bestC[j].alpha;
                update_if_better(stateC, accu_alpha, MANNER_C_eq_C_plus_P, -1);
            }
        }
        discrete_distribution<> d(stateC.alphalist.begin(), stateC.alphalist.end());
        stateC.distribution = d;
    
    }

    backtrack(next_pair, prev_pair, 0, j, result, stateC);
}

void BeamCKYParser::backtrack_beamP(int *next_pair, int *prev_pair, int i, int j, char* result){
    assert(i <= j - 4);
    auto &stateP = bestP[j][i];

    if (!stateP.visited){

        stateP.visited = true;

        float accu_alpha = 0.0;
        int newscore;

        int nuci = nucs[i];
        int nuci1 = (i+1) < seq_length ? nucs[i+1] : -1;
        int nucj = nucs[j];
        int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

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
                        accu_alpha = bestP[q][p].alpha + score_single/kT - bestP[j][i].alpha;

                        if (((p - i) == 1) && ((j - q) == 1)){
                            update_if_better(stateP, accu_alpha, MANNER_HELIX);
                        }
                        else{
                            update_if_better(stateP, accu_alpha, MANNER_SINGLE, static_cast<char>(p - i), j - q);
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
        accu_alpha = newscore/kT - bestP[j][i].alpha;
        update_if_better(stateP, accu_alpha, MANNER_HAIRPIN);
        auto iterator = bestMulti[j].find (i);
        if(iterator != bestMulti[j].end()) {
            int score_multi = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length);
            accu_alpha = bestMulti[j][i].alpha + score_multi/kT - bestP[j][i].alpha;
            update_if_better(stateP, accu_alpha, MANNER_P_eq_MULTI);
        }
        discrete_distribution<> d(stateP.alphalist.begin(), stateP.alphalist.end());
        stateP.distribution = d;
    }


    backtrack(next_pair, prev_pair, i, j, result, stateP);

}

void BeamCKYParser::backtrack_beamMulti(int *next_pair, int *prev_pair, int i, int j, char* result){
    auto &stateMulti = bestMulti[j][i];
    if (!stateMulti.visited){
        stateMulti.visited = true;
        float accu_alpha = 0.0;
        int nuci = nucs[i];
        int nuci1 = nucs[i+1];
        int jprev = prev_pair[nuci * seq_length + j];
        auto iterator = bestMulti[jprev].find (i);
        if (jprev > i+10 && iterator != bestMulti[jprev].end()) { // no sharp turn 
            accu_alpha = bestMulti[jprev][i].alpha - bestMulti[j][i].alpha;
            update_if_better(stateMulti, accu_alpha, MANNER_MULTI_JUMP, jprev);
        }

        for (int q = j-1; q >= i + 10; --q){
            for(auto& item : bestM2[q]){
                int p = item.first;
                if(p > i && (p - i) + (j - q) - 2 <= SINGLE_MAX_LEN){
                    accu_alpha = bestM2[q][p].alpha - bestMulti[j][i].alpha;

                    update_if_better(stateMulti, accu_alpha, MANNER_MULTI, static_cast<char>(p - i), j - q);
                }
            }
            if(_allowed_pairs[nuci][nucs[q]]) break;

        }
        discrete_distribution<> d(stateMulti.alphalist.begin(), stateMulti.alphalist.end());
        stateMulti.distribution = d;
    }
    backtrack(next_pair, prev_pair, i, j, result, stateMulti);

}

void BeamCKYParser::backtrack_beamM2(int *next_pair, int *prev_pair, int i, int j, char* result){

    auto &stateM2 = bestM2[j][i];

    if (!stateM2.visited){
        stateM2.visited = true;
        float accu_alpha = 0.0;

        int nuci = nucs[i];
        int nuci1 = nucs[i+1];
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

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
                
                    accu_alpha = bestM[m][i].alpha + bestP[j][k].alpha + M1_score/kT - bestM2[j][i].alpha;
                    update_if_better(stateM2, accu_alpha, MANNER_M2_eq_M_plus_P, m);
                }
            }
        }

        discrete_distribution<> d(stateM2.alphalist.begin(), stateM2.alphalist.end());
        stateM2.distribution = d;
    }
    backtrack(next_pair, prev_pair, i, j, result, stateM2);

}

void BeamCKYParser::backtrack_beamM(int *next_pair, int *prev_pair, int i, int j, char* result){
    auto &stateM = bestM[j][i];

    if (!stateM.visited){

        stateM.visited = true;

        float accu_alpha = 0.0;

        int nuci = nucs[i];
        int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        // M = M + U
        auto iterator = bestM[j-1].find(i);
        if(j > i+1 && iterator != bestM[j-1].end()) {
            accu_alpha = bestM[j-1][i].alpha - bestM[j][i].alpha;
            update_if_better(stateM, accu_alpha, MANNER_M_eq_M_plus_U);
        }

        iterator = bestP[j].find(i);
        if(iterator != bestP[j].end()) {
            int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
            accu_alpha = bestP[j][i].alpha + M1_score/kT - bestM[j][i].alpha;
            update_if_better(stateM, accu_alpha, MANNER_M_eq_P);
        }

        iterator = bestM2[j].find(i);
        if(iterator != bestM2[j].end()) {
            accu_alpha = bestM2[j][i].alpha - bestM[j][i].alpha;
            update_if_better(stateM,accu_alpha, MANNER_M_eq_M2);
        }
        discrete_distribution<> d(stateM.alphalist.begin(), stateM.alphalist.end());
        stateM.distribution = d;

    }
    backtrack(next_pair, prev_pair, i, j, result, stateM);

}

unsigned long choice(State& state){
    return state.distribution(generators[0]);
}



void BeamCKYParser::backtrack(int *next_pair, int *prev_pair, int i, int j, char* result, State& state){
    unsigned long index = choice(state);
    BackPointer backpointer = state.tracelist.at(index);
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
                backtrack_beamP(next_pair, prev_pair, p, q, result);
            }
            break;
        case MANNER_HELIX:
            {
                result[i] = '(';
                result[j] = ')';
                p = i + 1;
                q = j - 1;
                backtrack_beamP(next_pair, prev_pair, p, q, result);
            }
            break;
        case MANNER_MULTI: 
            p = i + backpointer.trace.paddings.l1;
            q = j - backpointer.trace.paddings.l2;
            backtrack_beamM2(next_pair, prev_pair, p, q, result);
            break;
        case MANNER_MULTI_JUMP:
            k = backpointer.trace.split;
            backtrack_beamMulti(next_pair, prev_pair, i, k, result);
            break;
        case MANNER_P_eq_MULTI:
            result[i] = '(';
            result[j] = ')';
            backtrack_beamMulti(next_pair, prev_pair, i, j, result);
            break;
        case MANNER_M2_eq_M_plus_P:
            k = backpointer.trace.split;
            backtrack_beamM(next_pair, prev_pair, i, k, result);
            backtrack_beamP(next_pair, prev_pair, k+1, j, result);
            break;
        case MANNER_M_eq_M2:
            backtrack_beamM2(next_pair, prev_pair, i, j, result);
            break;
        case MANNER_M_eq_M_plus_U:
            backtrack_beamM(next_pair, prev_pair, i, j-1, result);
            break;
        case MANNER_M_eq_P:
            backtrack_beamP(next_pair, prev_pair, i, j, result);
            break;
        case MANNER_C_eq_C_plus_U:
            k = j - 1;
            if (k != -1) backtrack_beamC(next_pair, prev_pair, k, result);
            break;
        case MANNER_C_eq_C_plus_P:
            {
                k = backpointer.trace.split;
                if (k != -1) {
                    backtrack_beamC(next_pair, prev_pair, k, result);
                    backtrack_beamP(next_pair, prev_pair, k+1, j, result);

                }
                else {
                    backtrack_beamP(next_pair, prev_pair, i, j, result);
                }
            }
            break;
        default:  // MANNER_NONE or other cases
            printf("wrong manner at %d, %d: manner %d\n", i, j, backpointer.manner); fflush(stdout);
            assert(false);
    }
    return;
}

void BeamCKYParser::sample(int *next_pair, int *prev_pair){
    char result[seq_length+1];

    srand (std::chrono::system_clock::now().time_since_epoch().count()); // lhuang: not time(NULL)!
    generators[0] = default_random_engine(rand());

    struct timeval parse_starttime, parse_endtime;
    gettimeofday(&parse_starttime, NULL);
    for(int i=0; i<sample_number; i++){
        memset(result, '.', seq_length);
        result[seq_length] = 0;

        backtrack_beamC(next_pair, prev_pair, seq_length-1, result);
        printf("%s\n", string(result).c_str());
    }   
    if(is_verbose){
        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        printf("Sequence_length: %d Sample Number: %d Sample Time: %f\n", seq_length, sample_number, parse_elapsed_time);
    }
    return;
}

