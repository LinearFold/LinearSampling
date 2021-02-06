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

#ifdef non_saving
// non-saving
BackPointer BeamCKYParser::recover_hyperedges(int i, int j, Type type){
  BackPointer backpointer;
  float sampled_alpha = ((float) rand() / (RAND_MAX));

  State & state = get_state(i, j, type);
  float localZ = state.alpha;
  
  float accu_alpha = 0.0;
  float temp_alpha;

  switch(type) {
    case TYPE_C: {
      int nucj = nucs[j];
      int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
    
      // C = C + U
      if (j == 0) return BackPointer(MANNER_C_eq_C_plus_U); // hzhang: N.B. j == 0

      temp_alpha = bestC[j-1].alpha - localZ;
      if (temp_alpha > th){
        accu_alpha += Fast_Exp(temp_alpha);
        if(accu_alpha > sampled_alpha) {
            return BackPointer(MANNER_C_eq_C_plus_U);
        }
      }
      // C = C + P
      for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm 
        int i = item.first;
        State& state = item.second;
        int nuci = nucs[i];
        int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

        int k = i - 1;
        int nuck = nuci_1;
        int nuck1 = nuci;
        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                    nucj, nucj1, seq_length);
        temp_alpha = ((k >= 0) ? bestC[k].alpha : 0) + state.alpha + score_external_paired/kT - localZ;
        if (temp_alpha > th){
          accu_alpha += Fast_Exp(temp_alpha);
          if(accu_alpha > sampled_alpha) {
            backpointer.set(MANNER_C_eq_C_plus_P, k);
            return backpointer;
          }
        }
      }
    }
    break;
    case TYPE_P: {
    int newscore;
    int nuci = nucs[i];
    int nuci1 = (i+1) < seq_length ? nucs[i+1] : -1;
    int nucj = nucs[j];
    int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

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
          temp_alpha = bestP[q][p].alpha + score_single/kT - localZ;
          if (temp_alpha > th){
            accu_alpha += Fast_Exp(temp_alpha);
            if(accu_alpha > sampled_alpha) {
              if (((p - i) == 1) && ((j - q) == 1))
                return BackPointer(MANNER_HELIX);
              else{
                backpointer.set(MANNER_SINGLE, static_cast<char>(p - i), j - q);
                return backpointer;
              }
            }
          }
        }
        p = next_pair[nucq * seq_length + p];
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
      temp_alpha = newscore/kT - localZ;
      if (temp_alpha > th){
        accu_alpha += Fast_Exp(temp_alpha);
        if(accu_alpha > sampled_alpha) {
          return BackPointer(MANNER_HAIRPIN);
        }
      }

      // Multiloop
      auto iterator = bestMulti[j].find (i);
      if(iterator != bestMulti[j].end()) {
        int score_multi = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length);
        temp_alpha = bestMulti[j][i].alpha + score_multi/kT - localZ;
        if (temp_alpha > th){
          accu_alpha += Fast_Exp(temp_alpha);
          if(accu_alpha > sampled_alpha) {
              return BackPointer(MANNER_P_eq_MULTI);
          }
        }
      }
    }
    break;
    case TYPE_M: {
      int nuci = nucs[i];
      int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;
      int nucj = nucs[j];
      int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

      // M = M + U
      auto iterator = bestM[j-1].find(i);
      if(j > i+1 && iterator != bestM[j-1].end()) {
        temp_alpha = bestM[j-1][i].alpha - localZ;
        if (temp_alpha > th){
          accu_alpha += Fast_Exp(temp_alpha);
          if(accu_alpha > sampled_alpha)
            return BackPointer(MANNER_M_eq_M_plus_U);
        }
      }

      // M = P
      iterator = bestP[j].find(i);
      if(iterator != bestP[j].end()) {
          int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
          temp_alpha = bestP[j][i].alpha + M1_score/kT - localZ;
          if (temp_alpha > th){
            accu_alpha += Fast_Exp(temp_alpha);
            if(accu_alpha > sampled_alpha)
              return BackPointer(MANNER_M_eq_P);
          }
      }

      // M = M2
      iterator = bestM2[j].find(i);
      if(iterator != bestM2[j].end()) {
        temp_alpha = bestM2[j][i].alpha - localZ;
          if (temp_alpha > th){
          accu_alpha += Fast_Exp(temp_alpha);
          if(accu_alpha > sampled_alpha)
            return BackPointer(MANNER_M_eq_M2);
        }
      }
    }
    break;
    case TYPE_M2: {
      int nuci = nucs[i];
      int nuci1 = nucs[i+1];
      int nucj = nucs[j];
      int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

      // M2 = M + P
      //for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm
      if(sortedP[j].size() == 0){
        for(auto const& item : bestP[j])
          sortedP[j].push_back(-item.first);
        sort(sortedP[j].begin(), sortedP[j].end());
      }
      
      for (auto & item : sortedP[j]) { // map not unorderd_map
        int k = -item;
        if (k > i + 4) { // lhuang: +4
            int m = k - 1;
            auto iterator = bestM[m].find(i);
            if(iterator != bestM[m].end()) {
                int nuck = nucs[k];
                int nuck_1 = (k-1>-1) ? nucs[k-1] : -1;
                value_type M1_score = - v_score_M1(k, j, j, nuck_1, nuck, nucj, nucj1, seq_length);
                temp_alpha = bestM[m][i].alpha + bestP[j][k].alpha + M1_score/kT - localZ;
                if (temp_alpha > th){
                  accu_alpha += Fast_Exp(temp_alpha);
                  if(accu_alpha > sampled_alpha) {
                      backpointer.set(MANNER_M2_eq_M_plus_P, m);
                    return backpointer;
                  }
                }
            }
        }
        else break;
      }
    }
    break;
    case TYPE_MULTI: {
      int nuci = nucs[i];
      int nuci1 = nucs[i+1];
      int jprev = prev_pair[nuci * seq_length + j];
      auto iterator = bestMulti[jprev].find (i);
      if (jprev > i+10 && iterator != bestMulti[jprev].end()) { // no sharp turn 
        temp_alpha = bestMulti[jprev][i].alpha - localZ;
        if (temp_alpha > th){
          accu_alpha += Fast_Exp(temp_alpha);
          if(accu_alpha > sampled_alpha) {
            backpointer.set(MANNER_MULTI_JUMP, jprev);
            return backpointer;
          }
        }
      }

      int p, q;
      for (q = j - 1; q >= jprev; q--) { 
        for (p = i+1; p <= q - 9 && (p - i) + (j - q) - 2 <= SINGLE_MAX_LEN; p++){
         // if (bestM2[q][p].alpha > -1e6) {
         auto bestM2_iter = bestM2[q].find(p);
         if(bestM2_iter != bestM2[q].end()){
          temp_alpha = bestM2[q][p].alpha - localZ;
          if (temp_alpha > th){
            accu_alpha += Fast_Exp(temp_alpha);
            if(accu_alpha > sampled_alpha) {
              backpointer.set(MANNER_MULTI, static_cast<char>(p - i), j - q);
              return backpointer;
            }
          }
         }  
       }
      }
    }
    break;
  }
  return BackPointer();
}
#else
void BeamCKYParser::recover_hyperedges(int i, int j, Type type, SampleState & samplestate) {
  samplestate.visited = true;
  uniq_visited ++; 

  State & state = get_state(i, j, type);
  float localZ = state.alpha;
  
  vector <float> alphalist;
  float edge_alpha;

  switch (type) {
  case TYPE_C: {
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
  
    // C = C + U
    if (j == 0) samplestate.append(alphalist, 0.0, MANNER_C_eq_C_plus_U); // hzhang: N.B. j == 0

    else{
      edge_alpha = bestC[j-1].alpha - localZ; // lhuang: j == 0?
      samplestate.append(alphalist, edge_alpha, MANNER_C_eq_C_plus_U);
      // C = C + P
      for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm 
        int i = item.first;
        State& state = item.second;
        int nuci = nucs[i];
        int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

        int k = i - 1;
        int nuck = nuci_1;
        int nuck1 = nuci;
        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                    nucj, nucj1, seq_length);
        edge_alpha = ((k >= 0) ? bestC[k].alpha : 0) + state.alpha + score_external_paired/kT - localZ; // lhuang?
        samplestate.append(alphalist, edge_alpha, MANNER_C_eq_C_plus_P, k);
      }
    }
  }
  break;
  case TYPE_P: {
    int newscore;
    int nuci = nucs[i];
    int nuci1 = (i+1) < seq_length ? nucs[i+1] : -1;
    int nucj = nucs[j];
    int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

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
          edge_alpha = bestP[q][p].alpha + score_single/kT - localZ;
          if (((p - i) == 1) && ((j - q) == 1))
            samplestate.append(alphalist, edge_alpha, MANNER_HELIX);
          else
            samplestate.append(alphalist, edge_alpha, MANNER_SINGLE, static_cast<char>(p - i), j - q);
        }
        p = next_pair[nucq * seq_length + p];
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
    edge_alpha = newscore/kT - localZ;
    samplestate.append(alphalist, edge_alpha, MANNER_HAIRPIN);

    // Multiloop
    auto iterator = bestMulti[j].find (i);
    if(iterator != bestMulti[j].end()) {
        int score_multi = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length);
        edge_alpha = bestMulti[j][i].alpha + score_multi/kT - localZ;
        samplestate.append(alphalist, edge_alpha, MANNER_P_eq_MULTI);
    }
  }
  break;
  case TYPE_M: {
    int nuci = nucs[i];
    int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

    // M = M + U
    auto iterator = bestM[j-1].find(i);
    if(j > i+1 && iterator != bestM[j-1].end()) {
      edge_alpha = bestM[j-1][i].alpha - localZ;
      samplestate.append(alphalist, edge_alpha, MANNER_M_eq_M_plus_U);
    }

    iterator = bestP[j].find(i);
    if(iterator != bestP[j].end()) {
        int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
        edge_alpha = bestP[j][i].alpha + M1_score/kT - localZ;
        samplestate.append(alphalist, edge_alpha, MANNER_M_eq_P);
    }

    iterator = bestM2[j].find(i);
    if(iterator != bestM2[j].end()) {
      edge_alpha = bestM2[j][i].alpha - localZ;
      samplestate.append(alphalist, edge_alpha, MANNER_M_eq_M2);
    }
  }
  break;
  case TYPE_M2: {
    int nuci = nucs[i];
    int nuci1 = nucs[i+1];
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

    if(sortedP[j].size() == 0){
      for(auto const& item : bestP[j])
        sortedP[j].push_back(-item.first);
      sort(sortedP[j].begin(), sortedP[j].end());
    }
        
    for (auto & item : sortedP[j]) { // map not unorderd_map
      int k = -item;
      if (k > i + 4) { // lhuang: +4
          int m = k - 1;
          auto iterator = bestM[m].find(i);
          if(iterator != bestM[m].end()) {
              int nuck = nucs[k];
              int nuck_1 = (k-1>-1) ? nucs[k-1] : -1;
              value_type M1_score = - v_score_M1(k, j, j, nuck_1, nuck, nucj, nucj1, seq_length);
              edge_alpha = bestM[m][i].alpha + bestP[j][k].alpha + M1_score/kT - localZ;
              samplestate.append(alphalist, edge_alpha, MANNER_M2_eq_M_plus_P, m);
          }
      }
      else break;
    }
  }
  break;
  case TYPE_MULTI: {
    int nuci = nucs[i];
    int nuci1 = nucs[i+1];
    int jprev = prev_pair[nuci * seq_length + j];
    auto iterator = bestMulti[jprev].find (i);
    if (jprev > i+10 && iterator != bestMulti[jprev].end()) { // no sharp turn 
      edge_alpha = bestMulti[jprev][i].alpha - localZ;
      samplestate.append(alphalist, edge_alpha, MANNER_MULTI_JUMP, jprev);
    }

    for (int q = j - 1; q >= jprev; q--) { 
      for (int p = i+1; p <= q - 9 && (p - i) + (j - q) - 2 <= SINGLE_MAX_LEN; p++){
       // if (bestM2[q][p].alpha > -1e6) {
       auto bestM2_iter = bestM2[q].find(p);
       if(bestM2_iter != bestM2[q].end()){
         edge_alpha = bestM2[q][p].alpha - localZ;
         samplestate.append(alphalist, edge_alpha, MANNER_MULTI, static_cast<char>(p - i), j - q);    
       }  
     }
    } 
   }
  }
  samplestate.distribution = discrete_distribution<> (alphalist.begin(), alphalist.end());
}
#endif

State & BeamCKYParser::get_state(int i, int j, Type type) {
  switch (type) {
  case TYPE_C:
    return bestC[j];
  case TYPE_P:
    return bestP[j][i];
  case TYPE_M:
    return bestM[j][i];
  case TYPE_M2:
    return bestM2[j][i];
  case TYPE_MULTI:
    return bestMulti[j][i];
  }
}
unordered_map<int, State> * BeamCKYParser::get_states(Type type) {
  switch (type) {
  case TYPE_C:
    assert(false);
  case TYPE_P:
    return bestP;
  case TYPE_M:
    return bestM;
  case TYPE_M2:
    return bestM2;
  case TYPE_MULTI:
    return bestMulti;
  }
}

int BeamCKYParser::backtrack(int i, int j, char* result, Type type){

#ifdef non_saving
  BackPointer backpointer = recover_hyperedges(i, j, type);
#else
  SampleState& samplestate = samplestates[type][j][i]; //get_sample_state(i, j, type);
  visited ++;
  if (!samplestate.visited)
    recover_hyperedges(i, j, type, samplestate);
  BackPointer backpointer = samplestate.tracelist.at(samplestate.distribution(generator)); // lhuang: buggy: vector index out of range for empty tracelist
#endif

  int p, q, k;
  switch(backpointer.manner) {
    case MANNER_HAIRPIN:
      result[i] = '(';
      result[j] = ')';
      break;
    case MANNER_SINGLE:
      result[i] = '(';
      result[j] = ')';
      p = i + backpointer.trace.paddings.l1;
      q = j - backpointer.trace.paddings.l2;
      backtrack(p, q, result, TYPE_P);
      break;
    case MANNER_HELIX:
      result[i] = '(';
      result[j] = ')';
      p = i + 1;
      q = j - 1;
      backtrack(p, q, result, TYPE_P);
      break;
    case MANNER_MULTI: 
      p = i + backpointer.trace.paddings.l1;
      q = j - backpointer.trace.paddings.l2;
      backtrack(p, q, result, TYPE_M2);
      break;
    case MANNER_MULTI_JUMP:
      k = backpointer.trace.split;
      backtrack(i, k, result, TYPE_MULTI);
      break;
    case MANNER_P_eq_MULTI:
      result[i] = '(';
      result[j] = ')';
      backtrack(i, j, result, TYPE_MULTI);
      break;
    case MANNER_M2_eq_M_plus_P:
      k = backpointer.trace.split;
      backtrack(i, k, result, TYPE_M);
      backtrack(k+1, j, result, TYPE_P);
      break;
    case MANNER_M_eq_M2:
      backtrack(i, j, result, TYPE_M2);
      break;
    case MANNER_M_eq_M_plus_U:
      backtrack(i, j-1, result, TYPE_M);
      break;
    case MANNER_M_eq_P:
      backtrack(i, j, result, TYPE_P);
      break;
    case MANNER_C_eq_C_plus_U:
      k = j - 1;
      if (k != -1) backtrack(0, k, result, TYPE_C);
      break;
    case MANNER_C_eq_C_plus_P:
      k = backpointer.trace.split;
      if (k != -1) //lhuang
        backtrack(0, k, result, TYPE_C);
      backtrack(k+1, j, result, TYPE_P);
      break;
    default:  // MANNER_NONE or other cases
      return -1;
  }
  return 0;
}

void BeamCKYParser::sample(int sample_number){

#ifndef non_saving
  samplestates = new unordered_map<int, SampleState> * [5]; //[seq_length];
  for (int i = 0; i < 5; i++)
    samplestates[i] = new unordered_map<int, SampleState>[seq_length];
#endif

  char result[seq_length+1];

  visited = 0, uniq_visited = 0;

  generator = default_random_engine(rand());

  int all_nodes = seq_length; // C's
  for (int j = 0; j < seq_length; j ++)
    all_nodes += bestH[j].size() + bestP[j].size() + bestM[j].size() + bestM2[j].size() + bestMulti[j].size();

  struct timeval starttime, endtime;

  // gettimeofday(&starttime, NULL);
  
  // gettimeofday(&endtime, NULL);
  // double recover_time = endtime.tv_sec - starttime.tv_sec + (endtime.tv_usec-starttime.tv_usec)/1000000.0;
  // fflush(stdout);

  gettimeofday(&starttime, NULL);
  // int count= 0;
  for(int i = 0; i < sample_number; i++){
    // count ++;
    memset(result, '.', seq_length);
    result[seq_length] = 0;

    try {
      int flag = backtrack(0, seq_length-1, result, TYPE_C);
      if(flag == -1) i--;
      else printf("%s\n", string(result).c_str());
    }
    catch (const out_of_range & err) {
      i--; // NB: hit vector index out of range exception
    }
  } 
  if(is_verbose){
      gettimeofday(&endtime, NULL);
      double sampling_time = endtime.tv_sec - starttime.tv_sec + (endtime.tv_usec-starttime.tv_usec)/1000000.0;
      printf("Sequence_length: %d Sample Number: %d Sample Time: %f secs  uniq_nodes: %d (%.2f%% of visits, %.2f%% of all nodes)\n",
       seq_length, sample_number, sampling_time, uniq_visited, uniq_visited * 100. / visited, uniq_visited * 100. / all_nodes);
  }
  fflush(stdout);
  cleanup();    
  return;
}
