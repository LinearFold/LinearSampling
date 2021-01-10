/*
 *LinearSampling.h*
 header file for LinearSampling.cpp.

 author: He Zhang, Liang Zhang
 edited by: 12/2019
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include <random>
#include <set>

#include "Utils/logspace.h"

// #define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755


#ifdef FAST_FLOAT
  typedef float pf_type;
#else
  typedef double pf_type;
#endif

#define th -9.91152

using namespace std;

typedef int value_type;
#define VALUE_MIN numeric_limits<double>::lowest()


enum Manner {
  MANNER_NONE = 0,              // 0: empty
  MANNER_H,                     // 1: hairpin candidate
  MANNER_HAIRPIN,               // 2: hairpin
  MANNER_SINGLE,                // 3: single
  MANNER_HELIX,                 // 4: helix
  MANNER_MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
  MANNER_MULTI_JUMP,            // 6: multi = multi + ...
  MANNER_P_eq_MULTI,            // 7: P = (multi)
  MANNER_M2_eq_M_plus_P,        // 8: M2 = M + P
  MANNER_M_eq_M2,               // 9: M = M2
  MANNER_M_eq_M_plus_U,         // 10: M = M + U
  MANNER_M_eq_P,                // 11: M = P
  MANNER_C_eq_C_plus_U,     // 12: C = C + U
  MANNER_C_eq_C_plus_P,     // 13: C = C + P
};

union TraceInfo {
    int split;
    struct {
        char l1;
        int l2;
    } paddings;
};


struct BackPointer {
  Manner manner;
  TraceInfo trace;

  BackPointer(): manner(MANNER_NONE) {};
  BackPointer(Manner m): manner(m) {};

  void set(Manner manner_, int split_) {
      manner = manner_; 
      trace.split = split_;
  }

  void set(Manner manner_, char l1_, int l2_) {
      manner = manner_;
      trace.paddings.l1 = l1_; trace.paddings.l2 = l2_;
  }
};


struct State {
    float alpha;
    bool visited;

    discrete_distribution<> distribution;
    vector<BackPointer> tracelist;
    vector<float> alphalist;

    State(): alpha(VALUE_MIN), visited(false), distribution{}, tracelist{}, alphalist{} {};

};



class BeamCKYParser {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;
    int sample_number;
  bool read_forest;

    struct DecoderResult {
        State& viterbi;
        unsigned long num_states;
        double time;
    };

    BeamCKYParser(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false,
                  int sample_number=10,
		  bool read_forest=false);

    DecoderResult parse(string& seq);

private:

  void load_forest();

    void get_parentheses(char* result, string& seq);

    unsigned seq_length;

    unordered_map<int, State> *bestH, *bestP, *bestM2, *bestMulti, *bestM;

    State *bestC;

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    int *nucs;

    void prepare(unsigned len);

    void postprocess();

    void update_if_better(State &state, float alpha_, Manner manner_) {
      if (alpha_ > th) {
        state.tracelist.push_back(BackPointer(manner_));
        state.alphalist.push_back(Fast_Exp(alpha_));
      }
    };

    void update_if_better(State &state, float alpha_, Manner manner_, int split_) {
      if (alpha_ > th) {
        BackPointer backpointer;
        backpointer.set(manner_, split_);
        state.tracelist.push_back(backpointer);
        state.alphalist.push_back(Fast_Exp(alpha_));
      }
    };

    void update_if_better(State &state, float alpha_, Manner manner_, char l1_, int l2_) {
      if (alpha_ > th) {
        BackPointer backpointer;
        backpointer.set(manner_, l1_, l2_);
        state.tracelist.push_back(backpointer);
        state.alphalist.push_back(Fast_Exp(alpha_));
      }
    };

    float beam_prune(unordered_map<int, State>& beamstep);
    vector<pair<float, int>> scores;

    // sampling
    void sample(int *next_pair, int *prev_pair);
    void backtrack(int *next_pair, int *prev_pair, int i, int j, char* result, State& state);
    void fill_alphalist(State& state);

    void backtrack_beamC(int *next_pair, int *prev_pair, int j, char* result);
    void backtrack_beamP(int *next_pair, int *prev_pair, int i, int j, char* result);
    void backtrack_beamMulti(int *next_pair, int *prev_pair, int i, int j, char* result);
    void backtrack_beamM2(int *next_pair, int *prev_pair, int i, int j, char* result);
    void backtrack_beamM(int *next_pair, int *prev_pair, int i, int j, char* result);

};

#endif //FASTCKY_BEAMCKYPAR_H
