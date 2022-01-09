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

#ifdef FAST_FLOAT
  typedef float pf_type;
#else
  typedef double pf_type;
#endif

#define kT 61.63207755
#define th -9.91152

using namespace std;

typedef int value_type;
#define VALUE_MIN numeric_limits<float>::lowest()

enum Type {
  TYPE_C = 0,
  TYPE_P,
  TYPE_M,
  TYPE_M2,
  TYPE_MULTI,
};

// enum Saving {
//   SAVING_NON = 0,
//   SAVING_LAZY = 1,
//   SAVING_FULL = 2,
// };

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
  MANNER_C_eq_C_plus_U,         // 12: C = C + U
  MANNER_C_eq_C_plus_P,         // 13: C = C + P
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

//#define TEST 
struct State {
  float alpha;
  State(): alpha(VALUE_MIN)
  {};
};

#ifndef non_saving
struct SampleState {
  bool visited;
  discrete_distribution<> distribution;
  vector<BackPointer> tracelist;
  SampleState(): visited(false), distribution{}, tracelist{}
  {};

void append(vector<float> & alphalist, float alpha_, Manner manner_) {
  if (alpha_ > th) {
    tracelist.push_back(BackPointer(manner_));
    alphalist.push_back(Fast_Exp(alpha_));
  }
};

void append(vector<float> & alphalist, float alpha_, Manner manner_, int split_) {
  if (alpha_ > th) {
    BackPointer backpointer;
    backpointer.set(manner_, split_);
    tracelist.push_back(backpointer);
    alphalist.push_back(Fast_Exp(alpha_));
  }
};

void append(vector<float> & alphalist, float alpha_, Manner manner_, char l1_, int l2_) {
    if (alpha_ > th) {
      BackPointer backpointer;
      backpointer.set(manner_, l1_, l2_);
      tracelist.push_back(backpointer);
      alphalist.push_back(Fast_Exp(alpha_));
    }
  };
};
#endif

class BeamCKYParser {
public:
  int beam;
  bool no_sharp_turn;
  bool is_verbose;
  int sample_number;
  bool read_forest;
  bool is_fasta;

  struct DecoderResult {
      State& viterbi;
      unsigned long num_states;
      double time;
  };

  BeamCKYParser(int beam_size=100,
                bool nosharpturn=true,
                bool is_verbose=false,
	              bool read_forest=false,
                bool is_fasta=false);

  DecoderResult parse(string& seq);
  void sample(int sample_number);

private:
  default_random_engine generator; // for choice
  //int *next_pair, *prev_pair;
  vector<int> next_pair;
  vector<int> prev_pair;

  void load_forest();
  void get_parentheses(char* result, string& seq);

  unsigned seq_length;

  unordered_map<int, State> *bestH, *bestP, *bestM2, *bestMulti, *bestM;
  vector<int> *sortedP;
  State *bestC;

  vector<int> if_tetraloops;
  vector<int> if_hexaloops;
  vector<int> if_triloops;

  int *nucs;

  void prepare(unsigned len);

  void cleanup();

  float beam_prune(unordered_map<int, State>& beamstep);
  vector<pair<float, int>> scores;

  // sampling
  int backtrack(int i, int j, char* result, Type type);
  
  State & get_state(int i, int j, Type type);
  unordered_map<int, State> * get_states(Type type);
  
#ifdef non_saving
  BackPointer recover_hyperedges(int i, int j, Type type);
#else
  SampleState & get_sample_state(int i, int j, Type type);
  void recover_hyperedges(int i, int j, Type type, SampleState & samplestate);
  unordered_map<int, SampleState> ** samplestates;
#endif

  int visited = 0, uniq_visited = 0;
  // int saving_option = SAVING_FULL;

};

#endif //FASTCKY_BEAMCKYPAR_H
