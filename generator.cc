/*
 *
 *
 */

#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <time.h>
#include <unordered_set>
#include <vector>

using namespace std;

//-----------------------------------------------------------------------------

static vector<char> kNucleotides  = { 'A', 'C', 'G', 'T' };

//-----------------------------------------------------------------------------

template<typename T>
T Shuffle(T permutation) {
  for (int ii = 0; ii < permutation.size(); ++ii) {
    int rand_offset = rand() % (permutation.size() - ii);
    swap(permutation[ii], permutation[ii + rand_offset]);
  }
  return permutation;
}

//-----------------------------------------------------------------------------

string CreateRandomSequence(int sequence_length);

//-----------------------------------------------------------------------------

void PrintSequenceWithMutations(stringstream& ss,
                                const string& motif_string,
                                int sequence_length,
                                int num_mutations);

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cout << "Incorret usage, expected: "
         << "./a.out <motif_string> <sequence_length> <num_mutations>" << endl;
    return -1;
  }

  string motif_string(argv[1]);
  const int sequence_length = atoi(argv[2]);
  if (sequence_length < motif_string.size()) {
    cout << "Invalid input, sequence must be at least the size of the motif"
         << endl;
    return -1;
  }
  const int num_mutations = atoi(argv[3]);

  cout << ">Motif-" << motif_string << ":" << "size="
       << motif_string.size() << endl;
  int motif_offset = sequence_length / 2;
  assert((sequence_length / 2) > motif_string.size());

  cout << ">Offset: " << sequence_length / 2<< endl;
  stringstream ss;
  srand(time(NULL));
  //int motif_offset = rand() % ((sequence_length - motif_string.length()) + 1);
  static int kNumSequencesGenerated = 20;
  for (int ii = 0; ii < 100; ++ii) {
    ss << ">Sequence-" << ii << ":\n";
    PrintSequenceWithMutations(ss, motif_string, sequence_length,
                               num_mutations);
  }
  cout << ss.str() << endl;
}

//-----------------------------------------------------------------------------

char RandomMutation(char c) {
  char m = c;
  while (m == c) {
    int rand_index = rand() % 4;
    m = kNucleotides[rand_index];
  }
  return m;
}

//-----------------------------------------------------------------------------

string CreateRandomSequence(int sequence_length) {
  stringstream ss;
  for (int ii = 0; ii < sequence_length; ++ii) {
    if (ii % 60 == 0 && ii != 0) {
      ss << '\n';
    }
    int rand_index = rand() % 4;
    ss << kNucleotides[rand_index];
  }
  return ss.str();
}

//-----------------------------------------------------------------------------

void PrintSequenceWithMutations(stringstream& ss,
                                const string& motif_string,
                                const int sequence_length,
                                const int num_mutations) {
  assert(num_mutations < motif_string.size());

  int motif_offset = sequence_length / 2;

  vector<int> shuffled_indices;
  for (int ii = 0; ii < motif_string.size(); ++ii) {
    shuffled_indices.push_back(ii);
  }
  shuffled_indices = Shuffle(shuffled_indices);

  unordered_set<int> mutated_indices;
  for (int ii = 0; ii < num_mutations; ++ii) {
    mutated_indices.insert(shuffled_indices[ii]);
  }
  string new_sequence = CreateRandomSequence(sequence_length);

  // Now overwrite the section where the motif is.
  int kk = motif_offset;
  int ii = 0;
  while (ii < motif_string.size()) {
    if (new_sequence[kk] == '\n') {
      ++kk;
      continue;
    }

    if (mutated_indices.count(ii) > 0) {
      new_sequence[kk] = RandomMutation(ii);
    } else {
      new_sequence[kk] = motif_string[ii];
    }
    ++kk;
    ++ii;
  }
  ss << new_sequence << '\n';
}

//-----------------------------------------------------------------------------
