/*
 *
 *
 */

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <time.h>
#include <unordered_set>
#include <vector>

using namespace std;

const string kTestFileName = "test.fasta";
const string kTrainFileName = "train.fasta";

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

string CreateRandomSequence(int sequence_length,
                            bool with_break_lines = true);

//-----------------------------------------------------------------------------

void WriteTestDataFile(const string& filename,
                       const string& motif_string,
                       int num_mutations,
                       int sequence_length,
                       int num_sequences);

//-----------------------------------------------------------------------------

void PrintSequenceWithMutations(ofstream& ofs,
                                const string& motif_string,
                                int sequence_length,
                                int num_mutations);

//-----------------------------------------------------------------------------

int main(int argc, char *argv[]) {
  if (argc != 5) {
    cout << "Incorret usage, expected: "
         << "./a.out <motif_length> <sequence_length> <num_mutations> "
         << "<num_sequences>" << endl;
    return -1;
  }

  int motif_length = atoi(argv[1]);
  const int sequence_length = atoi(argv[2]);
  if (sequence_length < motif_length) {
    cout << "Invalid input, sequence must be at least the size of the motif"
         << endl;
    return -1;
  }
  string motif_string =
    CreateRandomSequence(motif_length, false /* with_line_breaks */);
  const int num_mutations = atoi(argv[3]);
  const int num_sequences = atoi(argv[4]);

  int motif_offset = sequence_length / 2;
  assert((sequence_length / 2) > motif_string.size());

  srand(time(NULL));
  WriteTestDataFile(kTestFileName,
                    motif_string,
                    num_mutations,
                    sequence_length,
                    num_sequences);
/*
  WriteTestDataFile(kTrainFileName,
                    motif_string,
                    num_mutations,
                    sequence_length,
                    num_sequences);
*/
}

//-----------------------------------------------------------------------------

void WriteTestDataFile(const string& filename,
                       const string& motif_string,
                       const int num_mutations,
                       const int sequence_length,
                       const int num_sequences) {

  ofstream ofs;
  ofs.open(filename.c_str(), std::ofstream::out);
  if (!ofs.good()) {
    cout << "Failed to open " << filename << endl;
  }
  ofs << ">Motif-" << motif_string << ":" << "size="
      << motif_string.size() << ":" << "num_mutations=" << num_mutations
      << '\n';
  ofs << ">Offset: " << sequence_length / 2<< endl;
  for (int ii = 0; ii < num_sequences; ++ii) {
    ofs << ">Sequence-" << ii << ":\n";
    PrintSequenceWithMutations(ofs, motif_string, sequence_length,
                               num_mutations);
    if (!ofs.good()) {
      cout << "Failed to write to " << filename << endl;
    }
  }
  ofs.close();
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

string CreateRandomSequence(int sequence_length,
                            bool with_line_breaks) {
  stringstream ss;
  for (int ii = 0; ii < sequence_length; ++ii) {
    if (with_line_breaks &&
        ii % 60 == 0 && ii != 0) {
      ss << '\n';
    }
    int rand_index = rand() % 4;
    ss << kNucleotides[rand_index];
  }
  return ss.str();
}

//-----------------------------------------------------------------------------

void PrintSequenceWithMutations(ofstream& ofs,
                                const string& motif_string,
                                const int sequence_length,
                                const int num_mutations) {
  assert(num_mutations < motif_string.size());

  int motif_offset = sequence_length / 2;

  // Add space to account for added endlines.
  motif_offset += motif_offset / 60;

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
    assert(kk < new_sequence.size());
    if (new_sequence[kk] == '\n') {
      ++kk;
      continue;
    }

    if (mutated_indices.count(ii) > 0) {
      new_sequence[kk] = RandomMutation(motif_string[ii]);
    } else {
      new_sequence[kk] = motif_string[ii];
    }
    ++kk;
    ++ii;
  }
  ofs << new_sequence << '\n';
}

//-----------------------------------------------------------------------------
