/*
 * Implementation for sketch book. A streaming motif finder.
 *
 */

#include <iostream>
#include <fstream>
#include <random>
#include <sstream>

#include "sketch_book.h"
#include "util.h"

using namespace std;

namespace motif {

// We will buffer in 1 MB chunks.
static constexpr size_t kBufferSize = 1024 * 1024;
static constexpr size_t kSketchVecSize = 100;
static constexpr int kSketchNumRows = 4;
static constexpr int kSketchRowSize = 100000;

//-----------------------------------------------------------------------------

SketchBook::SketchState::SketchState(vector<int>&& p_projection_indices_vec,
                                     CountMinSketch&& p_count_min) :
  projection_indices_vec(p_projection_indices_vec),
  count_min(p_count_min),
  current_best_min_count(0) {

}

//-----------------------------------------------------------------------------

SketchBook::SketchBook(const string& filename,
                       const int motif_length,
                       const int num_mutations) :
  filename_(filename),
  motif_length_(motif_length),
  num_mutations_(num_mutations) {

  k_num_ = min(15, motif_length_ - num_mutations_  - 3);
  srand(time(NULL));

  vector<int> unshuffled_indices_vec;
  for (int ii = 0; ii < motif_length_; ++ii) {
    unshuffled_indices_vec.push_back(ii);
  }

  for (int ii = 0; ii < kSketchVecSize; ++ii) {
    vector<int> shuffled_indices_vec = unshuffled_indices_vec;
    Shuffle(shuffled_indices_vec);

    vector<int> projection_indices_vec;
    for (int ii = 0; ii < k_num_; ++ii) {
      projection_indices_vec.push_back(shuffled_indices_vec[ii]);
    }
    // Sort to get a readable order.
    sort(projection_indices_vec.begin(), projection_indices_vec.end());

    sketch_state_vec_.emplace_back(
      move(projection_indices_vec),
      CountMinSketch(kSketchNumRows, kSketchRowSize));
  }
}

//-----------------------------------------------------------------------------

// Scan FASTA file and look for motifs of length 'k_num_'.
// Returns the consensus string found during scan.
string SketchBook::Scan() {
  unique_ptr<char[]> buffer(new char[kBufferSize + 1]);

  // Add null at end of buffer so we can manipulate as a c-string.
  ifstream ifs(filename_.c_str());
  int num_blocks_read = 0;
  while (ifs) {
    //memset(buffer.get(), 0, kBufferSize + 1);
    // Read into buffer.
    ifs.read(buffer.get(), kBufferSize);
    ++num_blocks_read;
    cout << "Reading block " << num_blocks_read << endl;

    int index = 0;
    string sequence;
    while (index < kBufferSize) {
      string line;
      index = GetLine(buffer.get(), index, &line);
      //cout << line << endl;
      if (line.size() == 0) {
        // We must be at the end of a chunk.
        break;
      }

      if (line[0] == '>') {
        if (!sequence.empty()) {
          // The assumption is each sequence can fit into memory.
          ProcessSequence(sequence);
          sequence = "";
        }
        continue;
      }
      sequence += line;
    }
  }

  // We will now create a count matrix using our top performing strings.
  CountMatrix count_matrix;
  for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides); ++ii) {
    count_matrix.push_back(vector<double>(motif_length_, 0));
  }

  for (int ii = 0; ii < sketch_state_vec_.size(); ++ii) {
    const string& best_motif_word = sketch_state_vec_[ii].best_motif_word;
    cout << "Best motif " << ii <<  "  best_motif_word " <<  best_motif_word
         << endl;
    const vector<int>& projection_indices =
      sketch_state_vec_[ii].projection_indices_vec;
    assert(best_motif_word.size() == motif_length_);
    int projection_index = 0;
    for (int ii = 0; ii < motif_length_; ++ii) {
      const int n_index = static_cast<int>(
        CharToNucleotide(best_motif_word[ii]));
      if (projection_index < projection_indices.size() &&
          ii == projection_indices[projection_index]) {
        // We will double the weight of projection indices.
        count_matrix[n_index][ii] += 2;
        ++projection_index;
        continue;
      }
      count_matrix[n_index][ii] += 1;
    }
  }

  string consensus;
  for (int ii = 0; ii < count_matrix[0].size(); ++ii) {
    int best_count = -1;
    char best_char = 'x';
    for(int kk = 0; kk < count_matrix.size(); ++kk) {
      if (count_matrix[kk][ii] > best_count) {
        best_count = count_matrix[kk][ii];
        best_char = NucleotideToChar(static_cast<Nucleotide>(kk));
      }
    }
    assert(best_char != 'x');
    consensus += best_char;
  }
  return consensus;
}

//-----------------------------------------------------------------------------

void SketchBook::ProcessSequence(const string& sequence) {
  for (int kk = 0; kk < (int)sequence.size() - motif_length_ + 1; ++kk) {
    for (int ii = 0; ii < sketch_state_vec_.size(); ++ii) {
      const string motif_word = sequence.substr(kk, motif_length_);
      const string projection =
        GetProjection(motif_word,
                      sketch_state_vec_[ii].projection_indices_vec);
      sketch_state_vec_[ii].count_min.Increment(projection);
      int min_count = sketch_state_vec_[ii].count_min.GetCount(projection);
      if (min_count > sketch_state_vec_[ii].current_best_min_count) {
        // Since we are check the best count every update we should only ever
        // have to increment by 1 when the best improves.
        ++sketch_state_vec_[ii].current_best_min_count;
        sketch_state_vec_[ii].best_motif_word = motif_word;
      }
    }
  }
}

//-----------------------------------------------------------------------------

int SketchBook::GetLine(const char *const buffer,
                        const int index,
                        string *const line) {

  stringstream ss;
  int ii = index;

  // We will read until we hit our null terminator.
  while (buffer[ii] != '\0' && buffer[ii] != EOF) {
    if (buffer[ii] == '\n') {
      break;
    }
    ss << buffer[ii];
    ++ii;
  }
  *line = ss.str();

  // Return the next starting index;
  return ii + 1;
}

//-----------------------------------------------------------------------------

} // namespace
