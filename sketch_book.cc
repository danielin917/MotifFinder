/*
 * Implementation for sketch book. A streaming motif finder.
 *
 */

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <unordered_set>

#include "kruskals.h"
#include "motif_finder.h"
#include "sketch_book.h"
#include "util.h"
#include "weight_matrix_model.h"

using namespace std;

namespace motif {

// We will buffer in 1 MB chunks.
static constexpr size_t kBufferSize = 1024 * 1024;
static constexpr size_t kSketchVecSize = 10000;
static constexpr int kSketchNumRows = 4;
static constexpr int kSketchRowSize = 1000;
static constexpr int kSampleSequenceVecSize = 200;

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

  sequence_vec_ = make_shared<vector<DNASequence>>();

  struct stat stat_buf;
  int rc = stat(filename.c_str(), &stat_buf);
  assert(rc != -1);
  long file_size = rc == 0 ? stat_buf.st_size : -1;
  cout << "file size: " << file_size << endl;
  const int desired_k = log(file_size) / log(4) + 1.0;
  cout << "Desired k: " << desired_k << endl;
  const int max_k = motif_length_ - num_mutations_  - 1;
  k_num_ = min(max_k, desired_k);
  cout << "Selected k: " << k_num_ << endl;

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
    string identifier;
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
          if (sequence_vec_->size() <  kSampleSequenceVecSize) {
            sequence_vec_->push_back(DNASequence(move(sequence), identifier));
          }
          sequence = "";
        }
        identifier = line;
        continue;
      }
      sequence += line;
    }
  }


  // Here we are going to run 1 EM step to filter good seeds.
  vector<CountMatrix> count_matrices;
  for (int ii = 0; ii < sketch_state_vec_.size(); ++ii) {
    CountMatrix count_matrix;
    for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides);
           ++ii) {
      // Add pseudo counts.
      count_matrix.push_back(vector<double>(motif_length_, 1));
    }

    for (int kk = 0; kk < sketch_state_vec_[ii].best_motif_word_deque.size();
          ++kk) {
      const string motif_string =
        sketch_state_vec_[ii].best_motif_word_deque[kk];

      for (int ii = 0; ii < motif_string.size(); ++ii) {
        const int n_index = static_cast<int>(
          CharToNucleotide(motif_string[ii]));
        count_matrix[n_index][ii] += 1;
      }
    }
    count_matrices.push_back(count_matrix);
  }

  MotifFinder mf(sequence_vec_, motif_length_, num_mutations_);
  vector<WeightMatrixModel> wmm_vec = BuildWMMVec(count_matrices);
  cout << "Performing seed EM on " << wmm_vec.size() << " models" << endl;
  for (int ii = 0; ii < wmm_vec.size(); ++ii) {
    wmm_vec[ii] = mf.RunEM(move(wmm_vec[ii]), 1);
  }
  sort(wmm_vec.begin(), wmm_vec.end());
  reverse(wmm_vec.begin(), wmm_vec.end());

  // Get the half of seeds with the best entropy.
/*
  vector<string> motif_vec;
  for (int ii = 0; ii < sketch_state_vec_.size(); ++ii) {
    for (int kk = 0; kk < sketch_state_vec_[ii].best_motif_word_deque.size();
          ++kk) {
      const string motif_string =
        sketch_state_vec_[ii].best_motif_word_deque[kk];
      motif_vec.push_back(motif_string);
    }
  }
*/
  vector<string> motif_vec;
  for (int ii = 0; ii < wmm_vec.size() / 2; ++ii) {
    motif_vec.push_back(wmm_vec[ii].GetConsensusString());
  }

  // Now we will perform clustering based on Hamming weight as the next level
  // of filtering.
  cout << "Clustering " << motif_vec.size() << " filtered motifs" << endl;
  Kruskals k(motif_vec.size());
  // First add all edges.
  for (int ii = 0; ii < motif_vec.size(); ++ii) {
    for (int kk = ii + 1; kk < motif_vec.size(); ++kk) {
      Edge e;
      e.start = ii;
      e.end = kk;
      e.weight = HammingDistance(motif_vec[ii], motif_vec[kk]);
      k.AddEdge(e);
    }
  }

  // Cluster to bring the number of EM models down by 80%
  static double kPercentThreshold = 0.5;
  cout << sketch_state_vec_.size() * kPercentThreshold << "++" << endl;
  const int num_clusters = motif_vec.size() * kPercentThreshold;
  const int avg_cluster_size = motif_vec.size() / num_clusters;
  assert(k.Cluster(motif_vec.size() * kPercentThreshold));


/*
  vector<int> best_nodes =
    k.GetNodesFromLargestClusters(motif_vec.size() * kPercentThreshold);
*/
  // Get mapping of each cluster.
  unordered_map<int, vector<int>> motif_forest;
  for (int ii = 0; ii < motif_vec.size(); ++ii) {
    int leader = k.GetLeader(ii);
    motif_forest[leader].push_back(ii);
  }

  static const int kMinClusterSize = 1;

  // We will now create a count matrix for each of our clusters.
  //vector<CountMatrix> count_matrices;
  count_matrices.clear();
  for (auto iter = motif_forest.begin(); iter != motif_forest.end(); ++iter) {
    const vector<int>& motif_index_vec = iter->second;
    if (motif_index_vec.size() < kMinClusterSize) {
      continue;
    }

    CountMatrix count_matrix;
    for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides);
           ++ii) {
      // Add pseudo counts.
      count_matrix.push_back(vector<double>(motif_length_, 0.25));
    }
    for (int ii = 0; ii < motif_index_vec.size(); ++ii) {
      const string& best_motif_word = motif_vec[motif_index_vec[ii]];
      cout << "Best motif " << ii <<  "  best_motif_word " <<  best_motif_word
           << endl;
      assert(best_motif_word.size() == motif_length_);
      for (int ii = 0; ii < best_motif_word.size(); ++ii) {
        const int n_index = static_cast<int>(
          CharToNucleotide(best_motif_word[ii]));
        count_matrix[n_index][ii] += 1;
      }
    }
    count_matrices.push_back(count_matrix);
  }

  cout << count_matrices.size() << " models being analyzed" << endl;


  //MotifFinder mf(sequence_vec_, motif_length_, num_mutations_);
  /*vector<WeightMatrixModel> */wmm_vec = BuildWMMVec(count_matrices);
  for (int ii = 0; ii < wmm_vec.size(); ++ii) {
    wmm_vec[ii] = mf.RunEM(move(wmm_vec[ii]), 10 /* max_iterations */);
  }

  cout << "B: " << endl;
  WeightMatrixModel best_wmm;
  assert(!wmm_vec.empty());
  sort(wmm_vec.begin(), wmm_vec.end());
  for (int ii = 0; ii < min(5, (int)wmm_vec.size()); ++ii) {
    cout << wmm_vec[wmm_vec.size() - 1 - ii].GetConsensusString() << ": "
         << wmm_vec[wmm_vec.size() - 1 - ii].entropy_history_vec.back()
         << endl;
  }
  best_wmm = wmm_vec.back();
  best_wmm.PrintFrequencyMatrix();
  cout << "Relative entropy: " << best_wmm.entropy_history_vec.back() << endl;
  cout << (1.0 * best_wmm.entropy_history_vec.back()) / motif_length_
       << " bits per letter" << endl;

  return best_wmm.GetConsensusString();
}

//-----------------------------------------------------------------------------

vector<WeightMatrixModel> SketchBook::BuildWMMVec(
  const vector<CountMatrix>& count_matrices) {

  // Build weight matrices.
  static const vector<double> kDefaultBackgroundFrequencyVec =
    { 0.25, 0.25, 0.25, 0.25 };
  vector<WeightMatrixModel> wmm_vec;
  for (int ii = 0; ii < count_matrices.size(); ++ii) {
    vector<vector<double>> frequency_matrix =
      MakeFrequencyMatrix(count_matrices[ii]);
    wmm_vec.push_back(
      WeightMatrixModel(frequency_matrix, kDefaultBackgroundFrequencyVec));
  }
  return wmm_vec;
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
      long min_count = sketch_state_vec_[ii].count_min.GetCount(projection);

      deque<string> *const best_motif_word_deque =
        &sketch_state_vec_[ii].best_motif_word_deque;
      if (min_count >= sketch_state_vec_[ii].current_best_min_count &&
           find(best_motif_word_deque->begin(),
                best_motif_word_deque->end(), motif_word) ==
             best_motif_word_deque->end()) {
        best_motif_word_deque->push_front(motif_word);
      }

      if (min_count > sketch_state_vec_[ii].current_best_min_count) {
        sketch_state_vec_[ii].current_best_min_count = min_count;
      }

      static const int kMaxSize = 3;
      if (best_motif_word_deque->size() > kMaxSize) {
        best_motif_word_deque->pop_back();
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
