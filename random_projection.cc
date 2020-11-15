/*
 *
 *
 */

#include <algorithm>
#include <iostream>

#include "random_projection.h"

using namespace std;

namespace motif {

//-----------------------------------------------------------------------------

RandomProjection::RandomProjection(
  const std::shared_ptr<std::vector<DNASequence>>& sequence_vec,
  const int motif_length,
  const int k_num) :
  sequence_vec_(sequence_vec),
  motif_length_(motif_length),
  k_num_(k_num) {
  assert(sequence_vec_);
  assert(k_num < motif_length);

  for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides); ++ii) {
    count_matrix_.push_back(vector<double>(motif_length_, 1));
  }
}

//-----------------------------------------------------------------------------

void RandomProjection::RunRandomProjectionScans(const int num_iterations) {
  for (int ii = 0; ii < num_iterations; ++ii) {
    Scan();
  }
}

//-----------------------------------------------------------------------------

void RandomProjection::Scan() {
  vector<int> shuffled_indices;
  for (int ii = 0; ii < motif_length_; ++ii) {
    shuffled_indices.push_back(ii);
  }

  Shuffle(shuffled_indices);

  vector<int> projection_indices;
  for (int ii = 0; ii < k_num_; ++ii) {
    projection_indices.push_back(shuffled_indices[ii]);
  }
  // Sort projection indices so we can store projection strings in=
  // understandable order.
  sort(projection_indices.begin(), projection_indices.end());
  unordered_map<string, vector<pair<int, int>>> projection_bucket_map;

  for (int ii = 0; ii < sequence_vec_->size(); ++ii) {
    const string& sequence = (*sequence_vec_)[ii].sequence();
    // Loop all starting indexes.
    for (int kk = 0; kk < sequence.size() - motif_length_ + 1; ++kk) {
      const string projection =
        GetProjection(sequence.substr(kk, motif_length_), projection_indices);
      projection_bucket_map_[projection].push_back(pair<int, int>(ii, kk));
    }
  }

  static int kDefaultThreshold = 2;
  for (auto iter = projection_bucket_map_.cbegin();
        iter != projection_bucket_map_.cend(); ++iter) {
    const vector<pair<int, int>>& location_vec = iter->second;
    if (location_vec.size() < kDefaultThreshold) {
      continue;
    }
    cout << "Locations of projection: " << iter->first << endl;
    for (int ii = 0; ii < location_vec.size(); ++ii) {
      cout << location_vec[ii].second << endl;
    }
  }
  UpdateCountMatrix();
}

//-----------------------------------------------------------------------------

void RandomProjection::UpdateCountMatrix() {

  static int kDefaultThreshold = 2;
  for (auto iter = projection_bucket_map_.cbegin();
        iter != projection_bucket_map_.cend(); ++iter) {
    const vector<pair<int, int>>& location_vec = iter->second;
    if (location_vec.size() < kDefaultThreshold) {
      continue;
    }
    for (int ii = 0; ii < location_vec.size(); ++ii) {
      const int sequence_no = location_vec[ii].first;
      const int motif_index = location_vec[ii].second;
      for (int kk = 0; kk < motif_length_; ++kk) {
        const int n_index = static_cast<int>(
          CharToNucleotide((*sequence_vec_)[sequence_no].
            sequence()[motif_index + kk]));
        count_matrix_[n_index][kk] += 1;
      }
    }
  }
}

//-----------------------------------------------------------------------------

WeightMatrixModel RandomProjection::GetSeedWMM() {
  static const vector<double> kDefaultBackgroundFrequencyVec =
    { 0.25, 0.25, 0.25, 0.25 };
  vector<vector<double>> frequency_matrix = MakeFrequencyMatrix(count_matrix_);
  return WeightMatrixModel(frequency_matrix, kDefaultBackgroundFrequencyVec);
}

//-----------------------------------------------------------------------------

string RandomProjection::GetProjection(const string& search_motif,
                                       const vector<int>& projection_indices) {
  string projection;
  for (int ii = 0; ii < projection_indices.size(); ++ii) {
    assert(projection_indices[ii] < search_motif.size());
    projection += toupper(search_motif[projection_indices[ii]]);
  }
  return projection;
}

//-----------------------------------------------------------------------------

}
