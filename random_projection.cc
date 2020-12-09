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
  const std::shared_ptr<const std::vector<DNASequence>>& sequence_vec,
  const int motif_length,
  const int num_mutations) :
  sequence_vec_(sequence_vec),
  motif_length_(motif_length),
  num_mutations_(num_mutations) {
  assert(sequence_vec_);


  const int desired_k =
    log(sequence_vec_->size()*
      ((*sequence_vec_)[0].sequence().size() - motif_length_ + 1)) /
      log(4) + 1.0;
  cout << "Desired k: " << desired_k << endl;
  const int max_k = motif_length_ - num_mutations_  - 1;
  k_num_ = min(max_k, desired_k);
  cout << "Selected k: " << k_num_ << endl;
}

//-----------------------------------------------------------------------------

void RandomProjection::Scan() {
  const double num_motifs = sequence_vec_->size();

  // Set the bucket threshold to be twice the average bucket size.
  int bucket_threshold =
    2 * num_motifs *
      ((*sequence_vec_)[0].sequence().size() - motif_length_ + 1) /
        (pow(4.0, k_num_));
  bucket_threshold = max(4, bucket_threshold);
  //const int bucket_threshold = 3;
  cout << "Bucket threshold " << bucket_threshold << endl;

  static const double kBucketEnrichedProbability = 0.95;

  // We calculate the number of trials so the probability that in at least one
  // of the trials we will get an enriched bucket his is 0.95.
  const long double probability_planted_motif_hits =
    (1.0 * BinomialCoeff(motif_length_ - num_mutations_, k_num_)) /
      (1.0 * BinomialCoeff(motif_length_, k_num_));
  cout << "Prob given plant hits: " << probability_planted_motif_hits << endl;

  // Use Binomial distrbution to get probability that the threshold is not
  // reached by aggregating the probability of all hits less than threshold.
  double probability_threshold_miss = 0;
  for (int ii = 0; ii < bucket_threshold; ++ii) {
    // P(X = k) = Bin(n, k) * p^k * (1-p)^(n-k);
    probability_threshold_miss +=
      (1.0 * BinomialCoeff(num_motifs, ii)) *
        pow(probability_planted_motif_hits, ii) *
        pow(1.0 - probability_planted_motif_hits, num_motifs - ii);
  }
  cout << "Prob threshold miss: " << probability_threshold_miss << endl;

  // 1 - ProbThreshodMiss^m >= q
  // m = log(1 - q) / log(ProbThresholdMiss)
  const long num_trials =  probability_threshold_miss == 0 ?
                             1 : log(1.0 - kBucketEnrichedProbability) /
                                   log(probability_threshold_miss);

  cout << "Running " << num_trials << " trials for random projection" << endl;
  for (int ii = 0; ii < num_trials; ++ii) {
    projection_bucket_map_.clear();
    vector<int> shuffled_indices;
    for (int ii = 0; ii < motif_length_; ++ii) {
      shuffled_indices.push_back(ii);
    }

    Shuffle(shuffled_indices);

    vector<int> projection_indices;
    for (int ii = 0; ii < k_num_; ++ii) {
      projection_indices.push_back(shuffled_indices[ii]);
    }
    // Sort projection indices so we can store projection strings in
    // understandable order.
    sort(projection_indices.begin(), projection_indices.end());
    unordered_map<string, vector<pair<int, int>>> projection_bucket_map;

    for (int ii = 0; ii < sequence_vec_->size(); ++ii) {
      const string& sequence = (*sequence_vec_)[ii].sequence();
      assert(sequence.size() >  motif_length_);
      // Loop all starting indices.
      for (int kk = 0; kk < (int)sequence.size() - motif_length_ + 1; ++kk) {
        const string projection =
          GetProjection(sequence.substr(kk, motif_length_),
                        projection_indices);
        projection_bucket_map_[projection].push_back(pair<int, int>(ii, kk));
      }
    }
    GenerateCountMatrices(bucket_threshold);
  }
}

//-----------------------------------------------------------------------------

void RandomProjection::GenerateCountMatrices(const int bucket_threshold) {
  for (auto iter = projection_bucket_map_.cbegin();
        iter != projection_bucket_map_.cend(); ++iter) {
    const vector<pair<int, int>>& location_vec = iter->second;
    static int largest_bucket = 0;
    if (location_vec.size() > largest_bucket) {
      largest_bucket = location_vec.size();
      cout << "Largest bucket: " << largest_bucket << endl;
      continue;
    }

    if (location_vec.size() < bucket_threshold) {
      continue;
    }
    CountMatrix count_matrix;
    // We initialize with pseudocounts.
    for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides); ++ii) {
      count_matrix.matrix.push_back(vector<double>(motif_length_, 0.25));
    }
    count_matrix.total_count = location_vec.size();
    for (int ii = 0; ii < location_vec.size(); ++ii) {
      const int sequence_no = location_vec[ii].first;
      const int motif_index = location_vec[ii].second;
      for (int kk = 0; kk < motif_length_; ++kk) {
        const int n_index = static_cast<int>(
          CharToNucleotide((*sequence_vec_)[sequence_no].
            sequence()[motif_index + kk]));
        assert(kk < count_matrix.matrix[n_index].size());
        count_matrix.matrix[n_index][kk] += 1;
      }
    }
    count_matrices_.push_back(count_matrix);
  }
  sort(count_matrices_.begin(), count_matrices_.end());
  reverse(count_matrices_.begin(), count_matrices_.end());
}

//-----------------------------------------------------------------------------

vector<WeightMatrixModel> RandomProjection::GetWMMVec(const int max) {
  static const vector<double> kDefaultBackgroundFrequencyVec =
    { 0.25, 0.25, 0.25, 0.25 };
  vector<WeightMatrixModel> wmm_vec;
  for (int ii = 0; ii < count_matrices_.size(); ++ii) {
    vector<vector<double>> frequency_matrix =
      MakeFrequencyMatrix(count_matrices_[ii].matrix);
    wmm_vec.push_back(
      WeightMatrixModel(frequency_matrix, kDefaultBackgroundFrequencyVec));

    if (ii == max) {
      break;
    }
  }
  return wmm_vec;
}

//-----------------------------------------------------------------------------

}
