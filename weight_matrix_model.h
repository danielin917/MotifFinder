/*
 *
 *
 */

#ifndef WEIGHT_MATRIX_MODEL_H_
#define WEIGHT_MATRIX_MODEL_H_

#include <vector>

namespace motif {

struct WeightMatrixModel {
  // Constructor.
  WeightMatrixModel();

  // Constructor, constructs weight matrix model based on the provided
  // 'frequency_matrix' and a 'background_frequency_vec'. The background
  // frequency std::vector should be of length 4 with one column for each
  // nucleotide. The relative entropy for the sequence is the sum of each
  // column's relative entropy.
  WeightMatrixModel(
    const std::vector<std::vector<double>>& frequency_matrix,
    const std::vector<double>& background_frequency_vec);

  // Less than operator used for caomparison.
  bool operator<(const WeightMatrixModel& other) const {
    assert(entropy_history_vec.size() > 0);
    assert(other.entropy_history_vec.size() > 0);
    return entropy_history_vec.back() < other.entropy_history_vec.back();
  }

  // Print the matrix model to stdout.
  void PrintWMM() const;

  // Print the frequency matrix to stdout.
  void PrintFrequencyMatrix() const;

  // Returns the string representing the most probablye nucleotide at each
  // index.
  std::string GetConsensusString() const;

  // Weight matrix model representing a 4 x K matrix of log2(f/b) where f
  // is the foreground frequency of a particular nucleotide at a particular
  // index and b is the expected background frequency.
  std::vector<std::vector<double>> matrix_model;

  // Frequency matrix associated with the WMM.
  std::vector<std::vector<double>> frequency_matrix;

  //  Relative entropy between the probability distribution of this matrix
  //  model compared to background probability distribution.
  std::vector<double> entropy_history_vec;
};

}

#endif // WEIGHT_MATRIX_MODEL_H_
