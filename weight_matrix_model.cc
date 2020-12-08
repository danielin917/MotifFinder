/*
 *
 *
 */

#include <iostream>
#include <string>

#include "weight_matrix_model.h"
#include "util.h"

using namespace std;

namespace motif {

//-----------------------------------------------------------------------------

WeightMatrixModel::WeightMatrixModel() {
}

//-----------------------------------------------------------------------------

WeightMatrixModel::WeightMatrixModel(
  const std::vector<std::vector<double>>& p_frequency_matrix,
  const std::vector<double>& background_frequency_vec) :
  frequency_matrix(p_frequency_matrix) {

  // Build weight matrix using frequency matrix.
  for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
    matrix_model.push_back(vector<double>(frequency_matrix[0].size(), 0));
    for (int kk = 0; kk < frequency_matrix[ii].size(); ++kk) {
      assert(background_frequency_vec[ii] > 0);
      matrix_model[ii][kk] =
        log2(frequency_matrix[ii][kk] / background_frequency_vec[ii]);
    }
  }

  // TODO(DAN): Maybe do the entropy calculation inline above.
  entropy_history_vec.push_back(
    Entropy(frequency_matrix, background_frequency_vec));
}

//-----------------------------------------------------------------------------

void WeightMatrixModel::PrintWMM() const {
  assert(matrix_model.size() == 4);
  static const vector<char> kNucleotideChars = { 'A', 'C', 'G', 'T' };
  for (int ii = 0; ii < matrix_model.size(); ++ii) {
    cout << kNucleotideChars[ii] << ": ";
    for (int kk = 0; kk < matrix_model[ii].size(); ++kk) {
      cout << matrix_model[ii][kk] << " ";
    }
    cout << endl;
  }
}

//-----------------------------------------------------------------------------

void WeightMatrixModel::PrintFrequencyMatrix() const {
  assert(frequency_matrix.size() == 4);
  static const vector<char> kNucleotideChars = { 'A', 'C', 'G', 'T' };
  for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
    cout << kNucleotideChars[ii] << ": ";
    for (int kk = 0; kk < frequency_matrix[ii].size(); ++kk) {
      string frequency_string = to_string(frequency_matrix[ii][kk]);
      // Pad the back of short entropy numbers.
      if (frequency_string.size() < 8) {
        frequency_string += string(8 - (int)frequency_string.size(), ' ');
      }
      // Trim the end of larger entropy numbers.
      if (frequency_string.size() > 8) {
        frequency_string = frequency_string.substr(0, 8);
      }
      cout << frequency_string << ", ";
    }
    cout << endl;
  }
}

//-----------------------------------------------------------------------------

string WeightMatrixModel::GetConsensusString() const {
  string consensus_string;
  for (int kk = 0; kk < frequency_matrix[0].size(); ++kk) {
    double max_frequency = -1;
    Nucleotide n;
    for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
      if (frequency_matrix[ii][kk] > max_frequency) {
        n = static_cast<Nucleotide>(ii);
        max_frequency = frequency_matrix[ii][kk];
      }
    }
    consensus_string += NucleotideToChar(n);
  }
  return consensus_string;
}

//-----------------------------------------------------------------------------

}
