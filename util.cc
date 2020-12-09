/*
 *
 *
 */

#include <cassert>
#include <iostream>
#include <vector>

#include "util.h"

using namespace std;

namespace motif {

//-----------------------------------------------------------------------------

Nucleotide CharToNucleotide(const char c) {
  const char upperc = toupper(c);
  switch (upperc) {
    case 'A':
      return Nucleotide::kA;
    case 'C':
      return Nucleotide::kC;
    case 'G':
      return Nucleotide::kG;
    case 'T':
      return Nucleotide::kT;
    default:
      cerr << "Invalid nucleotide char " << c << endl;
      assert(false);
      break;
  }
  __builtin_unreachable();
}

//-----------------------------------------------------------------------------

char NucleotideToChar(Nucleotide n) {
  switch (n) {
    case Nucleotide::kA:
      return 'A';
    case Nucleotide::kC:
      return 'C';
    case Nucleotide::kG:
      return 'G';
    case Nucleotide::kT:
      return 'T';
    default:
      cerr << "Invalid nucleotide " << static_cast<int>(n) << endl;
      assert(false);
      break;
  }
  __builtin_unreachable();
}

//-----------------------------------------------------------------------------

vector<vector<double>> MakeFrequencyMatrix(
  const vector<vector<double>>& count_matrix) {
  assert(count_matrix.size() == static_cast<int>(Nucleotide::kNumNucleotides));

  // Fill in the new frequency_matrix with the same dimensions as the count
  // matrix.
  vector<vector<double>> frequency_matrix;
  for (int ii = 0; ii < count_matrix.size(); ++ii) {
    frequency_matrix.push_back(vector<double>(count_matrix[ii].size(), 0));
  }

  for (int kk = 0; kk < count_matrix[0].size(); ++kk) {
    double total_count = 0;
    for (int ii = 0; ii < count_matrix.size(); ++ii) {
      total_count += count_matrix[ii][kk];
    }
    assert(total_count > 0);

    for (int ii = 0; ii < count_matrix.size(); ++ii) {
      frequency_matrix[ii][kk] = count_matrix[ii][kk] / total_count;
    }
  }
  return frequency_matrix;
}

//-----------------------------------------------------------------------------

double Entropy(const vector<vector<double>>& frequency_matrix,
                     const vector<double>& background_frequency_vec) {
  // Iterate over each column summing up relative entropy.
  double entropy = 0;
  for (int kk = 0; kk < frequency_matrix[0].size(); ++kk) {
    for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
      // For each column we are adding entropy based on the model probability
      // distribution vs the background probability distribution over the same
      // alphabet.
      // Relative entropy is given by H(P||Q):
      //  for all i in alphabet A SUM { P(xi) * log2(P(xi) / Q(xi)) }
      assert(frequency_matrix[ii][kk] <= 1);
      assert(background_frequency_vec[ii] > 0);
      entropy +=
        frequency_matrix[ii][kk] *
          log2(frequency_matrix[ii][kk] / background_frequency_vec[ii]);
    }
  }
  return entropy;
}

//-----------------------------------------------------------------------------

double log2(const double x) {
  return log(x) / log(2);
}

//-----------------------------------------------------------------------------

int BinomialCoeff(int n, int k) {
    double res = 1;
    for (int i = 1; i <= k; ++i) {
        res = res * (n - k + i) / i;
    }
    return (int)(res + 0.01);
}

//-----------------------------------------------------------------------------

string GetProjection(const string& search_motif,
                     const vector<int>& projection_indices) {
  string projection;
  projection.reserve(projection_indices.size());
  for (int ii = 0; ii < projection_indices.size(); ++ii) {
    #ifdef DEBUG
    assert(projection_indices[ii] < search_motif.size());
    #endif
    projection += search_motif[projection_indices[ii]];
  }
  return projection;
}

//-----------------------------------------------------------------------------

int HammingDistance(const std::string& s1, const std::string& s2) {
  #ifdef DEBUG
  assert(s1.size() == s2.size());
  #endif
  int distance = 0;
  for (int ii = 0; ii < s1.size(); ++ii) {
    if (s1[ii] != s2[ii]) {
      ++distance;
    }
  }
  return distance;
}

//-----------------------------------------------------------------------------

} // namespace
