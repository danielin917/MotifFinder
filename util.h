/*
 *
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <algorithm>
#include <vector>

namespace motif {

//-----------------------------------------------------------------------------

enum class Nucleotide {
  kA = 0 /* Adenine */,
  kC = 1 /* Cytosine */,
  kG = 2 /* Guanine */,
  kT = 3 /* Thymine */,
  kNumNucleotides = 4
};

//-----------------------------------------------------------------------------

// Convert the character 'c' to the enum Nucleotide type if possible.
Nucleotide CharToNucleotide(char c);

//-----------------------------------------------------------------------------

// Returns the character form of nucleotide 'n'.
char NucleotideToChar(Nucleotide n);

//-----------------------------------------------------------------------------

// Knuth-Shuffle
template<typename T>
void Shuffle(T& permutation) {
  for (int ii = 0; ii < permutation.size(); ++ii) {
    int rand_offset = rand() % (permutation.size() - ii);
    std::swap(permutation[ii], permutation[ii + rand_offset]);
  }
}

//-----------------------------------------------------------------------------

//  Builds and returns a 4 x k frequency matrix given a valid 4 x k
//  'count_matrix'.
std::vector<std::vector<double>> MakeFrequencyMatrix(
  const std::vector<std::vector<double>>& count_matrix);

//-----------------------------------------------------------------------------

// Calculate the relative entropy H(P||Q) of a 'frequency_matrix' relative to
// the given 'background_frequency_vec'. The background frequency std::vector
// should be of length 4 with one column for each nucleotide. The relative
// entropy for the sequence is the sum of each column's relative entropy.
double Entropy(const std::vector<std::vector<double>>& frequency_matrix,
               const std::vector<double>& background_frequency_vec);

//-----------------------------------------------------------------------------

// Returns the base 2 logarithm of 'x'.
double log2(double x);

//-----------------------------------------------------------------------------

// Returns value of Binomial Coefficient C(n, k)
int BinomialCoeff(int n, int k);

//-----------------------------------------------------------------------------

// Returns the projection of 'search_motif' based on the 'projection_indices'.
std::string GetProjection(const std::string& search_motif,
                          const std::vector<int>& projection_indices);

//-----------------------------------------------------------------------------

int HammingDistance(const std::string& s1, const std::string& s2);

//-----------------------------------------------------------------------------

}

#endif // UTIL_H_
