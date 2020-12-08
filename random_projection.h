/*
 * This file describes the random projection class which implements the
 * k-random pojection method for finding motifs as described in:
 * Finding Motifs Using Random Projections
 *     Jeremy Buhler and Martin Tompa
 *
 */

#ifndef RANDOM_PROJECTION_H_
#define RANDOM_PROJECTION_H_

#include <unordered_map>
#include <vector>

//#include "seqan/sequence_journaled.h"
#include "dna_sequence.h"
#include "util.h"
#include "weight_matrix_model.h"

namespace motif {

class RandomProjection {
 public:
  typedef std::vector<std::vector<double>> CountMatrix;

  /* 'k' is the number of nucleotides that will be included in our proejction
   * scan. 'sequence_vec' contains the DNA sequences we will be scanning to
   * build our initial seed matrix.
   */
  RandomProjection(
    const std::shared_ptr<const std::vector<DNASequence>>& sequence_vec,
    int motif_length,
    int num_mutations);

  // Scan sequences using a random projection and update the 'count_matrix'
  // with heavy buckets.
  void Scan();

  // Return the weight matrix models associated with each bucket of l-mers.
  std::vector<WeightMatrixModel> GetWMMVec();

 private:
  // Create and add a count matrix per bucket.
  void GenerateCountMatrices(int bucket_threshold);

 private:
  // Vector of DNA sequences we will be traversing.
  std::shared_ptr<const std::vector<DNASequence>> sequence_vec_;

  // The 'k' number we will use to generate projections.
  int motif_length_;

  // Mapping from a k-projection string to it's occurrence in the search
  // sequences as a pair <sequence_no, index>.
  std::unordered_map<std::string, std::vector<std::pair<int, int>>>
    projection_bucket_map_;

  // The size of the random projection. If set so that 4^k > t(n - l + 1) then
  // the average bucket will contain less than 1 entry. This should also be set
  // less than l - d so that we may get hits on the planted motif even with d
  // mutations.
  int k_num_;

  // The number of expected mutations in the motif.
  int num_mutations_;

  // Count matrix that holds our running counts for nucleotides at different
  // indexes find from random projection scan.
  std::vector<CountMatrix> count_matrices_;
};

} // namespace

#endif // RANDOM_PROJECTION_H_
