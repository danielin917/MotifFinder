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
  /* 'k' is the number of nucleotides that will be included in our proejction
   * scan. 'sequence_vec' contains the DNA sequences we will be scanning to
   * build our initial seed matrix.
   */
  RandomProjection(
    const std::shared_ptr<std::vector<DNASequence>>& sequence_vec,
    int motif_length,
    int k_num);

  // Run random projection scans for 'num_iterations'.
  void RunRandomProjectionScans(int num_iterations);

  // Return a seed WMM based on the count matrix.
  WeightMatrixModel GetSeedWMM();

 private:
  // Scan sequences using a random projection and update the 'count_matrix'
  // with heavy buckets.
  void Scan();

  // Calculate the k-projection of the 'search_motif' string using the k
  // indices provided in 'projection_indices'.
  std::string GetProjection(const std::string& search_motif,
                            const std::vector<int>& projection_indices);

  // Update 'count_matrix' with found bucketed motifs.
  void UpdateCountMatrix();

 private:
  // Vector of DNA sequences we will be traversing.
  std::shared_ptr<std::vector<DNASequence>> sequence_vec_;

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

  // Count matrix that holds our running counts for nucleotides at different
  // indexes find from random projection scan.
  std::vector<std::vector<double>> count_matrix_;
};

} // namespace

#endif // RANDOM_PROJECTION_H_
