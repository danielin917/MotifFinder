/*
 * This file provides the interface for the MotifFinder class.
 *
 * Author: Dan Lin danielin@uw.edu
 */

#ifndef MOTIF_FINDER_H_
#define MOTIF_FINDER_H_

#include <unordered_map>

#include "dna_sequence.h"
#include "meme_graph.h"
#include "random_projection.h"
#include "weight_matrix_model.h"

namespace motif {

// Max number of iterations that EM will be run for.
static const int kMaxNumIterations = 1000;

class MotifFinder {
 public:
  // Construct MotifFinder class using provided 'sequence_vec' and the expected
  // 'k_length' of the motif we are looking for.
  MotifFinder(std::shared_ptr<std::vector<DNASequence>> sequence_vec,
              int motif_length,
              int num_mutations);

  // Test the current trained models against the sequences in
  // 'test_sequence_vec'.
  void Run(const std::vector<DNASequence>& test_sequence_vec);

 private:
  // Train the MEME model.
  void TrainMemeModel();

  // Train the random projection model.
  void TrainRandomProjectionModel();

  // Create a seed matrix model based on 'seed_kmer'.
  WeightMatrixModel CreateSeedMatrixModel(const std::string& seed_kmer);

  // Create an initial seed for EM using the random projection paradigm.
  WeightMatrixModel CreateRandomProjectionSeed();

  // Run expectation maximization using a staring 'wmm' as a starting model.
  // Run for 'num_iterations' steps.
  WeightMatrixModel RunEM(WeightMatrixModel&& wmm,
                          const int num_iterations = kMaxNumIterations);

  // Builds and returns a 4 x k count matrix based on counting the occurrences
  // of each nucleotide at each position in the k-mer for the given sequences.
  // k being the expected motif length. 'yhat_matrix' represents yhat values
  // that will weight each count in the count matrix. Pass in an empty
  // yhat_matrix if no weight is being provided.
  std::vector<std::vector<double>> makeCountMatrix(
    const std::vector<std::vector<double>>& yhat_matrix);

  // CSEP 527 MEME subroutines ////////////////////////////////////////////////

  // Updates the 4 x k 'count_matrix' with the k-mers present in the given
  // 'sequence_str' where k is the expected motif length. 'yhat_vec' contains
  // the yhat values associated with the k-mer starting at each corresponding
  // index in the sequence_str. If yhat_vec is empty then we will assume that
  // all indices are equally weighted when counting.
  void updateCountMatrix(const std::string sequence_str,
                         const std::vector<double>& yhat_vec,
                         std::vector<std::vector<double>> *count_matrix);

  // Add pseudo counts in 'pseudo_count_vec' to corresponding nucleotide counts
  // in the 'count_matrix'. 'pseud_count_vec' should be of length 4 with one
  // slot for each nucleotide.
  void addPseudo(const std::vector<double>& pseudo_count_vec,
                 std::vector<std::vector<double>> *const count_matrix);

  // Updates 'wmm' using 'frequency_matrix' and 'background_frequency_vec'. The
  // update will preserve the entropy history for the matrix.
  void updateMatrixModel(
    const std::vector<std::vector<double>>& frequency_matrix,
    const std::vector<double>& background_frequency_vec,
    WeightMatrixModel *wmm);

  // Scan/score each position of each sequence (excluding those < k from the
  // rightmost end) using the current weight matrix model 'wmm'.
  std::vector<std::vector<double>> scanWMM(
    const WeightMatrixModel& wmm,
    const std::vector<DNASequence>& sequence_vec);

  // Using the current weight matrix model run the E-step of MEME's EM
  // algorithm on each sequencs; i.e., what is E[zij], where zij is the
  // zero-one variable indicating whether the motif instance in sequence i
  // begins in position j. Matrix returned will be the corresponding
  // probabilities for each valid index in each provided sequence.
  std::vector<std::vector<double>> Estep(const WeightMatrixModel& wmm);

  // Mstep: given the Estep result, re-estimate and update the 'wmm'. (This is
  // similar to makeCountMatrix / addPseudo / makeFrequencyMatrix / makeWMM,
  // with the additional wrinkle that the k-mer inputs to makeCountMatrix each
  // have weights. 'yhat_matrix' is the probability at each index that that
  // index contains the motif calculated from the E-step. 'pseudo_count_vec' is
  // the baseline number of counts that will be given to each nucleotide.
  // 'background_frequency_vec' is the expected background frequency for each
  // nucleotide in the order A, C, G, T.
  void Mstep(const std::vector<std::vector<double>>& yhat_matrix,
             const std::vector<double>& pseudo_count_vec,
             const std::vector<double>& background_frequency_vec,
             WeightMatrixModel *wmm);

  /////////////////////////////////////////////////////////////////////////////

  // Scores the test sequences in 'test_sequence_vec' using 'wmm' and prints
  // the counted times an index represents the best score as a histogram.
  // Histogram data is written out to file 'histogram_filename'. Returns a
  // vector of plot points for an ROC graph based on wmm's predictions.
  std::vector<ROCPlotPoint> GetPlotData(
    const std::vector<DNASequence>& test_sequence_vec,
    const WeightMatrixModel& wmm,
    const std::string histogram_filename);

  // Returns the base 2 log of x.
  double log2(double x);

 private:
  // Vector of sequences that we are searching for the motif.
  std::shared_ptr<const std::vector<DNASequence>> sequence_vec_;

  // The expected length of the motif we are searching for.
  int motif_length_;

  // Number of potential mutations in the motif we are looking for.
  int num_mutations_;

  // High initial entropy.
  WeightMatrixModel meme_model_;

  // WeightMatrixModel that was initially seeded by random projection scans.
  WeightMatrixModel random_projection_model_;
};

} // namespace

#endif // MOTIF_FINDER_H_
