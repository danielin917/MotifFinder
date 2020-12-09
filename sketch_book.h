/*
 * Interface for 'SketchBook' a streaming algorithm for finding fixed length
 * motifs in large datasets using random projection.
 *
 * Author: danielin@umich.edu
 *
 */

#ifndef SKETCH_BOOK_H_
#define SKETCH_BOOK_H_

#include <string>
#include <vector>

#include "count_min_sketch.h"
#include "dna_sequence.h"

namespace motif {

class SketchBook {
  public:
    typedef std::vector<std::vector<double>> CountMatrix;

    // Constructor, takes in 'filename' for large file containing sequences
    // that are to be scanned for motifs.
    SketchBook(const std::string& filename,
               int motif_length,
               int num_mutations);

    // State containing all relevant information related to the contained
    // CountMinSketch.
    struct SketchState {
      // Constructor that takes int the projection indices and the backing
      // 'p_count_min' sketch.
      SketchState(std::vector<int>&& p_projection_indices_vec,
                  CountMinSketch&& p_count_min);

      bool operator<(const SketchState& other) const {
        return current_best_min_count < other.current_best_min_count;
      }

      // Projection indices used to get the projection 'hash' that was used to
      // key the CountMinSketch.
      std::vector<int> projection_indices_vec;

      // CountMinSketch backing this random projection.
      CountMinSketch count_min;

      // The full motif word associated with the most recent best key.
      std::deque<std::string> best_motif_word_deque;

      // The count associated with the current best string.
      long current_best_min_count;
    };

    // Scan FASTA file and look for motifs of length 'k_num_'.
    // Returns the consensus string found during scan.
    std::string Scan();

  // Process 'sequence' and look for motifs.
  void ProcessSequence(const std::string& sequence);

  // Read 'line' in from 'buffer' starting at 'index'. If null at end of buffer
  // reached break early. Returns the next index to be read from buffer after
  // reading line.
  int GetLine(const char *buffer, int index, std::string *line);

  private:
    // Filename for FASTA file.
    std::string filename_;

    // Expected motif length.
    int motif_length_;

    // Number of mutations expected on the motif.
    int num_mutations_;

    // The size we will use for out projections.
    int k_num_;

    // Vector containing sketches associated with all of our random
    // projections.
    std::vector<SketchState> sketch_state_vec_;

    // Sequence vector that will be used to perform EM.
    std::shared_ptr<std::vector<DNASequence>> sequence_vec_;
};

} // namespace

#endif // SKETCH_BOOK_H_
