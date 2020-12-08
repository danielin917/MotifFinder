/*
 * This file contains the interface for a count min sketch that
 * probabilistically increments counts. Also performs conservative update.
 *
 * Author: danielin@umich.edu
 */

#ifndef COUNT_MIN_SKETCH_H_
#define COUNT_MIN_SKETCH_H_

#include <random>
#include <string>
#include <vector>

class CountMinSketch {
 public:
  CountMinSketch(int num_rows, int row_size, int num_down_samples = 0);

  // Probabilistically increment the sketch based on 'input'.
  void Increment(const std::string& input);

  // Returns the current estimated count for this key.
  int GetCount(const std::string& key);

 private:
  // Down sample all counters and update sampling probabilities.
  void DownSample();

  // Returns the min short value contained in the sketch with no
  // modifications to account for sampling.
  short GetRawCount(const std::string& key);

 private:
  // Sketch matrix. Each row will contain counts that are accessed by pairwise
  // independent hash functions for incrementing.
  std::vector<std::vector<short>> sketch_;

  // Seeds that will be used for the hash funtion of each row.
  std::vector<int64_t> row_seeds_;

  // The geometric_distribution object associated with this sketch.
  std::geometric_distribution<int> geo_dist_;

  // Number of skips to be done before the next increment.
  int num_skips_remaining_;

  // Number of down samples that have occurred.
  int num_down_samples_;
};

#endif
