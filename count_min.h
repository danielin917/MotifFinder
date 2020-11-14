/*
 * This file contains the interface for a count min class.
 *
 * Author: danielin@umich.edu
 */

#include <string>
#include <vector>

class CountMin {
 public:
  CountMin(int num_rows, int row_size);

  // Increment the sketch based on 'input'.
  void Increment(const std::string& input);

  // Returns the current estimated count for this key.
  short GetCount(const std::string& key);

 private:
  // Sketch matrix. Each row will contain counts that are accessed by pairwise
  // independent hash functions for incrementing.
  std::vector<std::vector<short>> sketch_;

  // Seeds that will be used for the hash funtion of each row.
  std::vector<int64_t> row_seeds_;
};
