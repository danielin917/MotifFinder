/*
 *
 *
 *
 */

#include <chrono>
#include <limits>
#include <vector>

#include "count_min.h"
#include "xxhash.h"

using namespace std;

//-----------------------------------------------------------------------------

CountMin::CountMin(int num_rows, int row_size) {
  typedef chrono::time_point<chrono::system_clock> time_point;
  for (int ii = 0; ii < num_rows; ++ii) {
    const time_point tp = chrono::system_clock::now();
    row_seeds_.push_back(chrono::duration_cast<std::chrono::nanoseconds>(
                         tp.time_since_epoch()).count());
    sketch_.push_back(vector<short>(row_size, 0));
  }
}

//-----------------------------------------------------------------------------

void CountMin::Increment(const string& key) {
#ifdef DEBUG
  assert(row_seeds_.size() == sketch_.size());
#endif
  // TODO(Dan): Implement conservative update.
  for (int ii = 0; ii < row_seeds_.size(); ++ii) {
    XXH64_hash_t hash = XXH64((void*)key.c_str(), key.size(), row_seeds_[ii]);
    sketch_[ii][hash % sketch_[ii].size()] += 1;
  }
}
//-----------------------------------------------------------------------------

short CountMin::GetCount(const string& key) {
#ifdef DEBUG
  assert(row_seeds_.size() == sketch_.size());
#endif
  short min_count = numeric_limits<short>::max();
  for (int ii = 0; ii < row_seeds_.size(); ++ii) {
    XXH64_hash_t hash = XXH64((void*)key.c_str(), key.size(), row_seeds_[ii]);
    min_count = min(min_count, sketch_[ii][hash % sketch_[ii].size()]);
  }
  return min_count;
}

//-----------------------------------------------------------------------------
