/*
 *
 *
 *
 */

#include <chrono>
#include <iostream>
#include <limits>
#include <vector>

#include "count_min_sketch.h"
#include "xxhash.h"

using namespace std;

//-----------------------------------------------------------------------------

static int kSeed = std::chrono::system_clock::now().time_since_epoch().count();
static default_random_engine kRandomGenerator;

//-----------------------------------------------------------------------------

static const int kShortMax = numeric_limits<short>::max();

//-----------------------------------------------------------------------------

CountMinSketch::CountMinSketch(int num_rows, int row_size,
                               int num_down_samples) :
  // We start with sampling at 0.5.
  num_down_samples_(num_down_samples) {

  geo_dist_ = geometric_distribution<int>(1.0 / (1 << num_down_samples_));
  num_skips_remaining_ = geo_dist_(kRandomGenerator);

  typedef chrono::time_point<chrono::system_clock> time_point;
  for (int ii = 0; ii < num_rows; ++ii) {
    const time_point tp = chrono::system_clock::now();
    row_seeds_.push_back(chrono::duration_cast<std::chrono::nanoseconds>(
                         tp.time_since_epoch()).count());
    sketch_.push_back(vector<short>(row_size, 0));
  }
}

//-----------------------------------------------------------------------------

void CountMinSketch::Increment(const string& key) {
  if (num_skips_remaining_ > 0) {
    --num_skips_remaining_;
    return;
  }
#ifdef DEBUG
  assert(row_seeds_.size() == sketch_.size());
#endif

  #ifdef CONSERVATIVE_UPDATE
  // We need the RawCount for conservative update
  short current_val = GetRawCount(key);
  if (current_val == kShortMax) {
    current_val >>= 1;
    DownSample();
  }
  #endif

  for (int ii = 0; ii < row_seeds_.size(); ++ii) {
    XXH64_hash_t hash = XXH64((void*)key.c_str(), key.size(), row_seeds_[ii]);
    short *const inc = &sketch_[ii][hash % sketch_[ii].size()];
    #ifdef CONSERVATIVE_UPDATE
      if (*inc == current_val) {
        *inc += 1;
      }
    #else
      ++*inc;
      if (*inc == kShortMax) {
        DownSample();
      }
    #endif
  }

  // Set the num_skips_remaining_ now that we have incremented.
  num_skips_remaining_ = geo_dist_(kRandomGenerator);
}
//-----------------------------------------------------------------------------

long CountMinSketch::GetCount(const string& key) {
  return (long)GetRawCount(key) << num_down_samples_;
}

//-----------------------------------------------------------------------------

short CountMinSketch::GetRawCount(const std::string& key) {
#ifdef DEBUG
  assert(row_seeds_.size() == sketch_.size());
#endif
  short min_count = kShortMax;
  for (int ii = 0; ii < row_seeds_.size(); ++ii) {
    XXH64_hash_t hash = XXH64((void*)key.c_str(), key.size(), row_seeds_[ii]);
    min_count = min(min_count, sketch_[ii][hash % sketch_[ii].size()]);
  }
  return min_count;
}

//-----------------------------------------------------------------------------

void CountMinSketch::DownSample() {
  cout << "Down sampling" << endl;
  for (int ii = 0; ii < sketch_.size(); ++ii) {
    for (int kk = 0; kk < sketch_[0].size(); ++kk) {
      // Divide all values by 2.
      sketch_[ii][kk] >>= 1;
    }
  }

  // Increment the number of times we have down sampled.
  ++num_down_samples_;

  // Set the geometric distribution based on new geo_p_.
  geo_dist_ = geometric_distribution<int>(1.0 / (1 << num_down_samples_));
}

//-----------------------------------------------------------------------------
