/*
 * Simple C++ timer
 *
 * Author: danielin@umich.edu
 */

#include <chrono>

class Timer {
 public:
  enum class TimeUnit {
    kSecs,
    kMsecs,
    kUsecs
  };

  void Start() {
    t_ = std::chrono::high_resolution_clock::now();
  }

  // Stop the timer and return the time since start in 'unit' time units.
  int64_t Stop(const TimeUnit unit) {
    auto delta = std::chrono::high_resolution_clock::now() - t_;
    switch (unit) {
      case TimeUnit::kSecs:
        return
          std::chrono::duration_cast<std::chrono::seconds>(delta).count();
      case TimeUnit::kMsecs:
        return
          std::chrono::duration_cast<std::chrono::milliseconds>(delta).count();
      case TimeUnit::kUsecs:
        return
          std::chrono::duration_cast<std::chrono::microseconds>(delta).count();
    }
    throw -1;
  }
 private:
  std::chrono::high_resolution_clock::time_point t_;
};
