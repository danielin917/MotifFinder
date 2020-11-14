/*
 *
 *
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "meme_graph.h"

using namespace std;

namespace meme {

//-----------------------------------------------------------------------------

double CalculateAUC(vector<ROCPlotPoint> *const roc_plot_vec) {
  assert(roc_plot_vec);
  assert(!roc_plot_vec->empty());

  // First sort points based on x axis.
  sort(roc_plot_vec->begin(), roc_plot_vec->end());

  // We assert that the provided data starts with a fpr=0 and tpr=0 and ends
  // with 1, 1.
  assert(roc_plot_vec->at(0).fp_rate == 0);
  assert(roc_plot_vec->at(0).tp_rate == 0);
  assert(roc_plot_vec->back().fp_rate == 1);
  assert(roc_plot_vec->back().tp_rate == 1);

  // Initialize information for the first "left side" of the trapezoid.
  double height_sum = 0;
  double left_fpr = roc_plot_vec->at(0).fp_rate;
  int fpr_count = 0;
  int idx = 0;
  while ((*roc_plot_vec)[idx].fp_rate == left_fpr) {
    height_sum += (*roc_plot_vec)[idx].tp_rate;
    ++fpr_count;
    ++idx;
  }
  // Get average height.
  double left_height = height_sum / fpr_count;

  // We calculate using trapezoidal rule and find the average height at each
  // unique false positive rate.
  double auc_sum = 0;
  double right_fpr = -1;
  double right_height = -1;
  while (idx < roc_plot_vec->size()) {
    right_fpr = (*roc_plot_vec)[idx].fp_rate;
    fpr_count = 0;
    height_sum = 0;

    // Sum all true positive rates found with the same fpr.
    while(idx < roc_plot_vec->size() &&
          (*roc_plot_vec)[idx].fp_rate == right_fpr) {
      height_sum += (*roc_plot_vec)[idx].tp_rate;
      ++fpr_count;
      ++idx;
    }
    // Get average height.
    right_height = height_sum / fpr_count;

    // Now use trapezoidal rule to add to AUC sum.
    double avg_height = (left_height + right_height) / 2;
    const double trapezoid_area = (right_fpr - left_fpr) * avg_height;
    assert(trapezoid_area >= 0);
    auc_sum += trapezoid_area;

    //cout << "Adding " << trapezoid_area << " " << left_fpr << ", "
    //     << right_fpr << " " << avg_height << endl;

    // Now shift right stats over to the left.
    left_height = right_height;
    left_fpr = right_fpr;
  }
  return auc_sum;
}

//-----------------------------------------------------------------------------

};
