/*
 * This file contains utility datastructures and functions for graphing data
 * from MEME testing.
 *
 * Author: Dan Lin, danielin@uw.edu
 *
 */

#ifndef MEME_GRAPH_H_
#define MEME_GRAPH_H_

#include <vector>

namespace meme {

//-----------------------------------------------------------------------------

// A struct containing data associated with an ROC plot point in an ROC plot.
struct PlotScore {
  bool operator<(const PlotScore& other) const {
    return score < other.score;
  }

  // Whether this score represents a true positive hit on a motif.
  bool is_positive;

  // Associated score with this plot point.
  double score;
};


//-----------------------------------------------------------------------------

// Plot point for an ROC graph. Represents prediction on sample at a specific
// threshold.
struct ROCPlotPoint {

  ROCPlotPoint(const double fpr, const double tpr) :
    fp_rate(fpr),
    tp_rate(tpr) {

  }

  bool operator<(const ROCPlotPoint& other) const {
    return fp_rate < other.fp_rate;
  }
  // False positive rate.
  double fp_rate;

  // True positive rate.
  double tp_rate;
};

//-----------------------------------------------------------------------------

// Calculates an estaimate of the area under the curve of the ROC plot given
// by 'roc_plot_vec' using "trapezoidal rule". roc_plot_vec will be sorted as
// part of this routine.
double CalculateAUC(std::vector<ROCPlotPoint> *roc_plot_vec);

//-----------------------------------------------------------------------------

} // namespace

#endif  // MEME_GRAPH_H_
