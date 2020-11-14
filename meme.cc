/*
 *
 *
 */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "meme.h"
#include "meme_graph.h"

using namespace std;

namespace meme {

namespace {

//-----------------------------------------------------------------------------

static const vector<double> kDefaultBackgroundFrequencyVec =
  { 0.25, 0.25, 0.25, 0.25 };
static const vector<double> kDefaultPseudoCountVec =
  { 0.25, 0.25, 0.25, 0.25 };
static const int kDefaultSeedIterations = 3;
static const double kSeedPercent = 0.85;
static const double kAverageColumnThreshold = 1;
static const string kHistogramFileNameA = "histogram-A.dat";
static const string kHistogramFileNameB = "histogram-B.dat";
static const string kHistogramFileNameC = "histogram-C.dat";
static const string kHistogramFileNameD = "histogram-D.dat";
static const string kRocFilename = "roc.dat";

//-----------------------------------------------------------------------------

enum class Nucleotide {
  kA = 0 /* Adenine */,
  kC = 1 /* Cytosine */,
  kG = 2 /* Guanine */,
  kT = 3 /* Thymine */,
  kNumNucleotides = 4
};

//-----------------------------------------------------------------------------

static Nucleotide CharToNucleotide(const char c) {
  const char upperc = toupper(c);
  switch (upperc) {
    case 'A':
      return Nucleotide::kA;
    case 'C':
      return Nucleotide::kC;
    case 'G':
      return Nucleotide::kG;
    case 'T':
      return Nucleotide::kT;
    default:
      cerr << "Invalid nucleotide char " << c << endl;
      assert(false);
      break;
  }
  __builtin_unreachable();
}

//-----------------------------------------------------------------------------

} // Anonymous namespace

//-----------------------------------------------------------------------------

void Meme::WeightMatrixModel::PrintWMM() const {
  assert(matrix_model.size() == 4);
  static const vector<char> kNucleotideChars = { 'A', 'C', 'G', 'T' };
  for (int ii = 0; ii < matrix_model.size(); ++ii) {
    cout << kNucleotideChars[ii] << ": ";
    for (int kk = 0; kk < matrix_model[ii].size(); ++kk) {
      cout << matrix_model[ii][kk] << " ";
    }
    cout << endl;
  }
}

//-----------------------------------------------------------------------------

void Meme::WeightMatrixModel::PrintFrequencyMatrix() const {
  assert(frequency_matrix.size() == 4);
  static const vector<char> kNucleotideChars = { 'A', 'C', 'G', 'T' };
  for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
    cout << kNucleotideChars[ii] << ": ";
    for (int kk = 0; kk < frequency_matrix[ii].size(); ++kk) {
      string frequency_string = to_string(frequency_matrix[ii][kk]);
      // Pad the back of short entropy numbers.
      if (frequency_string.size() < 8) {
        frequency_string += string(8 - (int)frequency_string.size(), ' ');
      }
      // Trim the end of larger entropy numbers.
      if (frequency_string.size() > 8) {
        frequency_string = frequency_string.substr(0, 8);
      }
      cout << frequency_string << ", ";
    }
    cout << endl;
  }
}

//-----------------------------------------------------------------------------

Meme::Meme(vector<DNASequence>&& sequence_vec,
           const int k_length) :
  sequence_vec_(move(sequence_vec)),
  k_length_(k_length) {
}

//-----------------------------------------------------------------------------

void Meme::Test(const vector<DNASequence>& test_sequence_vec) {
  assert(wmm_a_.matrix_model.size() > 0);
  assert(wmm_b_.matrix_model.size() > 0);
  assert(wmm_c_.matrix_model.size() > 0);
  assert(wmm_d_.matrix_model.size() > 0);

  // Mapping from index to the number of times that index has the best score
  // within a sequence.
  vector<ROCPlotPoint> roc_plot_vec_a =
    GetPlotData(test_sequence_vec, wmm_a_, kHistogramFileNameA);
  vector<ROCPlotPoint> roc_plot_vec_b =
    GetPlotData(test_sequence_vec, wmm_b_, kHistogramFileNameB);
  vector<ROCPlotPoint> roc_plot_vec_c =
    GetPlotData(test_sequence_vec, wmm_c_, kHistogramFileNameC);
  vector<ROCPlotPoint> roc_plot_vec_d =
    GetPlotData(test_sequence_vec, wmm_d_, kHistogramFileNameD);

  // Make sure all ROC plot point vector sizes are equal.
  assert(roc_plot_vec_a.size() == roc_plot_vec_b.size());
  assert(roc_plot_vec_b.size() == roc_plot_vec_c.size());
  assert(roc_plot_vec_c.size() == roc_plot_vec_d.size());

  cout << "AUC A: " << CalculateAUC(&roc_plot_vec_a) << endl;
  cout << "AUC B: " << CalculateAUC(&roc_plot_vec_b) << endl;
  cout << "AUC C: " << CalculateAUC(&roc_plot_vec_c) << endl;
  cout << "AUC D: " << CalculateAUC(&roc_plot_vec_d) << endl;

  // Now right the ROC data out a .dat file. All roc_plot_vecs should be same
  // size. We will print data for each roc line into columns with 2 columns
  // per line.
  ofstream ofs(kRocFilename.c_str(), ofstream::out);
  for (int ii = 0; ii < roc_plot_vec_a.size(); ++ii) {
    ofs << roc_plot_vec_a[ii].fp_rate        /* Plot A */
        << " " << roc_plot_vec_a[ii].tp_rate
        << " " << roc_plot_vec_b[ii].fp_rate /* Plot B */
        << " " << roc_plot_vec_b[ii].tp_rate
        << " " << roc_plot_vec_c[ii].fp_rate /* Plot C */
        << " " << roc_plot_vec_c[ii].tp_rate
        << " " << roc_plot_vec_d[ii].fp_rate /* Plot D */
        << " " << roc_plot_vec_d[ii].tp_rate;
    if (ii == 0) {
      // Add "Line of no discrimination" start.
      ofs << " " << 0 << " " << 0;
    }

    if (ii == 1) {
      // Add "Line of no discrimination" end.
      ofs << " " << 1 << " " << 1;
    }
    ofs << '\n';
  }
  ofs.close();
}

//-----------------------------------------------------------------------------

vector<ROCPlotPoint> Meme::GetPlotData(
  const vector<DNASequence>& test_sequence_vec,
  const WeightMatrixModel& wmm,
  const string histogram_filename) {

  vector<vector<double>> score_matrix = scanWMM(wmm, test_sequence_vec);
  assert(score_matrix.size() > 0);

  // We will iterate and count the number of times an index is found to have
  // the best WMM score for a sequence.
  unordered_map<int, int> best_index_count_map;
  size_t greatest_sequence_length = 0;
  for (int ii = 0; ii < score_matrix.size(); ++ii) {
    assert(score_matrix[ii].size() > 0);
    double best_score = numeric_limits<double>::min();
    int best_index = -1;
    if (greatest_sequence_length < (score_matrix[ii].size() - 1)) {
      greatest_sequence_length = score_matrix[ii].size() - 1;
    }
    for (int kk = 0; kk < score_matrix[ii].size(); ++kk) {
      if (score_matrix[ii][kk] > best_score) {
        best_score = score_matrix[ii][kk];
        best_index = kk;
      }
    }
    assert(best_index >= 0);
    best_index_count_map[best_index] += 1;
  }

  // Send histogram data to file.
  ofstream ofs(histogram_filename.c_str(), ofstream::out);
  for (int ii = 0; ii < greatest_sequence_length; ++ii) {
    ofs << best_index_count_map[ii] << '\n';
  }
  ofs.close();

  pair<int, int> peak_pair =
    *max_element(best_index_count_map.begin(),
                 best_index_count_map.end(),
                 [](const pair<int, int>& l, const pair<int, int>& r){
                   return l.second < r.second;
                 });

  const int m_index = peak_pair.first;
  cout << histogram_filename << endl;
  cout << "m: " << m_index << endl;
  const double kmer_threshold = kAverageColumnThreshold * k_length_;

  // Get all plot scores along with information on whether score is attached
  // to a real positive.
  vector<PlotScore> plot_score_vec;
  int num_positives = 0;
  for (int ii = 0; ii < score_matrix.size(); ++ii) {
    for (int kk = 0; kk < score_matrix[ii].size(); ++kk) {
      PlotScore ps;
      ps.score = score_matrix[ii][kk];
      ps.is_positive = kk == m_index;
      plot_score_vec.push_back(ps);
      if (kk == m_index) {
        ++num_positives;
      }
    }
  }
  assert(num_positives == score_matrix.size());

  sort(plot_score_vec.begin(), plot_score_vec.end());

  // Now get all ROC plot points.
  vector<ROCPlotPoint> roc_plot_vec;

  // The initial threshold of zero will result in
  // tpr = num_positives / num_positives and
  // fpr = num_negatives / num_negatives resulting in ROC plot point <1, 1>.
  roc_plot_vec.push_back(ROCPlotPoint(1, 1));

  // At the highest threshold we should not get true or false positives.
  roc_plot_vec.push_back(ROCPlotPoint(0, 0));

  const int num_negatives = plot_score_vec.size() - num_positives;
  double true_positives = num_positives;
  double false_positives = num_negatives;

  // These are specific stats we are recording the largest threshold that
  // captures all true positives.
  double largest_true_positive_tau = -1;
  double largest_false_positives = -1;
  double largest_true_negatives = -1;
  double largest_false_negatives = -1;
  double largest_tpr = -1;
  double largest_fpr = -1;
  for (int ii = 0; ii < plot_score_vec.size(); ++ii) {
    const PlotScore& ps = plot_score_vec[ii];
    if (ps.is_positive) {
      if (true_positives == num_positives) {
        largest_false_positives = false_positives;
        // There should be no false negatives as this will be the first
        // positive sample that will become a false negative.
        largest_true_negatives =
          plot_score_vec.size() - false_positives - true_positives;
        largest_false_negatives = 0;
        largest_fpr = false_positives / num_negatives;
        largest_tpr = true_positives / num_positives;

        largest_true_positive_tau = ps.score;
      }
      true_positives -= 1;
      assert(true_positives >= 0);
    } else {
      false_positives -= 1;
      assert(false_positives >= 0);
    }
    roc_plot_vec.push_back(
      ROCPlotPoint(false_positives / num_negatives /* fp_rate */,
                   true_positives / num_positives /* tp_rate */));
  }

  cout << "Largest True Positive Tau: " << largest_true_positive_tau << endl;
  cout << "Corresponding stats for the above threshold" << endl;
  cout << "True Negatives: " << largest_true_negatives << endl;
  cout << "False Negatives: " << largest_false_negatives << endl;
  cout << "TPR: " << largest_tpr << endl;
  cout << "FPR: " << largest_fpr << endl;
  cout << endl;
  return roc_plot_vec;
}

//-----------------------------------------------------------------------------

Meme::WeightMatrixModel Meme::Run(const int num_final_iterations) {
  assert(sequence_vec_.size() > 0);
  vector<vector<double>> count_matrix;
  for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides); ++ii) {
    count_matrix.push_back(vector<double>(k_length_, 0));
  }

  vector<string> kmer_seeds;
  const string& seed_sequence = sequence_vec_[0].sequence();
  for (int ii = 0; (ii + k_length_) <= seed_sequence.size();
         ii += k_length_/ 2 ) {
    kmer_seeds.push_back(seed_sequence.substr(ii, k_length_));
  }

  // Special case and add the bookended kmer if it has not already been added.
  if (seed_sequence.size() % (k_length_ / 2) != 0) {
    kmer_seeds.push_back(
      seed_sequence.substr(seed_sequence.size() - k_length_, k_length_));
  }

  // Create the initial seed WMMs.
  vector<WeightMatrixModel> wmm_seeds;
  for (int ii = 0; ii < kmer_seeds.size(); ++ii) {
    wmm_seeds.push_back(CreateSeedMatrixModel(kmer_seeds[ii]));
  }

  // Now update each WMM after 3 iterations.
  for (int ii = 0; ii < wmm_seeds.size(); ++ii) {
    wmm_seeds[ii] = RunEM(move(wmm_seeds[ii]), kDefaultSeedIterations);
  }

  // Now sort the resulting WMM based on relative entropy.
  sort(wmm_seeds.begin(), wmm_seeds.end());

  wmm_a_ = wmm_seeds.back();
  wmm_b_ = wmm_seeds[wmm_seeds.size() / 2];
  wmm_c_ = wmm_seeds[0];
  for (int ii = 0; ii < wmm_seeds.size(); ++ii) {
    wmm_seeds[ii] = RunEM(move(wmm_seeds[ii]), num_final_iterations);
  }

  cout << "Seed Entropy Table" << endl;
  for (int ii = 0; ii < wmm_seeds.size(); ++ii) {
    cout << "Seed-" << ii << ": ";
    const vector<double>& entropy_history_vec =
      wmm_seeds[ii].entropy_history_vec;

    // Pad for different sized seed names.
    string seed_name = to_string(ii);
    for (int ii = 0; ii < 4 - seed_name.size(); ++ii) {
      cout << " ";
    }
    for (int kk = 0; kk < entropy_history_vec.size(); ++kk) {
      string entropy_string = to_string(entropy_history_vec[kk]);
      // Pad the back of short entropy numbers.
      if (entropy_string.size() < 8) {
        entropy_string += string(8 - (int)entropy_string.size(), ' ');
      }
      // Trim the end of larger entropy numbers.
      if (entropy_string.size() > 8) {
        entropy_string = entropy_string.substr(0, 8);
      }
      cout << entropy_string << " ";
    }
    cout << endl;
  }
  wmm_d_ = *max_element(wmm_seeds.begin(), wmm_seeds.end());

  cout << endl;
  cout << "Frequency Matrix A " << endl;
  wmm_a_.PrintFrequencyMatrix();
  wmm_a_.PrintWMM();

  cout << "Frequency Matrix B " << endl;
  wmm_b_.PrintFrequencyMatrix();

  cout << "Frequency Matrix C " << endl;
  wmm_c_.PrintFrequencyMatrix();

  cout << "Frequency Matrix D " << endl;
  wmm_d_.PrintFrequencyMatrix();
  cout << endl;

  return wmm_d_;
}

//-----------------------------------------------------------------------------

Meme::WeightMatrixModel Meme::RunEM(WeightMatrixModel&& wmm,
                                    const int num_iterations) {

  WeightMatrixModel next_wmm = move(wmm);
  for (int ii = 0; ii < num_iterations; ++ii) {
    const vector<vector<double>> yhat_matrix = Estep(next_wmm);
    Mstep(yhat_matrix, kDefaultPseudoCountVec,
          kDefaultBackgroundFrequencyVec, &next_wmm);
    assert(next_wmm.matrix_model.size() > 0);
  }
  return next_wmm;
}

//-----------------------------------------------------------------------------

Meme::WeightMatrixModel Meme::CreateSeedMatrixModel(const string& seed_kmer) {
  assert(seed_kmer.size() == k_length_);
  vector<vector<double>> count_matrix;
  for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides); ++ii) {
    count_matrix.push_back(vector<double>(k_length_, 0));
  }
  updateCountMatrix(seed_kmer, {} /* yhat_vec */, &count_matrix);

  for (int kk = 0; kk < k_length_; ++kk) {
    // Sanity check that our count matrix is made with just one kmer char count
    // in each column.
    double sum = 0;
    for (int ii = 0; ii < count_matrix.size(); ++ii) {
      sum += count_matrix[ii][kk];
    }
    assert(sum == 1);
  }

  // We want to add pseudocounts so that the original count matrix contributes
  // only kSeedPercent to the count in each column. Since there should only be
  // a sum of "1" count in each column we can get the total final count as
  // Total = 1 / kSeedPercent. The difference between total and the seed count
  // pseudo_count_sum = Total - 1 is the sum of all pseudocounts. We will then
  // distribute this pseudocount evenly across all nucleotide counts.
  const double kFinalColumnSum = 1 / kSeedPercent;
  const double kPseudoCount = (kFinalColumnSum - 1) / 4;
  vector<double> pseudo_count_vec(
    static_cast<int>(Nucleotide::kNumNucleotides),
    kPseudoCount);

  // Add the calculated pseudocounts.
  addPseudo(pseudo_count_vec, &count_matrix);

  vector<vector<double>> frequency_matrix = makeFrequencyMatrix(count_matrix);
  return makeWMM(frequency_matrix, kDefaultBackgroundFrequencyVec);
}

//-----------------------------------------------------------------------------

vector<vector<double>> Meme::makeCountMatrix(
  const vector<vector<double>>& yhat_matrix) {
  assert(!sequence_vec_.empty());

  vector<vector<double>> count_matrix;
  for (int ii = 0; ii < static_cast<int>(Nucleotide::kNumNucleotides); ++ii) {
    count_matrix.push_back(vector<double>(k_length_, 0));
  }

  for (int ii = 0; ii < sequence_vec_.size(); ++ii) {
    const string& sequence = sequence_vec_[ii].sequence();
    updateCountMatrix(sequence, yhat_matrix[ii], &count_matrix);
  }
  return count_matrix;
}

//-----------------------------------------------------------------------------

void Meme::updateCountMatrix(const string sequence_str,
                             const vector<double>& yhat_vec,
                             vector<vector<double>> *const count_matrix) {
  assert(count_matrix);
  assert(count_matrix->size() ==
           static_cast<int>(Nucleotide::kNumNucleotides));

  assert(k_length_ > 0);
  if (!yhat_vec.empty()) {
    assert(yhat_vec.size() == sequence_str.size() - k_length_ + 1);
  }

  // Iterate over starting indexes.
  for (int ii = 0; ii < (sequence_str.size() - k_length_ + 1); ++ii) {
    // Iterate over the corresponding k-mer columns given the starting index.
    for (int kk = 0; kk < k_length_; ++kk) {
      const Nucleotide n = CharToNucleotide(sequence_str[ii + kk]);
      assert(static_cast<int>(n) < count_matrix->size());
      // Add the count proportionally to how likely we are currently using the
      // starting index.
      (*count_matrix)[static_cast<int>(n)][kk] +=
        yhat_vec.empty() ? 1 : yhat_vec[ii];
    }
  }
}

//-----------------------------------------------------------------------------

void Meme::addPseudo(const vector<double>& pseudo_count_vec,
                     vector<vector<double>> *const count_matrix) {
  assert(count_matrix);
  assert(pseudo_count_vec.size() ==
           static_cast<int>(Nucleotide::kNumNucleotides));
  for (int ii = 0; ii < count_matrix->size(); ++ii) {
    for (int kk = 0; kk < (*count_matrix)[ii].size(); ++kk) {
      (*count_matrix)[ii][kk] += pseudo_count_vec[ii];
    }
  }
}

//-----------------------------------------------------------------------------

vector<vector<double>> Meme::makeFrequencyMatrix(
  const vector<vector<double>>& count_matrix) {
  assert(count_matrix.size() == static_cast<int>(Nucleotide::kNumNucleotides));
  assert(count_matrix[0].size() == k_length_);

  // Fill in the new frequency_matrix with the same dimensions as the count
  // matrix.
  vector<vector<double>> frequency_matrix;
  for (int ii = 0; ii < count_matrix.size(); ++ii) {
    frequency_matrix.push_back(vector<double>(count_matrix[ii].size(), 0));
  }

  for (int kk = 0; kk < k_length_; ++kk) {
    double total_count = 0;
    for (int ii = 0; ii < count_matrix.size(); ++ii) {
      total_count += count_matrix[ii][kk];
    }
    assert(total_count > 0);

    for (int ii = 0; ii < count_matrix.size(); ++ii) {
      frequency_matrix[ii][kk] = count_matrix[ii][kk] / total_count;
    }
  }
  return frequency_matrix;
}

//-----------------------------------------------------------------------------

double Meme::entropy(const vector<vector<double>>& frequency_matrix,
                     const vector<double>& background_frequency_vec) {
  // Iterate over each column summing up relative entropy.
  double entropy = 0;
  for (int kk = 0; kk < frequency_matrix[0].size(); ++kk) {
    for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
      // For each column we are adding entropy based on the model probability
      // distribution vs the background probability distribution over the same
      // alphabet.
      // Relative entropy is given by H(P||Q):
      //  for all i in alphabet A SUM { P(xi) * log2(P(xi) / Q(xi)) }
      assert(frequency_matrix[ii][kk] <= 1);
      assert(background_frequency_vec[ii] > 0);
      entropy +=
        frequency_matrix[ii][kk] *
          log2(frequency_matrix[ii][kk] / background_frequency_vec[ii]);
    }
  }
  return entropy;
}

//-----------------------------------------------------------------------------

void Meme::updateMatrixModel(
  const vector<vector<double>>& frequency_matrix,
  const vector<double>& background_frequency_vec,
  WeightMatrixModel *const wmm) {

  WeightMatrixModel new_wmm =
    makeWMM(frequency_matrix, background_frequency_vec);
  wmm->matrix_model = move(new_wmm.matrix_model);
  wmm->frequency_matrix = move(new_wmm.frequency_matrix);
  wmm->entropy_history_vec.push_back(new_wmm.entropy_history_vec.back());
}

//-----------------------------------------------------------------------------

Meme::WeightMatrixModel Meme::makeWMM(
  const vector<vector<double>>& frequency_matrix,
  const vector<double>& background_frequency_vec) {

  WeightMatrixModel wmm;
  vector<vector<double>>& weight_matrix = wmm.matrix_model;
  for (int ii = 0; ii < frequency_matrix.size(); ++ii) {
    weight_matrix.push_back(vector<double>(frequency_matrix[0].size(), 0));
    for (int kk = 0; kk < frequency_matrix[ii].size(); ++kk) {
      assert(background_frequency_vec[ii] > 0);
      weight_matrix[ii][kk] =
        log2(frequency_matrix[ii][kk] / background_frequency_vec[ii]);
    }
  }
  wmm.frequency_matrix = frequency_matrix;

  // TODO(DAN): Maybe do the entropy calculation inline above.
  wmm.entropy_history_vec.push_back(
    entropy(frequency_matrix, background_frequency_vec));
  return wmm;
}

//-----------------------------------------------------------------------------

vector<vector<double>> Meme::scanWMM(const WeightMatrixModel& wmm,
                                     const vector<DNASequence>& sequence_vec) {
  assert(wmm.matrix_model.size() > 0);
  vector<vector<double>> sequence_scores;
  const int klength = wmm.matrix_model[0].size();
  for (int ii = 0; ii < sequence_vec.size(); ++ii) {
    // We insert empty scores for all possible k-windows for a weight matrix of
    // width k.
    sequence_scores.push_back(
      vector<double>(sequence_vec[ii].sequence().size() - klength + 1, 0));
  }

  // Loop all sequences.
  for (int ii = 0; ii < sequence_vec.size(); ++ii) {
    const string& sequence = sequence_vec[ii].sequence();
    // Loop all starting indexes.
    for (int kk = 0; kk < sequence.size() - klength + 1; ++kk) {
      // Now scan this individual k-mer at this offset and score it.
      double score = 0;
      for (int xx = 0; xx < klength; ++xx) {
        const Nucleotide n = CharToNucleotide(sequence[kk + xx]);
        score += wmm.matrix_model[static_cast<int>(n)][xx];
      }
      sequence_scores[ii][kk] = score;
    }
  }
  return sequence_scores;
}

//-----------------------------------------------------------------------------

vector<vector<double>> Meme::Estep(const WeightMatrixModel& wmm) {
  vector<vector<double>> score_matrix = scanWMM(wmm, sequence_vec_);

  // The total vec maintains a total sum over all ~probabilities for each index
  // in the corresponding sequence.
  vector<long double> score_total_vec(score_matrix.size(), 0);

  // Yij = Product of the probabilities of each found char given model theta
  //       normalized by c so that the sume of Yij over all j for a given i
  //       is 1.
  // Our scores are the sum of the logs of (f/b). Which is equivalent to the
  // log of the product of (f/b). Therefore 2^score = Product(f/b) using f as
  // frequency and b as background. Since the background is a constant we can
  // say that for all scores the probability of the k-mer given the model is
  // proportional to 2^score.
  for (int ii = 0; ii < score_matrix.size(); ++ii) {
    // Iterate through the scores for each sequence and get the total sum of
    // all scores for the sequence. We do this to use this as the denominator
    // so the sum of all yhat probabilities is 1.
    for (int kk = 0; kk < score_matrix[ii].size(); ++kk) {
      score_total_vec[ii] += pow(2, score_matrix[ii][kk]);
    }
  }

  vector<vector<double>> yhat_matrix;
  for (int ii = 0; ii < score_matrix.size(); ++ii) {
    // We will only fill yhat with viable indices at least k length from the
    // end of a sequence.
    yhat_matrix.push_back(vector<double>(score_matrix[ii].size(), 0));
    for (int kk = 0; kk < score_matrix[ii].size(); ++kk) {
      assert(score_total_vec[ii] > 0);
      yhat_matrix[ii][kk] =
        pow(2, score_matrix[ii][kk]) / score_total_vec[ii];
    }
  }
  return yhat_matrix;
}

//-----------------------------------------------------------------------------

void Meme::Mstep(const vector<vector<double>>& yhat_matrix,
                 const vector<double>& pseudo_count_vec,
                 const vector<double>& background_frequency_vec,
                 WeightMatrixModel *const wmm) {

  vector<vector<double>> count_matrix = makeCountMatrix(yhat_matrix);
  addPseudo(pseudo_count_vec, &count_matrix);
  vector<vector<double>> frequency_matrix = makeFrequencyMatrix(count_matrix);
  updateMatrixModel(frequency_matrix, background_frequency_vec, wmm);
}

//-----------------------------------------------------------------------------

double Meme::log2(const double x) {
  return log(x) / log(2);
}

//-----------------------------------------------------------------------------

}
