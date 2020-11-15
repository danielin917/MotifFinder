/*
 * This file contains the driver for the MEME implementation.
 *
 * Author: Dan Lin, danielin@uw.edu
 */

#include <fstream>
#include <iostream>
#include <vector>

#include "em_tester.h"
#include "random_projection.h"

using namespace motif;
using namespace std;

static constexpr int kDefaultNumIterations = 1000;

bool ExtractSequencesFromFile(const string& filename,
                              vector<DNASequence> *dna_sequences);

int main(int argc, char *argv[]) {
  if (argc < 4) {
    cout << "Not enough args provided." << endl
         << "Format: ./a.out <motif_length> <training_file> <test_file>"
         << "<optional: num_iterations>" << endl;
    return -1;
  }
  const int k_length = atoi(argv[1]);

  // First extract the sequences from the training set.
  shared_ptr<vector<DNASequence>> train_sequences =
    make_shared<vector<DNASequence>>();
  string train_filename(argv[2]);
  if (!ExtractSequencesFromFile(train_filename, train_sequences.get())) {
    cerr << "Error reading file " << train_filename << endl;
    return -1;
  }

  // Now extract the sequences from the test set.
  string test_filename(argv[3]);
  vector<DNASequence> test_sequences;
  if (!ExtractSequencesFromFile(test_filename, &test_sequences)) {
    cerr << "Error reading file " << test_filename << endl;
    return -1;
  }

  int num_iterations= kDefaultNumIterations;
  if (argc > 4) {
    num_iterations = atoi(argv[5]);
  }

  int num_mutations = 3;
  RandomProjection rprojection(train_sequences, k_length,
                               k_length - num_mutations);
  rprojection.RunRandomProjectionScans(1 /* num_iterations */);

  cout << "Running " << num_iterations << " iterations on training set"
       << endl << endl;

  EMTester meme(train_sequences, k_length);
  meme.Run(num_iterations);
  meme.Test(test_sequences);
  return 0;
}

bool ExtractSequencesFromFile(const string& filename,
                              vector<DNASequence> *const dna_sequences) {

  ifstream ifs;
  ifs.open(filename.c_str(), std::ifstream::in);
  if (!ifs.good()) {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }

  string sequence;
  string identifier;
  while (true) {
    char peek_char = ifs.peek();
    if (peek_char == EOF) {
      // All input has been read.
      break;
    }

    string line_str;
    getline(ifs, line_str);
    if (peek_char == '>') {
      if (!sequence.empty()) {
        dna_sequences->push_back(DNASequence(sequence, identifier));
        sequence = "";
        identifier = "";
      }
      // Get the identifier set after '>' and prior to ':'.
      int end = line_str.find(":");
      identifier = line_str.substr(1, end - 1);
      continue;
    }

    line_str = DNASequence::MakeValidDNASequence(line_str);
    sequence +=line_str;
  }
  return true;
}
