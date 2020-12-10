/*
 * This file contains the driver for the MEME implementation.
 *
 * Author: Dan Lin, danielin@uw.edu
 */

#include <fstream>
#include <iostream>
#include <vector>

#include "motif_finder.h"
#include "random_projection.h"
#include "sketch_book.h"

using namespace motif;
using namespace std;

static constexpr int kDefaultNumIterations = 1000;

bool ExtractSequencesFromFile(const string& filename,
                              vector<DNASequence> *dna_sequences);

int main(int argc, char *argv[]) {
  srand(time(NULL));

  if (argc < 4) {
    cout << "Not enough args provided." << endl
         << "Format: ./a.out <motif_length> <test_file> <algorithm>"
         << "<optional: num_mutations>" << endl;
    return -1;
  }
  const int motif_length = atoi(argv[1]);
  string test_filename(argv[2]);
  int num_mutations = 0;
  if (argc > 4) {
    num_mutations = atoi(argv[4]);
  }
  const string algorithm(argv[3]);
  if (algorithm == "sketch_book") {
    SketchBook sk(test_filename, motif_length, num_mutations);
    cout << "Early SketchBook size:" << sk.GetSizeBytes() << endl;
    string consensus = sk.Scan();
    cout << "Final SketchBook size:" << sk.GetSizeBytes() << endl;
    cout << "Consensus: "<< consensus << endl;
    return 0;
  }


  // Now extract the sequences from the test set.
  shared_ptr<vector<DNASequence>> test_sequences =
    make_shared<vector<DNASequence>>();
  if (!ExtractSequencesFromFile(test_filename, test_sequences.get())) {
    cerr << "Error reading file " << test_filename << endl;
    return -1;
  }


  cout << num_mutations << " mutations expected in motif occurrences"
       << endl << endl;

  MotifFinder finder(test_sequences, motif_length, num_mutations);
  finder.Run(*test_sequences);
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
    if (!ifs.good()) {
      cerr << "Unable to read file: " << filename << endl;
      return false;
    }

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
