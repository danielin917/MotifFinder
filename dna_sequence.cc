/*
 *
 */

#include <iostream>

#include "dna_sequence.h"

using namespace std;

namespace motif {

//-----------------------------------------------------------------------------

DNASequence::DNASequence() {

}

//-----------------------------------------------------------------------------

DNASequence::DNASequence(const string& sequence, const string& identifier) :
  sequence_(sequence),
  identifier_(identifier) {
}

//-----------------------------------------------------------------------------

string DNASequence::MakeValidDNASequence(const string& line) {
  string sequence = line;
  for (int ii = 0; ii < sequence.size(); ++ii) {
    if (!IsValidNucleotide(sequence[ii])) {
      // Hack just make all invalid chars a T.
      sequence[ii] = 'T';
    }
  }
  return sequence;
}

//-----------------------------------------------------------------------------

bool DNASequence::IsValidDNASequence(const string& sequence) {
  for (int ii = 0; ii < sequence.size(); ++ii) {
    if (!IsValidNucleotide(sequence[ii])) {
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

bool DNASequence::IsValidNucleotide(const char c) {
  const char upper_c = toupper(c);
  if (upper_c != 'A' &&
      upper_c != 'C' &&
      upper_c != 'G' &&
      upper_c != 'T') {
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

void DNASequence::PrintWithColumns() const {
  cout << identifier_ << endl;
  for (int ii = 0; ii < sequence_.size(); ++ii) {
    if (ii % 60 == 0 && ii != 0) {
      cout << endl;
    }
    cout << sequence_[ii];
  }
  cout << endl;
}

//-----------------------------------------------------------------------------

}
