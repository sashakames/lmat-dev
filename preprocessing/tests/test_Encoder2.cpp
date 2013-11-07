#include <iostream>
#include "Encoder.hpp"
#include <string>
#include <cassert>
using namespace std;
using namespace metag;

string use("various tests for static Encoder methods\n";

int K = 20;

int main(int argc, char **argv) {
  string seq("GGGCTGATATTCTTAAAACCAAATATTTAAAAAATGTAAATATGTTAATA");
  for (size_t j=0; j<seq.size(); j++) {
    seq[j] = tolower(seq[j]);
  }
  string seq_rc;
  Encoder::rc(seq, seq_rc);
  cout << "seq: " << seq << endl;
  cout << "rc:  " << seq << endl << endl;

    Encoder e(seq, K);
    string mer, mer_rc, mer_rc2;
    uint64_t kmer, rc, kmer2;
    while (e.next(kmer)) {
      Encoder::decode(kmer, K, mer); //get kmer as string

      rc = Encoder::rc(kmer, K);  //get rc as uint64_t
      Encoder::decode(rc, K, mer_rc); //get rc as string, from rc as uint64_t

      Encoder::rc(mer, mer_rc2); //get rc as string -- as rc of kmer as string
      cout << "mer:    " << mer << endl;
      cout << "mer rc: " << mer_rc << endl;

      assert(mer_rc == mer_rc2);
      assert(seq.find(mer) != string::npos);
      assert(seq_rc.find(mer_rc) != string::npos);

      //rc of rc
      kmer2 = Encoder::rc(rc, K);  //get rc as uint64_t
      assert(kmer2 == kmer);
    }
}

