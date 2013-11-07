#include <iostream>
#include "all_headers.hpp"
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;
using namespace metag;

int main(int argc, char **argv) {
  if (argc != 4) {
    cerr << "usage: " << argv[0] << " fasta_fn K output_fn\n";
    exit(-1);
  }

  int K = atoi(argv[2]);

  ifstream in(argv[1]);
  assert(in);
  string seq;
  getline(in, seq);
  in.close();
  CanonicalEncoder e(seq, K);

  ofstream out(argv[3]);
  assert(out);
  uint64_t kmer;
  string mer;
  while (e.next(kmer)) {
    Encoder::decode(kmer, K, mer);
    out << mer << endl;
  }

  out.close();
}
