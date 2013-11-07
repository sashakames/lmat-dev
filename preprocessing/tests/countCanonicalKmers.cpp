#include <iostream>
#include "all_headers.hpp"
#include <fstream>
#include <cstdlib>
#include <string>
#include <ext/hash_map>


using namespace std;
using namespace metag;

int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " fasta_fn K\n";
    exit(-1);
  }

  int K = atoi(argv[2]);
  __gnu_cxx::hash_map<uint64_t, char> kmers;

  //open fasta file
  ifstream in(argv[1]);
  assert(in);

  string seq;
  char c = 'x';
  uint64_t kmer;
  uint64_t ct = 0;
  while (true) {
    seq = "";
    getline(in, seq);
    if (! seq.size()) break;
    assert(seq[0] == '>');
    getline(in, seq);
    CanonicalEncoder e(seq, K);
    while (e.next(kmer)) {
      kmers[kmer] = c;
      ++ct;
      if (ct % 1000000 == 0) cout << ct/1000000 << " M" << endl;
    }
  }
  cout << "kmer count: " << kmers.size() << endl;
  in.close();
}
