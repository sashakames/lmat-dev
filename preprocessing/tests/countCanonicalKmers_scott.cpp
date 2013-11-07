#include <iostream>
#include "all_headers.hpp"
#include <fstream>
#include <cstdlib>
#include <string>
#include <ext/hash_map>


using namespace std;
using namespace metag;

//typedef unsigned long kmer_t;

#define ENCODE(t, c, k) \
switch (c) { \
case 'a': case 'A': t = 0; break; \
case 'c': case 'C': t = 1; break; \
case 'g': case 'G': t = 2; break; \
case 't': case 'T': t = 3; break; \
default: k = 0; continue; \
}

__gnu_cxx::hash_map<uint64_t, char> kmers;

/*
 *    Given a sequence of nucleotide characters,
 *    break it into canonical k-mers in one pass.
 *    Nucleotides are encoded with two bits in
 *    the k-mer. Any k-mers with ambiguous characters
 *    are skipped.
 *    str:  DNA sequence (read)
 *    slen: DNA sequence length in nucleotides
 *    klen: k-mer length in nucleotides
 **/
static
void seq_lookup(const char *str, int slen, int klen) {
  int j; /* position of last nucleotide in sequence */
  int k = 0; /* count of contiguous valid characters */
  int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
  kmer_t mask = ((kmer_t)1 << klen*2)-1; /* bits covering encoded k-mer */
  kmer_t forward = 0; /* forward k-mer */
  kmer_t reverse = 0; /* reverse k-mer */
  kmer_t kmer; /* canonical k-mer */

  for (j = 0; j < slen; j++) {
    register int t;
    ENCODE(t, str[j], k);
    forward = ((forward << 2) | t) & mask;
    reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
    if (++k >= klen) {
      kmer = (forward < reverse) ? forward : reverse;
      kmers[kmer] = 'c';
      /* kmer_lookup(kmer); do k-mer lookup here... */
      /* zero based position of forward k-mer is (j-klen+1) */
      /* zero based position of reverse k-mer is (slen-j-1) */
    }
  }
}


int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " fasta_fn K\n";
    exit(-1);
  }

  int K = atoi(argv[2]);

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
    seq_lookup(seq.c_str(), seq.size(), K); 
    /*
    CanonicalEncoder e(seq, K);
    while (e.next(kmer)) {
      kmers[kmer] = c;
      ++ct;
      if (ct % 1000000 == 0) cout << ct/1000000 << " M" << endl;
    }
    */
  }
  cout << "kmer count: " << kmers.size() << endl;
  in.close();
}
