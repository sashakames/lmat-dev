#include "all_headers.hpp"
#include "metag_typedefs.hpp"
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;
using namespace metag;

void decode(uint64_t j, int kmer_len, string &mer) {
  static char *x = 0;
  if (!x) x = new char[kmer_len+1];
  static uint64_t mask = 3, tmp;
  static char c;
  x[kmer_len] = '\0';
  for (int k=0; k<kmer_len; k++) {
    tmp = j & mask;
    switch (tmp) {
      case 0 : c = 'a'; break;
      case 1 : c = 'c'; break;
      case 2 : c = 'g'; break;
      case 3 : c = 't'; break;
      default : cerr << "ERROR!!!!!!!!!\n"; exit(1);
    }
    x[kmer_len-k-1] = c;
    
    j = j >> 2;
  }
  mer = x;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains output from jellylist.bin\n";
    exit(-1);
  }

  //this will exit with a msg if the file can't be opened
  FILE *fp = Utils::openReadFile(argv[1]);

  //you MUST call getMetaData; after the call returns, the fp
  //will be placed at the beginning of the data for the first kmer.
  //get the kmer length and number of kmers in the file
  //get the kmer length and number of kmers in the file 
  KmerFileMetaData metadata;
  metadata.read(fp);
  metadata.write();
  uint32_t kmer_len = metadata.kmerLength();
  uint64_t kmer_count = metadata.size();
  cout << "kmer count: " << kmer_count << endl;
  string mer;

  KmerNode w;
  uint64_t test, sanity = ~0;
  for (uint64_t j=0; j<kmer_count; j++) {
      w.read(fp);
      decode(w.getKmer(), kmer_len, mer);
      const set<tid_t> & tids = w.getTaxIDs();
      cout << "j: " << j << " " << w.getKmer() << " " << mer << " tid count: " << tids.size() << " :: ";
      for (set<tid_t>::const_iterator t = tids.begin(); t != tids.end(); t++) {
        cout << *t << " ";
      }
      cout << endl;

      if ((1+j) % KMER_SANITY_COUNT == 0) {
        assert(fread(&test, sizeof(uint64_t), 1, fp) == 1);
        assert(test == sanity);
      }
  }

  fclose(fp);
}
