#include "all_headers.hpp"
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

  KmerFileMetaData metadata;
  metadata.read(fp);
  uint32_t kmer_len = metadata.kmerLength();
  uint64_t kmer_count = metadata.size();
  string mer;


  KmerNode w;
  uint64_t test, sanity = ~0;
  set<uint32_t> gids;
  for (uint64_t j=0; j<kmer_count; j++) {
      w.read(fp);
      decode(w.getKmer(), kmer_len, mer);
      const std::set<uint32_t> &g = w.getGenomes();
      cout << w.getKmer() << " " << mer << " ";
      gids.clear();
      for (set<uint32_t>::const_iterator t = g.begin(); t != g.end(); t++) {
        cout << *t << " ";
      }
      cout << endl;

      if ((1+j) % 1000 == 0) {
        assert(fread(&test, sizeof(uint64_t), 1, fp) == 1);
        assert(test == sanity);
      }
  }

  fclose(fp);
}
