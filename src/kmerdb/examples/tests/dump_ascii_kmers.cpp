#include "../../all_headers.hpp"
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains output from jellylist.bin\n";
    exit(-1);
  }

  //get the kmer length and number of kmers in the file
  KmerFileMetaData metadata;
  FILE *fp = Utils::openReadFile(argv[1]);
  metadata.read(fp);
  metadata.write();
  uint32_t K = metadata.kmerLength();
  uint64_t kmer_count = metadata.size();;

  uint64_t test, sanity = ~0;

  KmerNode w;
  for (uint64_t j=0; j<kmer_count; j++) {
      w.read(fp);
      w.printAscii(cout, K);
      if ((j+1) % 1000 == 0) {
        assert(fread(&test, sizeof(uint64_t), 1, fp) == 1);
        assert(test == sanity);
      }
  }

  fclose(fp);
}
