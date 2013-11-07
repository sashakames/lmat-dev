#include "../all_headers.hpp"
#include <iostream>
#include <vector>

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains output from jellylist.bin\n";
    exit(-1);
  }


  FILE *fp = Utils::openReadFile(argv[1]);

  //get the kmer length and number of kmers in the file
  KmerFileMetaData metadata;
  metadata.read(argv[1]);
  uint64_t kmer_count = metadata.size(); 

  KmerNode w;
  for (uint64_t j=0; j<kmer_count; j++) {
      w.read(fp);
      w.write();
  }

  fclose(fp);
}
