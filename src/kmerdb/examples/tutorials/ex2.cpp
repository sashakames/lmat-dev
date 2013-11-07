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

  //this will exit with a msg if the file can't be opened
  FILE *fp = Utils::openReadFile(argv[1]);

  //you MUST call getMetaData; after the call returns, the fp
  //will be placed at the beginning of the data for the first kmer.
  //get the kmer length and number of kmers in the file
  //get the kmer length and number of kmers in the file 
  KmerFileMetaData metadata;
  metadata.read(fp);
  uint32_t kmer_len = metadata.kmerLength();
  uint64_t kmer_count = metadata.size();
  cout << "kmer count: " << kmer_count << endl;

  //note the magical '1000' (bad style, I know) make sure you
  //put the "if ((1+j) % 1000 == 0)" block at the very end of the
  //for (uint64_t j=0; ... loop
  KmerNode w;
  uint64_t test, sanity = ~0;
  for (uint64_t j=0; j<kmer_count; j++) {
      w.read(fp);
      //w->print(cout, false); 
//      w.print(cout, true); 
      // at this point you can do whatever you want with the KmerNode,
      //e.g, gather statistics, etc.  Here, we're merely writing
      //the kmers to cout in ascii format.
 //     delete w;
      if ((1+j) % 1000 == 0) {
        assert(fread(&test, sizeof(uint64_t), 1, fp) == 1);
        assert(test == sanity);
      }
  }

  fclose(fp);
}
