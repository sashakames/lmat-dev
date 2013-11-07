#include <iostream>
#include <vector>
#include <cassert>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains output from jellylist.bin\n";
    cerr << "function: instantiates kmer from file, and from a\n";
    cerr << "          memory buffer, and tests for equality\n";
    exit(-1);
  }

  //read offsets file into memory
  char buf[1024];
  sprintf(buf, "%s.offsets", argv[1]);
  uint64_t *offsets = (uint64_t*)Utils::readFile(buf);

  //read kmer data file into memory
  char *kmer_data = Utils::readFile(argv[1]);

  //get the kmer length and number of kmers in the file
  KmerFileMetaData metadata;
  metadata.read(argv[1]);
  uint64_t kmer_count = metadata.size();

  FILE *fp = Utils::openReadFile(argv[1]);
  fseek(fp, metadata.dataStart(), SEEK_SET);
  assert((size_t)offsets[0] == (size_t)ftell(fp));

  KmerNode w, *ww;
  KmerNode *nds = new KmerNode[kmer_count];
  uint64_t sanity_mark;
  for (uint64_t j=0; j<kmer_count; j++) {
      w.read(fp);
      ww = &nds[j];
      ww->read(&kmer_data[0]+offsets[j]);

      if (! (w == *ww)) {
        cerr << "cur file offset: " << ftell(fp) << endl;
        cerr << "KmerNodes are not equal; test failed\n";
        cerr << endl;
        cerr << "j: " << j << " offsets[j]: " << offsets[j] << " ftell(fp): " << ftell(fp) << endl;
        exit(-1);
      }
      if ((j+1) % 1000 == 0) {
        assert(fread(&sanity_mark, 8, 1, fp) == 1);
      }
  }

  fclose(fp);
}
