#include "../all_headers.hpp"
#include <iostream>
#include <vector>

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains output from jellylist.bin\n";
    cerr << "output file is: input_fn.sorted\n";
    exit(-1);
  }

  cout << "getting metadata from filename\n";
  KmerFileMetaData metadata;
  metadata.read(argv[1]);
  metadata.write();
  cout << endl;

  cout << "getting metadata from opened file\n";
  FILE *fp = Utils::openReadFile(argv[1]);
  KmerFileMetaData metadata2;
  metadata2.read(fp);
  metadata2.write();

  char buf[1024];
  sprintf(buf, "%s.offsets", argv[1]);
  uint64_t *offsets = (uint64_t*)Utils::readFile(buf);
  cout << "offsets file says 1st kmer starts at offset: " << offsets[0] << endl;
  delete [] offsets;
}
