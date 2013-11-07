#include "../all_headers.hpp"
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

  //get the metata
  FILE *fp = Utils::openReadFile(argv[1]);
  KmerFileMetaData metadata;
  metadata.read(fp);
  uint64_t sz = metadata.size();
  cout << "current file offset: " <<ftell(fp) << endl;
  cout << "\nmetadata:\n";
  metadata.write();

  uint64_t sanity = ~0;
  uint64_t mark;
  KmerNode w;
  cerr << "starting to read kmers at " << ftell(fp) << endl;
  uint64_t k_ct = 0;
  uint64_t m_ct = 0;
  for (uint64_t j=0; j<sz; ++j) {
     w.read(fp);
     if ((j+1) % 1000 == 0) {
       ++m_ct;
       assert(fread(&mark, 8, 1, fp) == 1);
       if (mark != sanity) {
         cerr << "error check failed at j= " << j << " of " << sz << endl;
         cerr << "read " << mark << " should have been: " << sanity << endl;
         exit(-1);
       }
     }
     ++k_ct;
  }

  cout << "verified " << k_ct << " kmers and " << m_ct << " sanity checks\n";
}
