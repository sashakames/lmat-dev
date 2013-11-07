#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " gid2tid_fn gid_fn \n";
    exit(-1);
  }

  GenomeIdToTaxId genome_to_taxid(argv[1]);
  ifstream in(argv[2]);
  assert(in);

  uint32_t gid;
  int found = 0;
  int not_found = 0;
  while (!in.eof()) {
    in >> gid;
    if (genome_to_taxid.find(gid) != genome_to_taxid.end()) {
      ++found;
    } else {
      ++not_found;
    }
  }  
cout << "found: " << found << endl;
cout << "not found: " << not_found << endl;
}
