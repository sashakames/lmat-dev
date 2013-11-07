#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << endl;
    cerr << "function: tests that all mfe IDs have a corresponding node in a TaxTree";
    exit(-1);
}

int main(int argc, char *argv[]) {
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
    }
  }

  cout << "testing that all mfe IDs in microbe2.kpath_id.seqid \n"
       << "have a corresponding TaxTree node\n";

  char * gid_to_tid_fn = "../../../../../taxonomy/kpath/gid_to_kpath_tax_id.dat";

  TaxTree tax_tree("../../../../../taxonomy/kpath//kpath_taxonomy.dat", gid_to_tid_fn);
  GenomeIdToTaxId genome_to_taxid(gid_to_tid_fn);

  ifstream in("../../../../../taxonomy/kpath//microbe2.kpath_id.seqid");
  uint32_t kpath_id;
  uint32_t mfe_id;
  assert(in.is_open());
  uint32_t ct = 0;
  while (!in.eof()) {
    in >> kpath_id;
    in >> mfe_id;
    ++ct;
  }
  assert(tax_tree.find(mfe_id) != tax_tree.end());
  cout << "num mfe IDs read from file: " << ct << endl;
  cout << "tax_tree.size(): " << tax_tree.size() << endl;
  return 0;
}
