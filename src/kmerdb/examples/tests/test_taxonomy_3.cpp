#include <iostream>
#include <set>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << endl;
    cerr << "function: tests no interior node has genomes associated with it\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
    }
  }

  cerr << "testing that no interior node has genomes associated with it\n";

  char * gid_to_tid_fn = "../../../../../taxonomy/kpath/gid_to_kpath_tax_id.dat";

  TaxTree tax_tree("../../../../../taxonomy/kpath//kpath_taxonomy.dat", gid_to_tid_fn);
  GenomeIdToTaxId genome_to_taxid(gid_to_tid_fn);


  //build set of tax IDs that are associated with genome IDs
  set<uint32_t> taxid_with_genome_ids;
  for (GenomeIdToTaxId::const_iterator t = genome_to_taxid.begin(); t != genome_to_taxid.end(); t++) {
    uint32_t tid = t->second;

    //note: multiple genomes can be associated with the same tid
    taxid_with_genome_ids.insert(tid);
  }

  cout << "number of genomes: " <<  genome_to_taxid.size() << endl; 
  cout << "number of tax IDs that are associated with genomes: " << taxid_with_genome_ids.size() << endl;

  TaxTree::const_iterator iter;
  for (set<uint32_t>::const_iterator t = taxid_with_genome_ids.begin(); t != taxid_with_genome_ids.end(); t++) {

    iter = tax_tree.find(*t);

    //I think this is actually tested for elsewhere ... but can't 
    //be too safe
    assert(iter != tax_tree.end());
    assert(iter->second->isLeaf());
  }

  return 0;
}
