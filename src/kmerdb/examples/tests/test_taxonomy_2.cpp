#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << endl;
    cerr << "function: tests that paths from every TaxNode to root\n"
         << "terminate at the actual root node.\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
    }
  }

  cerr << "testing that paths from every TaxNode to root\n"
       << "terminate at the actual root node.\n";

  char * gid_to_tid_fn = "../../../../../taxonomy/kpath/gid_to_kpath_tax_id.dat";

  TaxTree tax_tree("../../../../../taxonomy/kpath//kpath_taxonomy.dat", gid_to_tid_fn);


  cout << "tax tree size: " << tax_tree.size() << endl;

  TaxNodeHash::const_iterator t2 = tax_tree.end(), tt = tax_tree.end();

  for (TaxNodeHash::const_iterator t = tax_tree.begin(); t != tax_tree.end(); t++) {
    t2 = t;
    while (true) {
      if (t2 == tax_tree.end()) {
        break;
      }

      if (t2->second->id() == t2->second->parent()) {
        break;
      }

      tt = tax_tree.find(t2->second->parent()); 
      if (tt == tax_tree.end()) {
        //ERROR: failed to find parent
        assert(false);
      }
      t2 = tt;
    }

      
    //verify that we're at the true root!
    if (t2->second->name() != "root") {
        cout << t2->second->id() << " " << t2->second->name() << endl;
        const set<uint32_t> &s = t2->second->getChildren();
        cout << "not root! child count: " << s.size() << endl;
        assert(false);
    }
  }

  return 0;
}
