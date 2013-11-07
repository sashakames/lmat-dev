#include <iostream>
#include <string>
#include "all_headers.hpp"
#include "test_common.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << " num_to_find -v\n";
    cerr << "function: finds LCA for num_to_find pairs of TaxTree leafs\n";
    exit(-1);
}

string base("../../../../../taxonomy/kpath/");

int main(int argc, char *argv[]) {
  if (argc < 2) {
      usage(argv);
      exit(-1);
  }
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
      exit(-1);
    }
  }

  int num_to_find = atoi(argv[1]);
  bool verbose = argc == 3 ? true : false;

  TaxTree tax_tree(tax_tree_fn);
  TaxNodeIDArray tax_id_data;
  cout << "tax tree size: " << tax_tree.size() << endl;

  GenomeIdToTaxId taxid_lookup(tax_lookup_fn);


/*
  KmerDB db;
  KmerFileMetaData metadata;
  metadata.read(argv[1]);
  uint64_t kmer_count = metadata.size();
  db.resize(kmer_count);
  db.ingest(argv[1], &taxid_lookup, &tax_id_data);
  */

  TaxNode *tn_1, *tn_2;
  const TaxNode *lca;
  set<uint32_t> s;
  vector<uint32_t> p_1, p_2;
  int ct = 0;
  StopWatch clock;
  clock.start();


  int failed = 0;
  bool stopme = false;
  for (TaxTree::const_iterator t = tax_tree.begin(); t != tax_tree.end(); t++) {
    if (stopme) break;
    tn_1 = t->second;
    tn_1 = t->second;
    if (tn_1->id() == 1) continue;
    if (! tn_1->isLeaf() == 1) continue;
    s.clear();
    s.insert(tn_1->id());
    tax_tree.getPathToRoot(tn_1->id(), p_1);

    for (TaxTree::const_iterator t2 = tax_tree.begin(); t2 != tax_tree.end(); t2++) {
      tn_2 = t2->second;
      if (tn_2->id() == 1) continue;
      if (! tn_2->isLeaf() == 1) continue;
      ++ct;
      tax_tree.getPathToRoot(tn_2->id(), p_2);

      //cout << ct << " " << num_to_find << endl;
      if (ct >= num_to_find) {
        stopme = true;
        break;
      }  

      if (tn_1->id() != tn_2->id()) {
        s.insert(tn_2->id());
        lca = tax_tree.lowestCommon(s);
        if (!lca) {
         ++failed;
        }
        if (verbose) {
        if (!lca) {
          cout << "failed to find lca\n";
        } else {
        cout << "found lca!\n";
        const set<uint32_t> &c_1 = lca->getChildren();
        cout << endl;
        cout << "first leaf id: " << tn_1->id() << " :: ";
        for (size_t j=0; j<p_1.size(); j++) {
          cout << p_1[j] << " -> ";
        }
        cout << endl;
        cout << "second leaf id: " << tn_2->id() << " :: ";
        for (size_t j=0; j<p_2.size(); j++) {
          cout << p_2[j] << " -> ";
        }
        cout << endl;

        cout << "lca: " << lca->id() << endl;
        cout << "lca children: ";

        for (set<uint32_t>::const_iterator t3 = c_1.begin(); t3 != c_1.end(); t3++) {
          cout << *t3 << " ";
        }
        cout << endl << endl;
        } //if verbose
        s.erase(tn_2->id());
      }
    }  
  }
  }
  cout << "time to find " << num_to_find << " LCAs: " << clock.stop() << endl;


  return 0;
}
