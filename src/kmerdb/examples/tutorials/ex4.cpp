#include "../../all_headers.hpp"
#include <iostream>
#include <string>

using namespace std;
using namespace metag;
using namespace __gnu_cxx;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << " kmerdb_data_file\n";
    cerr << "function: loads a KmerDB; loads a taxonomy; prints out\n";
    cerr << "          TaxNodes associated with the genomes in the KmerNodes\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
  }

  //load KmerDB
  KmerDB db;
  KmerFileMetaData metadata;
  metadata.read(argv[1]);
  uint64_t kmer_count = metadata.size();
  db.resize(kmer_count);

  GenomeIdToTaxId taxid_lookup("../../../../../taxonomy/kpath//gid_to_kpath_tax_id.dat");
  TaxNodeIDArray tax_id_data;
  db.ingest(argv[1], &taxid_lookup, &tax_id_data);
  cout << db.size() << endl << endl;

  set<uint32_t> tid;
  for (KmerDB::iterator t = db.begin(); t != db.end(); t++) {
    t->second->getTaxIDs(tid, &tax_id_data);
    if (tid.find(9404) != tid.end()) {
      cout << t->second->getKmer() << " has tid 9404\n";
    }
  }

exit(0);

  TaxTree tt("../../../../../taxonomy/kpath/kpath_taxonomy.dat");


      TaxNode *tn = tt.find(9404)->second;
      tn->write();

  //exit(0);

  //count the number of kmers associated with every TaxNode
  set<uint32_t> taxid;
  for (KmerDB::iterator t = db.begin(); t != db.end(); t++) {
    KmerNode *n = t->second;
    n->getTaxIDs(taxid, &tax_id_data);
    for (set<uint32_t>::iterator t2 = taxid.begin(); t2 != taxid.end(); t2++) {
      TaxNode *tn = tt.find(*t2)->second;
      tn->m_count += 1;
    }
  }

  for (TaxTree::iterator t = tt.begin(); t != tt.end(); t++) {
    TaxNode *tn = t->second;
    if (tn->m_count || !tn->isLeaf()) {
      cout << "found tax node " << tn->id() << " that has genomes, but is not a leaf\n";
    }
  }

  exit(0);

  cout << "id: 8912\n";
  if (tt.find(8912) == tt.end()) cout << "8912 not found\n";
  else tt.find(8912)->second->write();
  cout << endl;

  cout << "id: 366079\n";
  if (tt.find(366079) == tt.end()) cout << "366079 not found\n";
  else tt.find(366079)->second->write();
  cout << endl;

  cout << "id: 8897\n";
  if (tt.find(8897) == tt.end()) cout << "8897 not found\n";
  else tt.find(8897)->second->write();
  cout << endl;

  cout << "id: 474521\n";
  if (tt.find(474521) == tt.end()) cout << "474521 not found\n";
  else tt.find(474521)->second->write();
  cout << endl;

  cout << "id: 270225\n";
  if (tt.find(270225) == tt.end()) cout << "270225 not found\n";
  else tt.find(270225)->second->write();
  cout << endl;

  cout << "id: 290612\n";
  if (tt.find(290612) == tt.end()) cout << "290612 not found\n";
  else tt.find(290612)->second->write();
  cout << endl;

  cout << "8912\n";
  const std::vector<uint32_t> * p = tt.find(8912)->second->getPathToRoot();
  for (vector<uint32_t>::const_iterator t = (*p).begin(); t != (*p).end(); t++) {
    cout << *t << " ";
  }
  cout << endl << endl;

  cout << "8897\n";
  p = tt.find(366079)->second->getPathToRoot();
  for (vector<uint32_t>::const_iterator t = (*p).begin(); t != (*p).end(); t++) {
    cout << *t << " ";
  }
  cout << endl << endl;

  cout << "8897\n";
  p = tt.find(8897)->second->getPathToRoot();
  for (vector<uint32_t>::const_iterator t = (*p).begin(); t != (*p).end(); t++) {
    cout << *t << " ";
  }
  cout << endl << endl;

  cout << "474521\n";
  p = tt.find(474521)->second->getPathToRoot();
  for (vector<uint32_t>::const_iterator t = (*p).begin(); t != (*p).end(); t++) {
    cout << *t << " ";
  }

  cout << endl << endl;
  cout << "270225\n";
  p = tt.find(270225)->second->getPathToRoot();
  for (vector<uint32_t>::const_iterator t = (*p).begin(); t != (*p).end(); t++) {
    cout << *t << " ";
  }
  cout << endl << endl;

  cout << "290612\n";
  p = tt.find(290612)->second->getPathToRoot();
  for (vector<uint32_t>::const_iterator t = (*p).begin(); t != (*p).end(); t++) {
    cout << *t << " ";
  }

  exit(0);


/*
  TaxTree tt("/usr/mic/post1/metagenomics/taxonomy/kpath_taxonomy.dat");

  for (_kmer_hash::const_iterator t = db.begin(); t != db.end(); t++) {
    t->second->write(cout, &tt, &tax_id_data);

  }
  */

  exit(0);


  //iterate over the KmerNodes: print the taxid for each genome
  /*
  GidHash::const_iterator iter;
  for (hash_map<uint64_t, KmerNode*>::const_iterator t = table.begin(); t != table.end(); t++) {
    t->second->print();
    _gen_set gids = t->second->getGenomes();
    for (_gen_set::const_iterator t2 = gids.begin(); t2 != gids.end(); t2++) {
      iter = taxid_lookup.find(*t2);
      if (iter == taxid_lookup.end()) {
        cout << "tax node ID wasn't found for gid " << *t2 << endl;
      } else {
        cout << "    gid: " << iter->first << " tax node id: " << iter->second << endl;
      }
    }
  }
  */

  exit(0);
  /*
  for (TaxNodeHash::const_iterator t = tt.begin(); t != tt.end(); t++) {
    t->second->write();
  }
  */
}
