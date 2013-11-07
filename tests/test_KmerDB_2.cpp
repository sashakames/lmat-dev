#include "../include/all_headers.hpp"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace metag;
using namespace __gnu_cxx;


int main(int argc, char *argv[]) {
  //read binary encoded file into a vector of KmerNode objects
  string fn = "one.fa.bin.loc";
  GidToTid junk;
  vector<KmerNode*> kmers;
  KmerNode::binary_reader(fn.c_str(), kmers, junk);

  //construct and populate a hash table
  KmerDB_gnu_hash_map table(kmers);

  int j = -1;
  for (__gnu_cxx::hash_map<uint64_t, const KmerNode*>::const_iterator t = table.begin(); t != table.end(); t++) {
    ++j;
    if (j == 20) {
      break;
    }
    cout << t->second->getStringKmer() << endl;
  }

}
