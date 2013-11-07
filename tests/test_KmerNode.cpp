#include "../include/all_headers.hpp"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  string fn = "one.fa.bin.loc";
  GidToTid junk;
  vector<KmerNode*> kmers;
  KmerNode::binary_reader(fn.c_str(), kmers, junk);

  for (size_t h=0; h<kmers.size(); h++) {
    kmers[h]->print(cout);
  }
}
