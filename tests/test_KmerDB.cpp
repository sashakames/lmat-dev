#include "../include/all_headers.hpp"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  //read binary encoded file into a vector of KmerNode objects
  string fn = "one.fa.bin.loc";

  //construct and populate a hash table
  KmerDB table(fn.c_str(), 16);

  //perform some queries
  KmerNode * kmer_node = table.find("AAAAAAAAAAAAAAAA", 16);
  cout << "find AAAAAAAAAAAAAAAA returned: ";
  if (kmer_node == 0)  cout << " not found\n";
  else  cout << " FOUND\n";

  kmer_node = table.find("TGTAAAGCCCCCAAAA", 16);
  cout << "find TGTAAAGCCCCCAAAA returned: ";
  if (kmer_node == 0)  cout << " not found\n";
  else  cout << " FOUND\n";

  //this should not be in the table
  kmer_node = table.find("TGTAAATTTCCCAAAA", 16);
  cout << "find TGTAAATTTCCCAAAA returned: ";
  if (kmer_node == 0)  cout << " not found\n";
  else  cout << " FOUND\n";

  cout << endl << "constructing a new KmerNode with kmer TGTAAATTTCCCAAAA, and inserting it in the table\n";
  kmer_node = new KmerNode;
  kmer_node->setSize(16);
  uint64_t x = parse_dna::mer_string_to_binary("TGTAAATTTCCCAAAA", 16);
  kmer_node->setKmer(x);
  table.insert(kmer_node);

  //this should be in the table
  kmer_node = table.find("TGTAAATTTCCCAAAA", 16);
  cout << "find TGTAAATTTCCCAAAA returned: ";
  if (kmer_node == 0)  cout << " not found\n";
  else  cout << " FOUND\n";

  kmer_node = table.find("TGTAAAGCCCCCAAAA", 16);
  cout << "\ncontents of the kmer node for TGTAAAGCCCCCAAAA:\n";
  kmer_node->write();
}
