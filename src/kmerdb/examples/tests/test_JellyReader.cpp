#include <iostream>
#include <string>
#include "../../JellyReader.hpp"
#include "../../all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << "usage: " << argv[0] << " K fasta_fn hash_size\n";
    cerr << "function: outputs KmerNodes directly from jellylist code (class JellyReader)\n";
    exit(-1);
  }

  int K = atoi(argv[1]);
  uint64_t hash_size = atol(argv[3]);

  KmerDB *db = new KmerDB;

  StopWatch clock;
  clock.start();
  JellyReader jr(db, K, argv[2], hash_size, 32);
  cout << "db size: " << db->size() << endl;
  cout << "time to load: " << clock.stop() << endl;

//  db->dump("test.fa.JellyReader.dump");
}
