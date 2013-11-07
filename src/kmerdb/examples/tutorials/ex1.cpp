#include "../../all_headers.hpp"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;
using namespace metag;

void get_query(int K, string &q) {
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    cerr << "usage: " << argv[0] << " 0|1|2| 0|1 query_count kmer_fn\n";
    cerr << "function: loads a KmerDB and queries using randomly generated kmers;\n";
    cerr << "          kmers may be optionally sorted prior to querying\n";
    cerr << "second argument is: \n";
    cerr << "   0 - KmerDB_exp\n";
    cerr << "   1 - KmerDB_map\n";
    cerr << "   2 - KmerDB_gnu_hash\n";
    cerr << "third argument is: \n";
    cerr << "   0 - do not sort the queries\n";
    cerr << "   1 - sort the queries\n";
    exit(-1);
  }
    int which = atoi(argv[1]);
    int sortme = atoi(argv[2]);
    int q_ct = atoi(argv[3]);
    string fn = argv[4];

    int K;
    StopWatch clock;
    clock.start();
    KmerDB *h_map;
    switch (which) {
      case 0 : h_map = new KmerDB_exp;
               h_map->insert(fn.c_str());
               break;
      case 1 : h_map = new KmerDB_map;
               h_map->insert(fn.c_str());
               break;
               
      case 2 : h_map = new KmerDB_gnu_hash;
               h_map->insert(fn.c_str());
               break;
      default : cerr << "ilegal invocation\n";
                exit(-1);
    }
    K = h_map->getKmerLength();

    cout << "total time to load KmerDB: " << clock.stop() << endl;

    vector<uint64_t> queries(q_ct);
    string q(K, 'x');

    size_t mx = (size_t)pow(4.0, (double)K);
    for (int j=0; j<q_ct; j++) {
      queries[j] = random() % mx;
    }

    if (sortme) {
      sort(queries.begin(), queries.end());
    }

    clock.reset();
    clock.start();
    int found = 0;
    for (int j=0; j<q_ct; j++) {
      if (h_map->find(queries[j]) != 0) {
        ++found;
      }
    }

    double pct = 100.0*(double)found/(double)q_ct;
    cout << "query time: " << clock.stop() << endl;
    cout << "found " << pct << " of queries\n";
}
