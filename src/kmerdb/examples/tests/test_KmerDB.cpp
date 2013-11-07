#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " 1|2|3  kmer_fn\n";
    cerr << "function: tests time for loading a kmerDB\n";
    cerr << "second argument is: \n";
    cerr << "   1 - KmerDB_exp\n";
    cerr << "   2 - KmerDB_map\n";
    cerr << "   3 - KmerDB_gnu_hash\n";
    exit(-1);
  }

    assert(false);
#if 0
    int which = atoi(argv[1]);
    string fn = argv[2];


    StopWatch clock;
    clock.start();
    switch (which) {
      case 1 : {kmerdb<KmerDB_exp> kdb;
               kdb.ingest(fn.c_str());
               break;
               }
      case 2 : {kmerdb<KmerDB_map> kdb;
               kdb.ingest(fn.c_str());
               break;
               }
      case 3 : {kmerdb<KmerDB_gnu_hash> kdb;
               kdb.ingest(fn.c_str());
               break;
               }
      default : cerr << "illegal invocation\n";
                exit(-1);
    }

    cout << "total time: " << clock.stop() << endl;
#endif
}
