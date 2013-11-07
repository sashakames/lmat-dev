#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
//currently KmerDB is either typedefed or #defined as
//one of the three alternate KmerDB_XX classes, so I'm
//changing this test to reflect that.  

  if (argc != 3) {
    cerr << "usage: " << argv[0] << " kmer_fn output_fn\n";
    cerr << "function: loads a kmerDB, then dumps kmers to file in ascii format\n";
    exit(-1);
  }

    KmerDB kdb;
    kdb.ingest(argv[1]);
    kdb.dumpKmers(argv[2]);

#if 0
  if (argc != 4) {
    cerr << "usage: " << argv[0] << " 1|2|3   0|1 kmer_fn\n";
    cerr << "function: tests time for loading a kmerDB\n";
    cerr << "second argument is: \n";
    cerr << "   1 - KmerDB_exp\n";
    cerr << "   2 - KmerDB_map\n";
    cerr << "   3 - KmerDB_gnu_hash\n";
    cerr << "third argument is: \n";
    cerr << "   0 - load binary file to memory then read kmers\n";
    cerr << "   1 - read kmers directly from file\n";
    exit(-1);
  }
    int which = atoi(argv[1]);
    string fn = argv[3];

    bool read_from_file = false;
    if (strcmp(argv[2], "1")) {
      read_from_file = true;
    }

    KmerFileMetaData metadata;
    metadata.read(fn.c_str());
    //metadata.write();
    //cout << endl;

    StopWatch clock;
    clock.start();
    switch (which) {
      case 1 : {kmerdb<KmerDB_exp> kdb;
               kdb.ingest(fn.c_str());
               kdb.dump("data/test.fa.dump");
               cout << "loaded size: " << kdb.size() << endl;
               break;
               }
      case 2 : {kmerdb<KmerDB_map> kdb;
               kdb.ingest(fn.c_str());
               cout << "loaded size: " << kdb.size() << endl;
               kdb.dump("data/test.fa.dump");
               break;
               }
      case 3 : {kmerdb<KmerDB_gnu_hash> kdb;
               kdb.ingest(fn.c_str(), read_from_file);
               cout << "loaded size: " << kdb.size() << endl;
               kdb.dump("data/test.fa.dump");
               break;
               }
      default : cerr << "illegal invocation\n";
                exit(-1);
    }
    #endif
}
