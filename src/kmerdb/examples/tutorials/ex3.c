#include "../../all_headers.hpp"
#include <iostream>
#include <string>

using namespace std;
using namespace metag;
using namespace __gnu_cxx;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << " output_fn fasta_1.fa [fasta_2.fa ...]\n";
    cerr << "function: loads a KmerDB with one or more files,\n";
    cerr << "          then dumps kmers to output.fn is ascii format\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
  }
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
    }
  }

  string output_fn = argv[1];

  KmerDB_gnu_hash h_map;
  for (int j=2; j<argc; j++) {
    StopWatch cl;
    cl.start();
    KmerFileMetaData metadata;
    metadata.read(argv[j]);
    metadata.write();
    h_map.KmerDB::insert(argv[j], false);
    cout << "\nload time: " << cl.stop() << endl << endl;
  }
}
