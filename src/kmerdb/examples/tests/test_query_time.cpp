#include <iostream>
#include <fstream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

int K = 16;

void getReadKmers(set<uint64_t> &kmers, string &read) {
  cout << "------------------------------------- getReadKmers; read.size: " << read.size() << "\n";
  cout << read << endl;
  kmers.clear();
  char work[1024];
  uint64_t kmer, rc;
  for (size_t j=0; j<read.size()-K; j++) {
    sprintf(work, "%s", read.substr(j, K).c_str());
    kmer = parse_dna::mer_string_to_binary(work, K);
    rc = parse_dna::reverse_complement(kmer, 16);
    //if (kmer < rc) {
      kmers.insert(kmer);
    //} else {
      kmers.insert(rc);
    //}
  }
  cout << "\nkmer count: " << kmers.size() << endl;
}

template <class DB>
void get_percent(kmerdb<DB>& db, set<uint64_t> &kmers) {
  int ct = 0;
  cout << "starting get_percent\n";
  for (set<uint64_t>::const_iterator t = kmers.begin(); t != kmers.end(); t++) {
    cout << "next kmer: " << *t << endl;
    if (db.lookup(*t) != 0) {
      ++ct;
    }
  }
  cout << "percent: " << 100.0 * (double)ct/(double)kmers.size() << endl;
}

template <class DB>
void testme(kmerdb<DB>& db, string &line, int read_len) {
  string read;
  if (line.size() < read_len) return;
  set<uint64_t> kmers;
  cout << "line size: " << line.size() << endl;
  for (size_t j=0; j<line.size()-read_len-1; j++) {
    cout << "testme; j= " << j << endl;
    read = line.substr(j, read_len);
    cout << ">>>>>>>>> " << read.size() << " read len: " << read_len << endl;
    if (read.size() == read_len) {
      getReadKmers(kmers, read);
      cout << "got kmers!\n";
      get_percent(db, kmers);
    }
  }
}

template <class DB>
void query_test(kmerdb<DB>& db, char *argv[]) {
    db.dump();
    exit(-1);

    int len = atoi(argv[4]);
    int num_lines = atoi(argv[5]);
    ifstream in(argv[3]);
    assert(in.is_open());
    string line;
    while(true) {
      line = "";
      getline(in, line);  //header
      if (line.size() == 0) {
        break;
      }
      getline(in, line);  //header
      testme(db, line, len);
    }
    in.close();
}

int main(int argc, char *argv[]) {
  if (argc != 6) {
    cerr << "usage: " << argv[0] << " 1|2|3  kmer_fn query_fn read_len num_lines_to_process\n";
    cerr << "function: tests time for loading a kmerDB\n";
    cerr << "second argument is: \n";
    cerr << "   1 - KmerDB_exp\n";
    cerr << "   2 - KmerDB_map\n";
    cerr << "   3 - KmerDB_gnu_hash\n";
  }
    int which = atoi(argv[1]);
    string fn = argv[2];

    StopWatch clock;
    clock.start();
    switch (which) {
      case 1 : {kmerdb<KmerDB_exp> kdb;
               kdb.ingest(fn.c_str());
               query_test(kdb, argv);
               break;
               }
      case 2 : {kmerdb<KmerDB_map> kdb;
               kdb.ingest(fn.c_str());
               query_test(kdb, argv);
               break;
               }
      case 3 : {kmerdb<KmerDB_gnu_hash> kdb;
               kdb.ingest(fn.c_str());
               query_test(kdb, argv);
               break;
               }
      default : cerr << "illegal invocation\n";
                exit(-1);
    }

    //cout << "time to load KmerDB: " << clock.stop() << endl;
}
