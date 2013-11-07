#include <iostream>
#include <fstream>
#include <set>
#include <cassert>
#include <string>
#include <string.h>
#include <stdint.h>
#include "all_headers.hpp"
#include <ext/hash_set>

#define TID_T uint32_t

using namespace std;
using namespace metag;

void usage(char *argv[]) {
  cout << "Usage:\n"
       << " -i <string>  - sequence\n"
       << " -k <int>     - kmer length\n"
       << " -p <string>  - kmer base name prefix\n"
       << " -c <int>     - tax_histo file count\n"
       << " -t           - taxid\n"
       << " -h           - print help and exit\n\n"
       << "function: will test for every kmer (and it's rc)\n"
       << "          in the sequence; if found, prints the entry\n"
       << "sample invocation: verifyKmerEntries -i GGGCTGATATTCTTAAAACCAAATATTTAAAAAATGTAAATATGTTAATA -k 20 -p /p/lscratchd/hysom/m5_hera/m5. -c 256\n";

}

int main(int argc, char *argv[]) {
  const string opt_string="i:k:p:c:h";
  string seq, prefix;
  size_t kmer_len = 0; 
  int t_count = 0;
  bool prn_help = false;
  int count = 0;
  char c;
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
      case 'i':
        seq = optarg;
        ++count;
        break;
      case 'k':
        kmer_len = atoi(optarg);
        ++count;
        break;
      case 'p':
        prefix = optarg;
        ++count;
        break;
      case 'c':
        t_count = atoi(optarg);
        ++count;
        break;
      case 'h':
        prn_help = false;
        break;
      default:
       cerr << "Unrecognized option: "<<c<<", ignore."<<endl;
       prn_help = true;
       break;
    }
  }

  if (prn_help || count != 4) {
    usage(argv);
    exit(-1);
  }

  //build a vector s.t. starts[j] contains the first kmer in the j-th file;
  //since kmers are sorted, this will enable us to map a kmer to the tax_histo
  //file in which it exists (if it exists)
  string line;
  char buf[1024];
  vector<uint64_t> starts;
  uint64_t kmer;
  KmerFileMetaData metadata;
  for (int j=0; j<t_count; ++j) {
    sprintf(buf, "%s%d", prefix.c_str(), j);
    FILE *fp = fopen(buf, "r");
    if (!fp) {
      cerr << "failed to open " << buf <<" for reading\n";
      exit(9);
    }  
    metadata.read(fp);
    assert(fread(&kmer, 8,1,fp) == 1);
    starts.push_back(kmer);
  }

  string rev;
  string s(seq);
  uint64_t test, sanity = ~0;
  for (size_t j=0; j<(int)s.size()-(int)kmer_len+1; j++) {
    string km = s.substr(j,kmer_len);
    cout << "\nlooking for kmer: " << km << endl;
    uint64_t fwd = Encoder::encode(km);
    uint64_t rc = Encoder::rc(fwd, kmer_len);
    Encoder::decode(rc, kmer_len, rev);
    size_t idx;
    if (rc < fwd) {
      fwd = rc;
    }

    //find the kmer file that may contain the kmer
    int use_me = -1;
    for (idx = 0; idx < starts.size()-1; idx++) {
      if (fwd >= starts[idx] && fwd < starts[idx+1]) {
        use_me = idx;
        break;
      }
    }
    if (use_me == -1) idx = starts.size()-1;
    sprintf(buf, "%s%d", prefix.c_str(), use_me);
    FILE *in = fopen(buf, "r");
    assert(in);
    metadata.read(in);
    KmerNode<uint32_t>  nd;
    cout << "looking in file: " << buf << endl;
    bool found = false;
    for (size_t j=0; j<metadata.size(); j++) {
      nd.read(in);
      if (nd.getKmer() == fwd) {
        found = true;
        const set<TID_T> &tids = nd.getTaxIDs();
        cout << "kmer: " << nd.getKmer() << " tax IDs: ";
        for (set<TID_T>::const_iterator t = tids.begin(); t != tids.end(); t++) {
          cout << *t << " ";
        }
        cout << endl;
        break;
      }

      if ((j+1) % KMER_SANITY_COUNT == 0) {
        assert(fread(&test, 8, 1, in) == 1);
        assert(test == sanity);
      }
    }
    if (!found) cout << "  NOT FOUND\n";

  }
}


