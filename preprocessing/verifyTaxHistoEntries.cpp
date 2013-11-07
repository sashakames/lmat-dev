#include <iostream>
#include <fstream>
#include <set>
#include <cassert>
#include <string>
#include <string.h>
#include <stdint.h>
#include "all_headers.hpp"
#include <ext/hash_set>

using namespace std;
using namespace metag;

void usage(char *argv[]) {
  cout << "Usage:\n"
       << " -i <string>  - sequence\n"
       << " -k <int>     - kmer length\n"
       << " -p <string>  - tax_histo base name prefix\n"
       << " -s <string>  - tax_histo base name suffix\n"
       << "                if tax histo files are: <path>/m3.int.bin.no_loc.20.partition_##.tax_histo,\n"
       << "                then -p <path>/m3.int.bin.no_loc.20.partition_ and -s .tax_histo\n"
       << " -c <int>     - tax_histo file count\n"
       << " -t           - taxid\n"
       << " -h           - print help and exit\n\n"
       << "function: will tax_histo output for every kmer (and it's rc)\n"
       << "          in the sequence; reports if found, and, if found, if\n"
       << "          the taxid maps to the kmer\n\n"
       << "sample invocation: verifyTaxHistoEntries -i GGGCTGATATTCTTAAAACCAAATATTTAAAAAATGTAAATATGTTAATA -k 20 -p /p/lscratchd/hysom/m3_prefix/m3.fa.int.k=20. -s .tax_histo -c 256 -t 8122\n";

}

int main(int argc, char *argv[]) {
  const string opt_string="i:k:p:s:c:t:h";
  string seq, prefix, suffix;
  size_t kmer_len = 0, taxid, t_count = 0;
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
      case 's':
        suffix = optarg;
        ++count;
        break;
      case 'c':
        t_count = atoi(optarg);
        ++count;
        break;
      case 't':
        taxid = atoi(optarg);
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

  if (prn_help || count != 6) {
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
  for (int j=0; j<t_count; ++j) {
    sprintf(buf, "%s%d%s", prefix.c_str(), j,suffix.c_str());
    ifstream in2(buf);
    //assert(in2);
    if (!in2) {
      cout << "failed to read " << buf << endl;
      continue;
    }
    getline(in2, line);
    getline(in2, line);
    in2 >> kmer;
    starts.push_back(kmer);
    in2.close();
  }
  uint64_t m = ~0;
  starts.push_back(m);


  string rev;
  string s(seq);
  char buf2[1024];
  for (size_t j=0; j<(int)s.size()-(int)kmer_len+1; j++) {
    string km = s.substr(j,kmer_len);
    uint64_t fwd = Encoder::encode(km);
    uint64_t rc = Encoder::rc(fwd, kmer_len);
    Encoder::decode(rc, kmer_len, rev);
    size_t idx;
    bool kontinue = true;

    cout << "------------------------------------------------------\n";
    cout << "           " << fwd<<" "<<rc << endl;
    if (rc < fwd) cout << "rc is canonical\n";
    else cout << "fwd is canonical\n";
    /*
    if (rc < fwd) {
      fwd = rc;
      cout << "using: RC: " << rc << "\n";
    }  else {
      cout << "using: FWD: " << fwd << "\n";
    }
    */

    //find the tax_histo file that may contain the kmer
    int use_me = -1;
    for (idx = 0; idx < starts.size()-1; idx++) {
      if (fwd >= starts[idx] && fwd < starts[idx+1]) {
        use_me = idx;
        break;
      }
    }
    if (use_me == -1) idx = starts.size()-1;

    sprintf(buf, "%s%d%s", prefix.c_str(), (int)idx,suffix.c_str());
    Encoder::rc(km, rev);
    cout << "           kmer: " << km << " rc: " << rev << endl;
    cout << "           looking in: " << buf << endl;
    sprintf(buf2, "egrep %llu %s ", fwd, buf);
    cout << ">>>>>>> " << buf2 << endl;
    system(buf2);

    use_me = -1;
    for (idx = 0; idx < starts.size()-1; idx++) {
      if (!kontinue) break;
      if (rc >= starts[idx] && rc < starts[idx+1]) {
        use_me = idx;
        break;
      }
    }
    if (use_me == -1) idx = starts.size()-1;
    sprintf(buf, "%s%d%s", prefix.c_str(), (int)idx,suffix.c_str());
    cout << "           looking in: " << buf << endl;
    sprintf(buf2, "egrep %llu %s ", fwd, buf);
    //cout << ">>>>>>> " << buf2 << endl;
    system(buf2);
  }
}


