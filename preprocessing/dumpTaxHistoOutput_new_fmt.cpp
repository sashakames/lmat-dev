#include "all_headers.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <map>

using namespace std;
using namespace metag;

void usage() {
  cout << "usage:\n"
       "  -i <string>    - tax_histo filename\n"
       "  -f 16|32       - number of bits in tax IDs\n"
       "  -q <int>       - quit after printing this many entries [optional]\n"
       "  -m <string>    - mapping file for 32 - 16 bit tax IDs; optionally given if\n"
       "                   using 16 bit IDs; if given, tax IDs will be converted back\n"
       "                   to their true identities\n";

}

template<typename tid_T>
void doit(string &t_file, string &map_fn, int quit, tid_T junk); 


int main(int argc, char *argv[]) {
  string t_file, map_fn;
  bool prn_help = false;
  int count = 0;
  int which = 0;
  int quit = 0;
  char c;
  while ((c = getopt(argc, argv, "h i:f:q:m:")) != -1) {
    switch (c) {
    case 'm':
      map_fn = optarg;
      break;
    case 'h':
      prn_help = true;
      break;
    case 'q':
      quit = atoi(optarg);
      break;
    case 'i':
      ++count;
      t_file = optarg;
      break;
    case 'f':
      ++count;
      which = atoi(optarg);
      break;
    default:
      cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      prn_help = true;
      break;
    }
  }


  if (prn_help || count != 2) {
    usage();
    exit(-1);
  }

  if (which == 16) {
    uint16_t w = 0;
    doit(t_file, map_fn, quit, w);
  } else {
    uint32_t w = 0;
    doit(t_file, map_fn, quit, w);
  }
}

template<typename tid_T>
void doit(string &t_file, string &map_fn, int quit, tid_T junk) {
  FILE *fp = fopen(t_file.c_str(), "r");
  assert(fp);
  KmerFileMetaData metadata;
  metadata.read(fp);
  assert(metadata.version() == TAX_HISTO_VERSION);

  map<uint16_t, uint32_t> mp;
  if (map_fn.size()) {
    ifstream in(map_fn.c_str());
    assert(in);
    uint16_t a; 
    uint32_t b;
    while (in >> b >> a) {
      mp[a] = b;
    }
  }

  uint64_t test, sanity = ~0;
  kmer_t kmer;
  tid_T tid;
  uint16_t tid_count;
  assert(sizeof(kmer_t) == 8);

  for (size_t j=0; j<metadata.size(); j++) {
    if (quit && (int)j == quit) {
      break;
    }
    assert(fread(&kmer, sizeof(kmer_t), 1, fp) == 1);
    assert(fread(&tid_count, 2, 1, fp) == 1);
    cout << "j: " << j << " kmer: " << kmer << " tid count: " << tid_count << " :: ";

    for (uint16_t h=0; h<tid_count; h++) {
      assert(fread(&tid, sizeof(tid_T), 1, fp) == 1);
      if (mp.size()) {
        if (mp.find(tid) == mp.end()) {
          cout << "failed to find mapping for 16 bit tax ID: " << tid << " which is " << h << " of " << tid_count <<  endl;
          exit(-1);
        }
        cout << mp[tid] << " ";
      } else {
        cout << tid << " ";
      }  
    }
    cout << endl;
    if ((j+1) % TAX_HISTO_SANITY_COUNT == 0) {
      assert(fread(&test, 8, 1, fp) == 1);
      assert(test == sanity);
    }
  }

}
