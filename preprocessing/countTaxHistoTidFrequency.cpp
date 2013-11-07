#include "all_headers.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <map>

using namespace std;
using namespace metag;

void usage(const char *nm) {
  cout << "usage:\n"
       "  -b <int>       - first tax_histo file\n"
       "  -l <int>       - last tax_histo file\n"
       "  -p <string>    - prefix\n"
       "  -s <string>    - suffix\n"
       "\n"
       "example:\n"
       << nm << " -b 0 -l 256 -p /p/lscratchd/hysom/m5/m5. -s .tax_histo.bin\n";

}


int main(int argc, char *argv[]) {
  bool prn_help = false;
  int count = 0;
  string prefix, suffix;
  int start=0, end=0;

  char c;
  while ((c = getopt(argc, argv, "b:l:p:s:")) != -1) {
    switch (c) {
    case 'b':
      ++count;
      start = atoi(optarg);
      break;
    case 'l':
      ++count;
      end = atoi(optarg);
      break;
    case 'p':
      ++count;
      prefix = optarg;
      break;
    case 's':
      ++count;
      suffix = optarg;
      break;
    default:
      cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      prn_help = true;
      break;
    }
  }


  if (prn_help || count != 4) {
    usage(argv[0]);
    exit(-1);
  }

  map<int,int> f;
  for (int j=start; j<end; j++) {
    char buf[1024];
    sprintf(buf, "%s%d%s", prefix.c_str(), j, suffix.c_str());
    cerr << "opening for read: " << buf << endl;
    FILE *fp = fopen(buf, "r");
    assert(fp);
    KmerFileMetaData metadata;
    metadata.read(fp);
    assert(metadata.version() == TAX_HISTO_VERSION);

    uint64_t test, sanity = ~0;
    kmer_t kmer;
    uint32_t  tid;
    uint16_t tid_count;
    assert(sizeof(kmer_t) == 8);
    for (size_t h=0; h<metadata.size(); h++) {
      assert(fread(&kmer, sizeof(kmer_t), 1, fp) == 1);
      assert(fread(&tid_count, 2, 1, fp) == 1);
      if (f.find(tid_count) == f.end()) {
        f[tid_count] = 0;
      }
      f[tid_count] += 1;

      for (uint16_t i=0; i<tid_count; i++) {
        assert(fread(&tid, 4, 1, fp) == 1);
      }

      if ((h+1) % TAX_HISTO_SANITY_COUNT == 0) {
        assert(fread(&test, 8, 1, fp) == 1);
        assert(test == sanity);
      }
    }
    fclose(fp);
  }

  uint64_t total = 0;
  for (map<int,int>::iterator t = f.begin(); t != f.end(); t++) {
    total += t->second;
  }
  cout << "\n\ntotal tid count: " << total << endl << endl;

  cout << "tid_count   kmer_count   frequency   running _freq.\n";
  double z = 0;
  for (map<int,int>::iterator t = f.begin(); t != f.end(); t++) {
    z += 100.0*t->second/total;
    printf("%10d %13d   %3.4f%%    %3.4f%%\n", t->first, t->second, 100.0*t->second/total, z);
    //cout << t->first << " " << t->second << " " << 100.0*t->second/total << "%" << endl;
  }

}
