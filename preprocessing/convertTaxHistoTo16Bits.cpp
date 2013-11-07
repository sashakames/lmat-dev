#include "all_headers.hpp"
#include <iostream>
#include <cstdlib>
#include <stdint.h>
#include <fstream>
#include <map>
#include <cassert>

using namespace std;
using namespace metag;

void usage() {
  cout << "usage (all are required):\n"
       " -d <string>  - input directory\n"
       " -b <string>  - basename\n"
       " -s <string>  - suffix\n"
       " -f <int>     - first count\n"
       " -l <int>     - last count\n"
       " -o <string>  - output directory\n"
       " -m <string>  - 32 to 16 bit mapping file\n";
}

int main(int argc, char *argv[]) {
  bool prn_help = false;
  int count = 0;
  string input_dir, output_dir, base, suffix, map_fn;
  const string opt_string="d:b:s:f:l:o:m:h ";
  char c;
  int first = 0, last = 0;
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
    case 'h':
      prn_help = true;
      break;
    case 'b':
      ++count;
      base = optarg;
      break;
    case 's':
      ++count;
      suffix = optarg;
      break;
    case 'f':
      ++count;
      first = atoi(optarg);
      break;
    case 'l':
      ++count;
      last = atoi(optarg);
      break;
    case 'o':
      ++count;
      output_dir = optarg;
      break;
    case 'm':
      ++count;
      map_fn = optarg;
      break;
    case 'd':
      ++count;
      input_dir = optarg;
      break;
    default:
      cerr << "Unrecognized option: "<<c<<", ignore."<<endl;
      prn_help = true;
      break;
    }
  }
  cout << count << endl;

  if (count != 7 || prn_help) {
    usage();
    exit(1);
  }

  //read in tid mapping
  map<uint32_t, uint16_t> mp;
  uint32_t old;
  uint16_t new_tid;
  ifstream in(map_fn.c_str());
  assert(in);
  while (in >> old >> new_tid) {
    mp[old] = new_tid;
  }
  in.close();

  uint64_t kmer;
  uint16_t tid_count, tid;
  uint32_t t;
  uint64_t test, sanity = ~0;
  for (int j=first; j<last; j++) {
    char buf[1024];
    sprintf(buf, "%s/%s%d%s", input_dir.c_str(), base.c_str(), j, suffix.c_str());
    cout << "processing: " << buf << endl;
    FILE *in = fopen(buf, "r");
    assert(in);
    KmerFileMetaData metadata;
    metadata.read(in);
    sprintf(buf, "%s/%s%d%s", output_dir.c_str(), base.c_str(), j, suffix.c_str());
    cout << "opening " << buf << " for writing\n";
    FILE *out = fopen(buf, "w+");
    assert(out);
    metadata.write(out);

    for (int i=0; i<metadata.size(); i++) {
      assert(fread(&kmer, 8, 1, in) == 1);
      assert(fwrite(&kmer, 8, 1, out) == 1);

      assert(fread(&tid_count, 2, 1, in) == 1);
      assert(fwrite(&tid_count, 2, 1, out) == 1);

      for (int k=0; k<tid_count; k++) {
        assert(fread(&t, 4, 1, in) == 1);
        if (mp.find(t) == mp.end()) {
          cout << "failed to find tid mapping for " << t << " for j= " << j << " which is " << k << " of " << tid_count << endl;
          exit(-1);
        }
        tid = mp[t];
        assert(fwrite(&tid, 2, 1, out) == 1);
      }

      if ((i+1) % TAX_HISTO_SANITY_COUNT == 0) {
        assert(fread(&test, 8, 1, in) == 1);
        assert(test == sanity);
        assert(fwrite(&sanity, 8, 1, out) == 1);
      }
    }
    fclose(out);
    fclose(in);
  }
}
