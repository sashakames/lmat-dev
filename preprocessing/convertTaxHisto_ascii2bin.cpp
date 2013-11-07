#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include "all_headers.hpp"

using namespace metag;
using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " input_fn  output_dir\n";
    cerr << "output fn is: output_dir/input_fn.bin\n";
    exit(-1);
  }

  //open input file
  ifstream in(argv[1]);
  assert(in);
  Utils::skipHeaderLines(in);

  //open output file
  string f(argv[1]);
  int j = f.rfind('/');
  if (j == string::npos) {
    cerr << "rfind('/') failed for: " << f << endl;
  }
  assert(j != string::npos);
  string out_fn = f.substr(j+1);
  char buf[1024];
  sprintf(buf, "%s/%s.bin", argv[2], out_fn.c_str());
  FILE *out = fopen(buf, "wb");
  assert(out);
  
  string line;
  uint64_t kmer;
  uint16_t genome_count, tid_count, present, tuple_count;
  uint32_t tid;

  uint64_t total = 0;
  uint64_t test = ~0;
  while (true) {
    //if (in.eof() || !in.good()) break;

    //need special handling since there's probably a space after the
    //last entry on each line
    kmer = ~0;
    in >> kmer;
    if (kmer == test) {
      break;
    }

    ++total;
    in >> genome_count;
    in >> tid_count;
    in >> tuple_count;
    assert(fwrite(&kmer,8, 1, out) == 1);
    assert(fwrite(&genome_count,2, 1, out) == 1);
    assert(fwrite(&tid_count,2, 1, out) == 1);
    assert(fwrite(&tuple_count,2, 1, out) == 1);
    for (uint16_t k=0; k<tuple_count; k++) {
      in >> tid;
      in >> present;
      assert(fwrite(&tid,4, 1, out) == 1);
      assert(fwrite(&present,2, 1, out) == 1);
    }
  }
  cout << "num kmers read: " << total << endl;
  in.close();
}
