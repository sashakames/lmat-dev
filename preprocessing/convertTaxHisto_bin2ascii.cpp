#include "all_headers.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

using namespace std;
using namespace metag;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " input_fn output_fn\n";
    cerr << "where: input_fn contains tax_histo binary output\n";
    exit(-1);
  }

  cout << "opening for input: " << argv[1] << endl;
  FILE *in = fopen(argv[1], "rb");
  assert(in);

  uint64_t kmer;
  uint16_t genome_count, tid_count;
  map<uint32_t, uint16_t> tids;
  ofstream out(argv[2]);
  assert(out);

  fseek(in, 0, SEEK_END);
  long end = ftell(in);
  fseek(in, 0, SEEK_SET);

  while (true) {
    Utils::readTaxHisto_bin(kmer, genome_count, tid_count, tids, in);
    out << kmer<<" "<< genome_count << " " << tid_count << " " << tids.size();
    for (map<uint32_t, uint16_t>::const_iterator t = tids.begin(); t != tids.end(); t++) {
      out << " " << t->first<< " 1";
    }
    out << endl;
    if (ftell(in) == end) {
      break;
    }
  }
  fclose(in);
  out.close();
}
