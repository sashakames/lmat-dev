#include "all_headers.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

using namespace std;
using namespace metag;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains tax_histo binary output\n";
    exit(-1);
  }

//  FILE *fp = Utils::openReadFile(argv[1]);
  ifstream fp(argv[1], ios::in | ios::binary);
  /*
  fseek(fp, 0, SEEK_END);
  long end = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  */

  uint64_t kmer;
  uint16_t genome_count, tid_count;
  map<uint32_t, uint16_t> tids;
  while (true) {
  /*
    if (ftell(fp) == end) {
      break;
    }
    */
    if (fp.eof()) break;
    Utils::readTaxHisto_bin(kmer, genome_count, tid_count, tids, fp);
    cout << kmer<<" "<< " g_ct: " << genome_count << " leaf_ct: " << tid_count << " :: tids: ";
    for (map<uint32_t, uint16_t>::const_iterator t = tids.begin(); t != tids.end(); t++) {
      cout << t->first<< " ";
    }
    cout << endl;
  }

 // fclose(fp);
}
