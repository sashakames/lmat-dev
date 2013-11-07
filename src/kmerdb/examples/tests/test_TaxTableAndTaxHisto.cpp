#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " tax_hist_output_fn\n";
    exit(-1);
  }


#if USE_GENOME_LOCATIONS != 0
#error "Must compile with USE_GENOME_LOCATIONS=0"
#endif

  TaxTable tt;
  tt.registerFile(argv[1]);
  StopWatch c;
  c.start();
  tt.ingest(false);
  double d = c.stop();
  cout << "kmers/s: " << 100000/d << endl;

  vector<uint64_t> kmers_to_test;
  kmers_to_test.push_back(0);
  kmers_to_test.push_back(535);

  for (int j=0; j<kmers_to_test.size(); j++) {
    TaxNodeStat ts(tt);
    ts.begin(kmers_to_test[j]);
    while (ts.next()) {
      uint32_t taxid = ts.taxid();
      uint16_t present = ts.present();
      cout << "("<<taxid<<","<<present<<") ";
    }
  }
}
