#include "../../all_headers.hpp"
#include <iostream>
#include <string>

using namespace std;
using namespace metag;
using namespace __gnu_cxx;

void usage(char *argv[]) {
    exit(-1);
}

int main(int argc, char *argv[]) {

  //load data from tax_histo into a TaxTable
  TaxTable tt;
  tt.registerFile("/usr/mic/post1/metagenomics/tmp_kmer_dbs/all_virus.int.bin.no_loc.16.partition_0.db.r.bin");;
  tt.ingest();

  //we use a TaxNodeStat to access that statistical information
  TaxNodeStat tax_stat(tt);
  
  vector<uint64_t> kmers;
  kmers.push_back(3630760592);
  kmers.push_back(2413096856);

  cout << "correct answer for first kmer:\n3630760592 4 0 7 493171 4 128 8859 1 0 32820 2 1 287054 1 87 500241 1 0 20001253 1 0 121843 1 0\n";
  cout << endl << "correct answer for second kmer:\n2413096856 2 0 3 167365,2,172 485590,1,0 485591,1,0\n";
  cout << endl << "==================================================================================\n";


  //loop over the test kmers
  for (size_t j=0; j<kmers.size(); j++) {
    //access the data: call begin(), then loop over next()
    tax_stat.begin(kmers[j]);
    cout << kmers[j] << " " <<tax_stat.taxidCount()<<" "<<tax_stat.genomeCount() << " " << tax_stat.tupleCount() << " ";

    while (tax_stat.next()) {
      cout << tax_stat.taxid()<< " " << tax_stat.present() << " " << tax_stat.absent() << " ";
    }
    cout << endl << endl;
  }  
}
