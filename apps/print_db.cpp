#include <unistd.h>
#include <fstream>
#include <math.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <omp.h>
#include "kencode.hpp"
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"

#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

using namespace kencode_ns;
using namespace metag;

void proc_kmer(TaxTable *table, uint64_t kmer_id)
{
  TaxNodeStat *h = new TaxNodeStat(*table);

  h->begin(kmer_id);
  cout<<"kmer "<<kmer_id<<" ";
  const uint16_t ng = h->genomeCount();
  if ( ng > 0 ) {
    cout<<" genome count "<<ng<<" leaf taxid count "<<h->taxidCount();
    while( h->next() ) {  
      const uint32_t tid = h->taxid();
      //collect the unique set of tax ids
      const uint16_t pr_cnt = h->present();
      cout<<" taxid  "<<tid<<" pres-cnt "<<pr_cnt;
    }
  } 
  else {
    cout<<" none";
  }
  cout<<endl; 
  delete h;
}

int main(int argc, char* argv[]) 
{
  char c = '\0';

  omp_set_num_threads(1);

  int max_count = 50;
  int count = 0;
  int k_size = 0;

  string fn = "/mnt/virident/maya/marker_lib.20.v1.1";
  while ((c = getopt(argc, argv, "f:n:")) != -1) {
    switch (c) {
    case 'f':
      fn = optarg;
      break;
    case 'n':
      max_count = atoi(optarg);
      break;
    default:
      cout << "Unrecognized option: "<<c<<", ignore."<<endl;
    }
  }
  cout << "Start kmer DB load\n";

  INDEXDB *taxtable;

  cout << "starting restore for " << fn << endl;

  perm(&taxtable, sizeof(taxtable));
  mopen(fn.c_str(), "r", 0);

  if (k_size < 1)
    k_size = taxtable->get_kmer_length();

  cout << "There are " << taxtable->size() << " " << k_size << "-mers in the database."  <<   endl ;

  INDEXDB::iterator hm_iter;

  for (hm_iter = taxtable->begin(); hm_iter != taxtable->end(); hm_iter++)
    {
      proc_kmer(taxtable, hm_iter->first);
      count++;
      if (count >=max_count) break;
    }
  return 0;

}
