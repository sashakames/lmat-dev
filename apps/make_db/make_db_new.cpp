/**
Create a mmap kmer db file of the desired data structure
*/
#include "TwoLevelRW.hpp"


#include <fstream>
#include <omp.h>

using namespace std;
using namespace metag;


//typedef pair<size_t, int> kmer_info_t; 

void usage() {
  cout<<"Usage: \n"
    //      "  -f 16|32  - determines if 16 or 32 bit tax IDS are used [required]\n"
    "  -i <fn>   - input file -or- filename of file that contains a listing of the input binary files\n"
    
    "              output from tax_histo_fast_limited [required]\n"
    "  -l        - is file a list?        [default=no, optional]\n"       
    "  -o <fn>   - output filename                    [required]\n"
    "  -k <int>  - k-mer length used in input file    [required]\n"
    "  -h        - if given, assumes input files are from kmerPrefixCounter, instead of tax_histo [optional]\n"
    "  -s <int>  - size to reserve for the memory-mapped file [optional; default is -s 500, or 500G]\n"
    "  -q <int>  - stop at N k-mers (per input file) - helpful for troubleshooting \n"
    "  -g <int>  - taxid list cutoff \n"
    "  -m <fn>   - pruning taxid table file  \n"    
    "  -f <fn>   - use 32 to 16-it id conversion using <fn>  \n"    
    "  -w        - pruing uses strain-to-species mapping\n" 
    "guidance for setting -s: from our paper, our full reference DB required 619G;\n"
    "this was for a fasta file that was ~19G\n";
}

 
int main(int argc, char *argv[]) {

  string inputfn, outputfn;

  int c, n_threads;


  bool restore = false;

  // default database size in GB
  size_t mmap_size = 255;

  bool list = false;

  int kmer_len = 0;

  size_t hash_size = 1;
  size_t storage_size = 1;

  unsigned long long stopper = 0;
  n_threads = -1;

  int count = 0;  

  int tid_cut = 0;

  string species_map_fn, id_bit_conv_fn;

  bool strainspecies = false;
  //  while ((c = getopt(argc, argv, "t:g:q:k:i:o:s:t:h:r l a ")) != -1) {

  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;

  
  uint16_t in_taxid ;
 
 while ((c = getopt(argc, argv, "g:q:k:i:o:s: l h m:f:w ")) != -1) {
    switch(c) {
    case 'd':
      in_taxid = atoi(optarg);
      break;
      
    case 'f':
      id_bit_conv_fn = optarg;
      break;
    case 'q':
      stopper = strtoull(optarg, NULL, 10);
      break;
      /*    case 'h':
      hash_size = strtoull(optarg, NULL, 10);
      break; */
    case 'k':
      ++count;
      kmer_len = atoi(optarg);
      break;
    case 'l':
      list = true;
      break;
   case 'r':
      restore = true;
      break;
    case 'i':
      ++count;
      inputfn = optarg;
      break;
    case 'o':
      ++count;
      outputfn = optarg;
      break;
    case 't':
      n_threads = atoi(optarg);
       omp_set_num_threads(n_threads);
       break;
    case 's':
      mmap_size = atoi(optarg);
      break;
    default:
      usage();
      exit(1);
    }
  }

  

  // setup kmerdb
   TwoLevelRW *ttable;
  


  mmap_size = mmap_size * (1<<30);

  cout << "size requested: " << mmap_size << endl;


  perm(&ttable, sizeof(ttable));

  if (!restore) {
    mopen(outputfn.c_str(), "w+", mmap_size);
    ttable = PERM_NEW(TwoLevelRW)(kmer_len);
  }
  else
    mopen(outputfn.c_str(), "r+", 0);


  ttable->set_base_addr(reinterpret_cast<char*>(pje_get_base()));
  ttable->addData(inputfn.c_str(), in_taxid);


  return(0);
}
