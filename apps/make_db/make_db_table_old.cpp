/**
Create a mmap kmer db file of the desired data structure
*/
#include <unistd.h>

#include "all_headers.hpp"
#include <TaxTable.hpp>

#define DATA_DIR "../../src/kmerdb/examples/tests/data/"

#include <fstream>

using namespace std;
using namespace metag;

void usage();

#if 1
size_t perm_bytes_allocd = 0;
#endif

int main(int argc, char *argv[]) {

  string inputfn, outputfn;

  int c;

  inputfn = DATA_DIR "test.fa.int.bin.no_loc.12";
  outputfn = DATA_DIR "test.fa.mmap";

  size_t mmap_size = 500;
  bool list = true;
  int kmer_len = 0;
  size_t hash_size = 5000000000;
  int ct = 0;
  int which = 0;
  bool use_tax_histo_files = true;

  while ((c = getopt(argc, argv, "f:k:i:o:s:h:g ")) != -1) {
    switch (c) {
    case 'g':
      use_tax_histo_files = false;
      break;
    case 'f':
      ++ct;
      which  = atoi(optarg);
      break;
    case 'h':
      hash_size = strtoull(optarg, NULL, 10);
      break;
    case 'i':
      ++ct;
      inputfn = optarg;
      break;
    case 'o':
      ++ct;
      outputfn = optarg;
      break;
    case 's':
      mmap_size = atoi(optarg);
      break;
    case 'k':
      ++ct;
      kmer_len = atoi(optarg);
      break;
    default:
      usage();
      exit(1);
    }
  }

  if (ct != 4 || !(which == 16 || which == 32)) {
    usage();
    exit(1);
  }

  cout << "which: " << sizeof(which) << endl;

  // setup kmerdb
  TaxTable<uint16_t> *ttable_16;
  TaxTable<uint32_t> *ttable_32;

  mmap_size = mmap_size * (1<<30);
  cout << "size requested: " << mmap_size << endl;

#if USE_BOOST == 1
  bip::managed_mapped_file mfile(bip::create_only, (const char*) outputfn.c_str(
                                 ) , mmap_size );
  if (which == 16) {
    ttable_16 = mfile.construct<TaxTable<uint16_t> >("KmerDB")(mfile, hash_size);
  } else {
    ttable_32 = mfile.construct<TaxTable<uint32_t> >("KmerDB")(mfile, hash_size);
  }

#else

  if (which == 16) {
    perm(&ttable_16, sizeof(ttable_16));
    mopen(outputfn.c_str(), "w+", mmap_size);
    ttable_16 = PERM_NEW(TaxTable<uint16_t>)(hash_size);
  } else {
    perm(&ttable_32, sizeof(ttable_32));
    mopen(outputfn.c_str(), "w+", mmap_size);
    ttable_32 = PERM_NEW(TaxTable<uint32_t>)(hash_size);
  }


#endif


  if (list) {
    ifstream ifs(inputfn.c_str());
    assert(ifs);
    string line;
    while (ifs>>line) {
      if (which == 16) {
        ttable_16->registerFile(line.c_str());
        ttable_16->set_kmer_length(kmer_len);
      } else {
        ttable_32->registerFile(line.c_str());
        ttable_32->set_kmer_length(kmer_len);
      }
    }
    if (which == 16) {
      ttable_16->ingest(use_tax_histo_files);
      cout << "DB size is " << ttable_16->size() << endl;
    } else {
      ttable_32->ingest(use_tax_histo_files);
      cout << "DB size is " << ttable_32->size() << endl;
    }
  } else {
    assert(false);
    /*
    cout << "input fn: " << inputfn << endl;
    if (which == 16) {
      ttable_16->registerFile(inputfn.c_str());
      ttable_16->set_kmer_length(kmer_len);
      ttable_16->ingest();
      cout << "DB size is " << ttable_16->size() << endl;
    } else {
      ttable_32->registerFile(inputfn.c_str());
      ttable_32->set_kmer_length(kmer_len);
      ttable_32->ingest();
      cout << "DB size is " << ttable_32->size() << endl;
    }
    */
  }

  //mclose();

#if 1
  cout << "Perm hash alloc size: " << perm_bytes_allocd << endl;
#endif

#if USE_BOOST != 1
  mclose();
#endif

  return(0);
}

void usage() {
  cout<<"Usage: \n"
      "  -f 16|32  - determines if 16 or 32 bit tax IDS are used [required]\n"
      "  -i <fn>   - filename of file that contains a listing of the input binary files\n"
      "              output from tax_histo_fast_limited [required]\n"
      "  -o <fn>   - output filename                    [required]\n"
      "  -k <int>  - k-mer length used in input file    [required]\n"
      "  -g        - if given, assumes input files are from kmerPrefixCounter, instead of tax_histo [optional]\n"
      "  -s <int>  - size to reserve for the memory-mapped file [optional; default is -s 500, or 500G]\n"
      "  -h <int>  - expected number of input k-mers [optional; default is -h 5000000000]\n"
      "\n"
      "guidance for setting -s: from our paper, our full reference DB required 619G;\n"
      "this was for a fasta file that was ~19G\n";
}
