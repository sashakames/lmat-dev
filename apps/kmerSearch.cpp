/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory.
 * Written by LMAT development team. LLNL-CODE-605872 All rights reserved.
 * This file is part of LMAT Please read COPYRIGHT file in root directory.
 *
 * */

#include <iostream>
#include <cassert>
#include <stdlib.h>
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"


#define MMAP_SIZE 0

#define TID_T uint32_t

using namespace metag;

void usage(char *execname) {
  cout << "Usage:\n"
       << "  -d <string>  - database filename       [required]\n"
       << "  -i <string>  - query filename          [required]\n";
}

int main(int argc, char* argv[]) {
  string kmer_db_fn, query_fn;
  int count = 0;
  char c = '\0';
  while ((c = getopt(argc, argv, "i:d:")) != -1) {
    switch (c) {
    case 'i':
      ++count;
      query_fn = optarg;
      break;
    case 'd':
      ++count;
      kmer_db_fn = optarg;
      break;
    default:
      cout << "Unrecognized option: "<<c<<", ignore."<<endl;
    }
  }

  if (count != 2) {
    usage(argv[0]);
    return -1;
  }

  INDEXDB<TID_T> *taxtable;
  cout << "starting restore for " << kmer_db_fn << endl;

#if USE_BOOST == 1
  bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
  std::pair<INDEXDB<TID_T>*, std::size_t> ret = mfile.find<INDEXDB<TID_T> >("KmerDB");
  taxtable = ret.first;
  taxtable->conv_ptrs();
#else
  assert(false);
#endif

  cout << "num kmers: " << taxtable->size() << endl;

  //read kmers, and query against the taxtable
  FILE *f = fopen(query_fn.c_str(), "r");
  assert(f);
  uint64_t kmer;
  uint64_t kmers_read = 0;
  uint64_t kmers_found = 0;
  while (true) {
    if (fread(&kmer,8,1,f) != 1) {
      break;
    }
    ++kmers_read;
    if (taxtable->exists(kmer)) {
      ++kmers_found;
    }
  }
  cout << "kmers read:  " << kmers_read << endl;
  cout << "kmers found: " << kmers_found << endl;

  return 0;
}
