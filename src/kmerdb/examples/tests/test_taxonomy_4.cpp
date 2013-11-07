#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << " kmer_data_fn" << endl;
    cerr << "function: tests that genome IDs map to mfe IDs\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
    }
  }
  if (argc != 2) {
      usage(argv);
  }    

  GenomeIdToTaxId genome_to_taxid("../../../../../taxonomy/kpath/gid_to_kpath_tax_id.dat");

/*
  TaxNodeIDArray data;
  KmerDB db;
  GenomeIdToTaxId lookup("../../../../../taxonomy/kpath/gid_to_kpath_tax_id.dat");

  db.ingest(argv[1], &lookup, &data);
  //exit(0);
*/

  //read offsets file into memory
  char buf[1024];
  sprintf(buf, "%s.offsets", argv[1]);
  uint64_t *offsets = (uint64_t*)Utils::readFile(buf);

  //read kmer data file into memory
  char *kmer_data = Utils::readFile(argv[1]);

  //get the kmer length and number of kmers in the file
  KmerFileMetaData metadata;
  metadata.read(argv[1]);
  uint64_t kmer_count = metadata.size();
  metadata.write();

  set<uint32_t> found_gids; 
  set<uint32_t> not_found_gids; 
  uint32_t gid;
  uint32_t gid_count;

  uint32_t *p;
  uint64_t kmer;
  char *b;
  uint32_t g_ct;
  for (uint64_t j=0; j<kmer_count; j++) {
  b = kmer_data + offsets[j];
  kmer = *(kmer_t*)b;
  b += sizeof(kmer_t);
  uint32_t *bb = (uint32_t*)b;
  g_ct = bb[0];
  ++bb;


  GidHash::const_iterator iter;

  //loop over genome IDs, and insert their associated taxid into s_work
  for (uint32_t h=0; h<g_ct; h++) {
    gid = bb[0];
    ++bb;
      iter = genome_to_taxid.find(gid);
      //if (genome_to_taxid.find(gid) != genome_to_taxid.end()) {
      if (iter != genome_to_taxid.end()) {
         found_gids.insert(gid);
       } else {
         not_found_gids.insert(gid);
         //std::cout << "failed to find tax ID for genome with kpath ID: " << gid << " " << __FILE__ << " " << __LINE__ <<endl;
    }
  }
  /*
     c = kmer_data + offsets[j] + 8; //+8 to skip over kmer
     p = (uint32_t*)c;
     gid_count = p[0];
     p += 4;
     cout << "count: " << gid_count << endl;
     for (uint32_t h=0; h<gid_count; h++) {
       gid = p[0];
       p += 4;
       cout << "  gid: " << gid<< endl;
       if (genome_to_taxid.find(gid) == genome_to_taxid.end()) {
         not_found_gids.insert(gid);
       } else {
         found_gids.insert(gid);
       }
     }
     */
  }

  cout << "mapped " << found_gids.size() << " kpath IDs to mfe IDs\n";
  cout << "failed to map " << not_found_gids.size() << " kpath IDs to mfe IDs\n";

  for (set<uint32_t>::iterator t = not_found_gids.begin(); t != not_found_gids.end(); t++) {
         cout << "not found: " << *t<< endl;
  }
  cout << endl;
  for (set<uint32_t>::iterator t = found_gids.begin(); t != found_gids.end(); t++) {
         cout << "found: " << *t<< endl;
  }

  return 0;
}
