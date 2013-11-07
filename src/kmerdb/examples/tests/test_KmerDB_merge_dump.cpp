#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << " output_fn fasta_1.fa [fasta_2.fa ...]\n";
    cerr << "function: loads a KmerDB with one or more files,\n";
    cerr << "          then dumps kmers to output.fn is ascii format\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv);
  }
  for (int j=0; j<argc; j++) {
    if (strcmp(argv[j], "-h") == 0) {
      usage(argv);
    }
  }


  char * gid_to_tid_fn = "../../../../../taxonomy/kpath/gid_to_kpath_tax_id.dat";

  TaxTree tax_tree("../../../../../taxonomy/kpath//kpath_taxonomy.dat", gid_to_tid_fn);
  GenomeIdToTaxId genome_to_taxid(gid_to_tid_fn);


  TaxNodeIDArray data;
  KmerNode::s_tax_tree = &tax_tree;
  KmerNode::s_tax_data = &data;

  string output_fn = argv[1];

  KmerDB kdb;
  for (int j=2; j<argc; j++) {
    cout << "loading: " << argv[j] << endl;
    kdb.ingest(argv[j], &genome_to_taxid, &data);
    KmerFileMetaData metadata;
    metadata.read(argv[j]);
    //metadata.write();
    //cout << endl;
  }

  cout << "\nwriting kmers to: " << output_fn << endl;
  kdb.dumpKmers(output_fn.c_str());
}
