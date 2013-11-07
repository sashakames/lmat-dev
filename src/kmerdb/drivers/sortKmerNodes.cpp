#include "../all_headers.hpp"
#include "KmerNodeSorter_no_locations.hpp"
#include <iostream>
#include <vector>

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " input_fn\n";
    cerr << "where: input_fn contains output from jellylist.bin\n";
    cerr << "output file is: input_fn.sorted\n";
    exit(-1);
  }

  StopWatch clock;

  //get the kmer length and number of kmers in the file
  KmerFileMetaData metadata;
  metadata.read(argv[1]);

  //open output files
  char buf[1024];
  sprintf(buf, "%s.sorted.offsets", argv[1]);
  FILE *out_offsets = fopen(buf, "w");
  assert(out_offsets != NULL);

  sprintf(buf, "%s.sorted", argv[1]);
  FILE *out = fopen(buf, "w");
  assert(out != NULL);

  //read the data file to buffer
  clock.reset(); clock.start();
  sprintf(buf, "%s", argv[1]);
  cout << "calling Utils::readFile( " << buf << " )\n";
  char * data = Utils::readFile(buf);
  const  char *bb = data;
  cout << "time to read data file: " << clock.stop() << endl;

  //read the offsets file
  sprintf(buf, "%s.offsets", argv[1]);
  clock.reset(); clock.start();
  const uint64_t *offsets = (uint64_t*)Utils::readFile(buf);
  cout << "time to read offsets file: " << clock.stop() << endl;

  //allocate storage for KmerNodes
  clock.reset(); clock.start();
  KmerNodeSorter_no_locations *nds = new KmerNodeSorter_no_locations[metadata.size()];
  cout << "time to allocate storge for kmer nodes: " << clock.stop() << endl;

  vector<KmerNodeSorter_no_locations*> v;
  v.reserve(metadata.size());

  clock.reset(); clock.stop();
  KmerNodeSorter_no_locations *w;
  uint64_t kmer_count = metadata.size();
  for (uint64_t j=0; j<kmer_count; j++) {
    w = &nds[j];
    w->read(bb+offsets[j]);
    v.push_back(w);
  }
  cout << "time to read kmers: " << clock.stop() << endl;

  clock.reset(); clock.start();
  sort(v.begin(), v.end(), KmerNodeSorter_no_locations::NodeCmp());
  cout << "time to sort kmers: " << clock.stop() << endl;

  metadata.write(out);

  uint64_t offset;
  uint64_t sanity = ~0;
  cout << "kmer count: " << kmer_count << endl;
  cout << "reading kmers beginning at offset " << ftell(out) << endl;
  for (uint64_t j=0; j<kmer_count; j++) {
    //write the offset where data for the KmerNode appears in the 
    //kmer data file
    offset = ftell(out);
    assert(fwrite(&offset, sizeof(uint64_t), 1, out_offsets) == 1);

    //write the kmer
    v[j]->writeBinary(out);
    
    //write sanity marks
    if ((j+1) % 1000 == 0) {
      assert(fwrite(&sanity, sizeof(uint64_t), 1, out) == 1);
    }
  }

  fclose(out);
  fclose(out_offsets);
  delete [] offsets;
  delete [] nds;
}
