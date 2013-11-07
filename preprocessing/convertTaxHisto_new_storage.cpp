#include "mpi.h"
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include "all_headers.hpp"

using namespace metag;
using namespace std;

void doit(const char *input_dir, const char *fn, const char *output_dir); 

int main(int argc, char* argv[]) {
  if (argc != 5) {
    cerr << "usage: " << argv[0] << " input_dir base_fn count output_dir\n";
    exit(-1);
  }

  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const char *input_dir = argv[1];;
  const char *base = argv[2];
  int count = atoi(argv[3]);
  const char *output_dir = argv[4];

  char buf[1024];
  for (int j=0; j<count; j++) {
    if (j % size == rank) {
      sprintf(buf, "%s%d", base,j);
      StopWatch clck;
      clck.start();
      doit(input_dir, buf, output_dir);
      if (rank == 0) {
        cout << "time to process file: " << clck.stop() << endl;
      }
    }
  }

  MPI_Finalize();
}

void doit(const char *input_dir, const char *fn, const char *output_dir) {
  //open output files
  char buf[1024];
  sprintf(buf, "%s/%s.kmers", output_dir, fn);
  FILE *out_kmers = fopen(buf, "w");
  assert(out_kmers);
  sprintf(buf, "%s/%s.tids", output_dir, fn);
  FILE *out_tids = fopen(buf, "w");
  assert(out_tids);

  string line;
  uint64_t kmer;
  uint16_t genome_count, tid_count, present, tuple_count;
  uint32_t tid;
  uint64_t f;

  sprintf(buf, "%s/%s", input_dir, fn);
  ifstream in(buf);
  assert(in);

  //skip header lines
  Utils::skipHeaderLines(in);

  while (true) {
    //loop exit condition
    if (in.eof()) break;

    //read kmer, taxid count, and tuple count
    in >> kmer;
    in >> genome_count;
    in >> tid_count;
    in >> tuple_count;
    assert(fwrite(&kmer,8, 1, out_kmers) == 1);
    f = ftell(out_tids);
    assert(fwrite(&f,8, 1, out_kmers) == 1);

    for (uint16_t k=0; k<tuple_count; k++) {
      in >> tid;
      in >> present;
      assert(fwrite(&tid,4, 1, out_tids) == 1);
    }
  }
  fclose(out_tids);
  fclose(out_kmers);
  in.close();
}
