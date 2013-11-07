#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <string>
#include <cstdlib>
#include "mpi.h"

using namespace std;

char *base = "/g/g10/hysom/metag_repo/dev/apps/";

int main(int argc, char **argv) {
  int rank, size;

  if (argc != 2) {
    cout << "usage: " << argv[0] << " input\n"
         << "where: input contains a list of output files from jellylist\n";
    exit(1);     
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ifstream in(argv[1]);
  assert(in);


  int j = -1;
  string fn;
  while (in >> fn) {
    ++j;
    if (j%size == rank) {
      if (rank == 2) cout << j << " " << fn << endl;
    int x = fn.find(".bin");
    string table = fn.substr(0, x+4);
    table += ".table";

    stringstream s;
    cout << "1\n";

    char buf[1024];
    sprintf(buf, "%s.tax_histo", fn.c_str());
    s << "tax_histo_fast_limited -o " << buf <<  " -d " << fn
      << " -a " << table  << " -g " << base << "../runtime_inputs/gid_to_tid_all.txt -p 5 -t "
      << base << "../runtime_inputs/kpath_taxonomy.dat";
    cout << s << endl;
    system(s.str().c_str());
  cout << s.str() << endl;
  }
    }

  MPI_Finalize();
}
