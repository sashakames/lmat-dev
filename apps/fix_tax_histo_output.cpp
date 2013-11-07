#include <iostream>
#include <string>
#include "all_headers.hpp"

using namespace std;
using namespace metag;

void usage(char *argv[]) {
    cerr << "usage: " << argv[0] << " tax_histo_fn output_fn\n";
    cerr << "function: detects and corrects duplicate taxIDs\n";
    exit(-1);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    usage(argv);
    exit(-1);
  }

  const char *fn = argv[1];
  ifstream in(fn);
  assert(in);
  Utils::skipHeaderLines(in);
  ofstream out(argv[2]);
  assert(out);

  uint32_t tid, count, lca_count, tid_count, genome_count;
  uint64_t kmer, x = ~0;
  map<int, int> m; 
  uint64_t z = -1;
  uint64_t ct = 0;
  while (true) {
    m.clear();
    ++z;
    kmer = ~0;
    in >> kmer;
    if (kmer == x) {
      break;
    }  
    in >> genome_count;
    in >> tid_count;
    in >> lca_count;
    for (uint32_t j=0; j<lca_count; j++) {
      in >> tid;   
      in >> count;
      if (m.find(tid) != m.end()) {
        if (count > (int)m[tid]) {
          m[tid] = count;
          ++ct;
        }
      } else {
        m[tid] = count;
      }  
    }
    out << kmer << " " << genome_count << " " << tid_count << " " << m.size() << " ";
    for (map<int,int>::const_iterator t = m.begin(); t != m.end(); t++) {
      out << t->first << " " << t->second << " ";
    }
    out << endl;
  }
  in.close();
  out.close();
  cout << "processed " << z << " entries\n";
  cout << "num duplicates deleted: " << ct << endl;
}
