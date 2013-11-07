#include "GiToTid.hpp"
#include <iostream>
#include <fstream>
#include <cassert>

using namespace metag;
using namespace std;

GiToTid::GiToTid(const char *fn) {
  ifstream in(fn);
  if (!in.is_open()) {
    cerr << "failed to open " << fn << " for reading\n";
    exit(-1);
  }
  int gid = -1, tid = -1;
  while (in.good()) {
    in >> gid;
    in >> tid;
//cout << "adding: " << gid << " " << tid << endl;
    if (gid == -1) {
      break;
    }
    m_map[gid] = tid;
  }
  in.close();
}

int GiToTid::getTid(int genome_id) {
  if (m_map.find(genome_id) == m_map.end()) {
    return -1;
  }
  return m_map[genome_id];
}

