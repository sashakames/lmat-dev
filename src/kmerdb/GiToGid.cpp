#include "GiToGid.hpp"
#include <iostream>
#include <fstream>
#include <cassert>

using namespace metag;
using namespace std;

GiToGid::GiToGid(const char *fn) {
  ifstream in(fn);
  if (!in.is_open()) {
    cerr << "failed to open " << fn << " for reading\n";
    exit(-1);
  }
  int gi = -1, gid = -1;
  while (in.good()) {
    in >> gi;
    in >> gid;
//cout << "adding: " << gid << " " << tid << endl;
    if (gid == -1) {
      break;
    }
    m_map[gi] = gid;
  }
  in.close();
}

int GiToGid::getGid(int gi) {
  if (m_map.find(gi) == m_map.end()) {
    return -1;
  }
  return m_map[gi];
}

