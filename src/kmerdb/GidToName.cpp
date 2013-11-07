#include "GidToName.hpp"
#include <fstream>
#include <cassert>

using namespace metag;
using namespace std;

GidToName::GidToName(const char *fn) {
  ifstream in(fn);
  assert(in.is_open());
  int gid;
  string name;
  while (in.good()) {
    in >> gid;
    in >> name;
    m_map[gid] = name;
  }
  in.close();
}

string not_found;
const string & GidToName::getName(int genome_id) {
  if (m_map.find(genome_id) == m_map.end()) {
    return not_found;
  }
  return m_map[genome_id];
}

