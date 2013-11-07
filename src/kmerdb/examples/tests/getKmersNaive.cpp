#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <cassert>
#include <cstdlib>

using namespace std;

void RC(string &in, string &out) {
  char c;
  size_t k = 0;
  for (int j=in.size()-1; j>=0; j--, k++) {
    switch (in[j]) {
      case 'A' : c = 'T'; break;
      case 'C' : c = 'G'; break;
      case 'T' : c = 'A'; break;
      case 'G' : c = 'C'; break;
      default : cerr << "ERROR!: input string: " << in << " j: " << j << " k: " << k << " in[j]: " << in[j] << endl;
      assert(false);
    }
    out[k] = c;
  }
}


void getKmers(string &line, set<string> &kmers, int K) {
  string s;
  string rc;
  rc.resize(K);
  for (size_t j=0; j<line.size()-K; j++) {
    s = line.substr(j, K);
    bool good = true;
    for (size_t h=0; h<s.size(); h++) {
      char c = s[h];
      if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
        good = false;
        break;
      }
    }
    if (good) {
      RC(s, rc);
      if (s < rc) {
        kmers.insert(s);
      } else {
        kmers.insert(rc);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " fasta_fn K\n";
    return (-1);
  }

  set<string> kmers;
  int K = atoi(argv[2]);
  string line;
  ifstream in(argv[1]);
  assert(in.is_open());
  while (true) {
    getline(in, line);
    if (line.size() == 0) {
      break;
    }
    getline(in, line);
    getKmers(line, kmers, K);
  }

  for (set<string>::const_iterator t = kmers.begin(); t != kmers.end(); t++) {
    cout << *t << endl;
  }
}
