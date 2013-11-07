#include <iostream>
#include <stdint.h>
#include "all_headers.hpp"
#include <cassert>
#include <string>
using namespace std;
using namespace metag;

int main(int argc, char **argv) {
  string a("aacc");
  string a_rc("ggtt");
  string b("ggtta");
  string b_rc("taacc");
  string c("ccaat");
  string c_rc("attgg");

  //test encode, decode methods
  uint64_t kmer;
  string mer;
  kmer = Encoder::encode(a);
  Encoder::decode(kmer, (int)a.size(), mer);
  assert(mer == a);

  kmer = Encoder::encode(b);
  Encoder::decode(kmer, (int)b.size(), mer);
  assert(mer == b);

  kmer = Encoder::encode(c);
  Encoder::decode(kmer, (int)c.size(), mer);
  assert(mer == c);

  //test rc (reverse complement)
  uint64_t rc;
  kmer = Encoder::encode(a);
  rc = Encoder::rc(kmer, (int)a.size());
  cout << ">>> kmer: " << kmer << " rc: " << rc << endl;
  Encoder::decode(rc, (int)a.size(), mer);
  cout << "a:    "<<a<<endl;
  cout << "a_rc: "<<a_rc<<endl;
  cout << "mer : "<<mer<<endl;
  assert(mer == a_rc);

  kmer = Encoder::encode(b);
  rc = Encoder::rc(kmer, (int)b.size());
  cout << ">>> kmer: " << kmer << " rc: " << rc << endl;
  Encoder::decode(rc, (int)b.size(), mer);
  cout << "b:    "<<b<<endl;
  cout << "b_rc: "<<b_rc<<endl;
  cout << "mer : "<<mer<<endl;
  assert(mer == b_rc);

  kmer = Encoder::encode(c);
  rc = Encoder::rc(kmer, (int)c.size());
  cout << ">>> kmer: " << kmer << " rc: " << rc << endl;
  Encoder::decode(rc, (int)c.size(), mer);
  cout << "c:    "<<c<<endl;
  cout << "c_rc: "<<c_rc<<endl;
  cout << "mer : "<<mer<<endl;
  assert(mer == c_rc);
}
