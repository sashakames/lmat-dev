#include "KmerSingleBitHashTable.hpp"
#include "parse_dna.hpp"
#include <cassert>
#include <math.h>
#include <sstream>

using namespace metag;
using namespace std;


KmerSingleBitHashTable::~KmerSingleBitHashTable() { 
}


KmerSingleBitHashTable::KmerSingleBitHashTable(int k) :  m_setbit_count(0), m_set_fwd_count(0), m_set_rc_count(0), m_k(k) {

  //maximum possible number of kmers
  double m = pow(4.0, (double)(k));
  m_maxnumofkmers = (uint64_t) m;

  uint64_t tabsize = NUMOFINTSFORBITS(m_maxnumofkmers);
  cout << "tabsize: " << tabsize << " max kmer: " << m_maxnumofkmers << endl;

  m_hash_table.resize(tabsize);
  m_tab = const_cast<uint64_t*>(&m_hash_table[0]);
  memset(m_tab,0,sizeof (uint64_t) * m_hash_table.size());

  m_characters = "ACGT";
}

void KmerSingleBitHashTable::getCount(uint64_t &is_set, uint64_t &not_set) const {
  is_set = 0;
  not_set = 0;
  for (uint64_t code=0; code<m_maxnumofkmers; code++) {
    if (find(code)) {
      ++is_set;
    } else {
      ++not_set;
    }
  }
}


void KmerSingleBitHashTable::printTextKmers() const {
  char kmer[128];
  for (uint64_t code=0; code<m_maxnumofkmers; code++) {
    if (find(code)) {
      parse_dna::mer_binary_to_string(code, m_k, kmer);
      cout << kmer << endl;
    }
  }
}
#if 0
void KmerSingleBitHashTable::printTextKmers(const char *fn, bool print_kmers) const {
assert(false);
#if 0
  uint64_t ct = 0;
  stringstream ss;
  ss << fn << "_t";
  ofstream out(ss.str().c_str());
  assert(out.is_open());
  cerr << "opened for writing: " << ss.str() << endl;
  char kmer[32];

  if (m_rank == 0) {
    uint64_t s, n;
    getCount(s, n);
    cout << "set: " << s << " not: " << n << endl;
  }

  strncpy(kmer, m_prefix.c_str(), m_prefix.size());
  for (uint64_t code=0; code<m_maxnumofkmers; code++) {
    if (print_kmers) {
      //a set bit represents a kmer
      if (isBitSet(code)) {
        Coder::kmercode2string(kmer+m_prefix.size(), code, m_kk);
        out << kmer << endl;
        ++ct;
      }
    } else {
      //an unset but represents an unword
      if (!isBitSet(code)) {
        Coder::kmercode2string(kmer+m_prefix.size(), code, m_kk);
        out << kmer << endl;
        ++ct;
      }
    }
  }
  out.close();
  cout << "KmerSingleBitHashTable::printTextKmers, count: " << ct << endl;
#endif
}

void KmerSingleBitHashTable::printBinaryKmers(const char *fn, bool print_kmers) const {
  //for debugging
  uint64_t ct = 0;

  //open output file
  stringstream ss;
  ss << fn << "_b";
  FILE *out = fopen(ss.str().c_str(), "w");
  assert(out != NULL);

  static char kmer[32];
  kmer[m_k] = '\0';
  for (int j=0; j<m_k; j++) {
    kmer[j] = 'a';
  }

  uint64_t fw, rc;
  strncpy(kmer, m_prefix.c_str(), m_prefix.size());
  for (uint64_t code=0; code<m_maxnumofkmers; code++) {
    if (print_kmers) {
      if (isBitSet(code)) {
        Coder::decode(kmer+m_prefix.size(), code, m_kk);
        Coder::encode(kmer, fw, rc, m_k);
        fwrite(&fw, sizeof(uint64_t), 1, out);
        ++ct;
      }
    } else {
      if (! isBitSet(code)) {
        Coder::decode(kmer+m_prefix.size(), code, m_kk);
        Coder::encode(kmer, fw, rc, m_k);
        fwrite(&fw, sizeof(uint64_t), 1, out);
        ++ct;
      }
    }
  }
  fclose(out);
  cout << "KmerSingleBitHashTable::printBinaryKmers, count: " << ct << endl;
}



void KmerSingleBitHashTable::flipBits() {
  uint64_t x = 0;
  x = ~x;
  for (size_t j=0; j<m_hash_table.size(); j++) {
    m_hash_table[j] ^= x;
  }
}

void KmerSingleBitHashTable::printTable(const char *fn) const {
  FILE *fp = fopen(fn, "w");
  assert(fp != NULL);
  const char *t = (char*) &m_hash_table[0];
  uint64_t sz = m_hash_table.size()*8;
  uint64_t test = fwrite(t, sizeof(char), sz, fp);
  if (test != sz) {
    cout << "fwrite returned " << test << "  should be " << sz << endl;
    exit(-1);
  }
  assert(fclose(fp) == 0);
  cout << "\nKmerSingleBitHashTable::printTable wrote " << sz << " bytes\n";
}

void KmerSingleBitHashTable::loadTable(const char *fn) {
  uint64_t sz = Util::getFileLength(fn);
  cout << "KmerSingleBitHashTable::loadTable; sz= " << sz << endl;
  FILE *fp = fopen(fn, "r");
  assert(fp != NULL);
  m_hash_table.resize(sz/8);
  cout << "m_hash_table.resize(" << sz/8 << ")\n";
  char *t = (char*)&m_hash_table[0];
  uint64_t len = fread(t, sizeof(char), sz, fp);
  if (len != sz) {
    cout << "fread returned " << len << "  should be " << sz << endl;
  }
  assert(fclose(fp) == 0);
}

/*
void KmerSingleBitHashTable::checkKmerOcurrences(uint64_t &code)  {
    static char kmer[32];
    char *kmer_ptr = kmer + m_prefix.size();
    static char kmer2[32];
    static uint64_t ct = 0;
    ++ct;
    static uint64_t fwd, rc, cc, tmpcode;
    Coder::kmercode2string(kmer, code, m_k);
    for (size_t i=0; i<m_prefix.size(); i++) {
      if (m_prefix[i] != kmer[i]) {
        return;
      }
    }
    Coder::encode(kmer_ptr, fwd, rc, m_kk);
    setBit(fwd);
    setBit(rc);
  }
*/

#if 0
void KmerSingleBitHashTable::setBit(unsigned long code) {
    m_tab[DIVWORDSIZE(code)] |= ITHBIT(MODWORDSIZE(code));
  }

bool KmerSingleBitHashTable::isBitSet(unsigned long code) const {
    return m_tab[DIVWORDSIZE(code)] & ITHBIT(MODWORDSIZE(code));
}

#  define __INT64_C(c)  c ## L
# define INT64_MAX    (__INT64_C(9223372036854775807))


void KmerSingleBitHashTable::checkKmerOcurrences(uint64_t &code)  {
++total_codes;
static uint64_t ctr = 0;
++ctr;

static char kmer[33];
static string test("TTTTTTTTTTTTTTTTT");
Coder::decode(kmer, code, m_kk);
if (test == kmer) cout << endl << "got one: " << ctr << "  " << CM(DIVWORDSIZE(code)) << endl;


if (DIVWORDSIZE(code) > m_hash_table.size()) {
  cerr << "error; code= " << Util::commify(code) << " DIVWORDSIZE(code): " << Util::commify((uint64_t)DIVWORDSIZE(code)) << " table size: " << Util::commify(m_hash_table.size()) << " max uint64_t: " << Util::commify(INT64_MAX) << " m_k: " << m_k << " m_kk: " << m_kk << endl;
  cerr << "kmer number: " << ctr << endl;
  Coder::decode(kmer, code, m_kk);
  cerr << "decoded kmer: " << kmer << endl;
  ++bad_codes;
  return;
  //exit(0);
}
//cout << code << " table sz:" << m_hash_table.size() << endl;
    if (m_k == m_kk) {
      if (!isBitSet(code)) {
        setBit(code);
      }
    } else {
      static uint64_t fwd, rc, cc, tmpcode;
      Coder::kmercode2string(kmer, code, m_k);
      for (size_t i=0; i<m_prefix.size(); i++) {
        if (m_prefix[i] != kmer[i]) {
          return;
        }
      }
      Coder::encode(kmer+m_prefix.size(), fwd, rc, m_kk);
      if (!isBitSet(fwd)) {
        ++m_setbit_count;
        setBit(fwd);
      }
    }
  }
#endif
/*
void KmerSingleBitHashTable::printFastaFormatedKmers(const char) const {
  uint64_t ct = 0;
  stringstream ss;
  ss << fn << "_t";
  ofstream out(ss.str().c_str());
  assert(out.is_open());
  cerr << "opened for writing: " << ss.str() << endl;
  char kmer[32];

  if (m_rank == 0) {
    uint64_t s, n;
    getCount(s, n);
    cout << "set: " << s << " not: " << n << endl;
  }
  out << ">diff\n";

  strncpy(kmer, m_prefix.c_str(), m_prefix.size());
  for (uint64_t code=0; code<m_maxnumofkmers; code++) {
      //a set bit represents a kmer
      if (isBitSet(code)) {
        Coder::kmercode2string(kmer+m_prefix.size(), code, m_kk);
        out << kmer << "N";
        ++ct;
      }
  }
  out.close();
  cout << "KmerSingleBitHashTable::printTextKmers, count: " << ct << endl;
}
*/
#endif
