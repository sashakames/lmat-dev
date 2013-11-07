#ifndef __UNWORDS_INFO__
#define __UNWORDS_INFO__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>

//for uint64_t
#define LOGWORDSIZE 6

#define INTWORDSIZE (1U << LOGWORDSIZE)
#define DIVWORDSIZE(I) ((I) >> LOGWORDSIZE)
#define MODWORDSIZE(I) ((I) & (INTWORDSIZE-1))
#define NUMOFINTSFORBITS(N)\
        ((DIVWORDSIZE(N) == 0)\
           ? 1UL \
           : (1UL + DIVWORDSIZE((N) - 1UL)))

#define FIRSTBIT  (1UL << (INTWORDSIZE-1))

#define ITHBIT(I) (FIRSTBIT >> (I))

#define ISIBITSET(TAB,I)  ((TAB)[DIVWORDSIZE(I)] & ITHBIT(MODWORDSIZE(I)))

//#define SETIBIT(TAB,I)    (TAB)[DIVWORDSIZE(I)] |= ITHBIT(MODWORDSIZE(I))
#define SETIBIT(TAB,I,OFFSET)    (TAB)[DIVWORDSIZE(I)-OFFSET] |= ITHBIT(MODWORDSIZE(I))

/*
#define ALPHASIZE    4U

#define DIV2(N)      ((N) >> 1)
#define DIV4(N)      ((N) >> 2)
#define MULT4(N)     ((N) << 2)
*/

namespace metag {

class KmerSingleBitHashTable {
/**
Implements a hash table with fast lookup, where each kmer is represented
by a single bit.  Required memory is 4^k/8. 136G is needed for K=20,
550G for K=21
 */

public:

  //! ctor
  KmerSingleBitHashTable(int k);

  //! dtor
  ~KmerSingleBitHashTable();

#define SLOWER

#ifdef SLOWER
  //! insert a kmer into the hash table
  void insert(uint64_t &code)  {
    if (!find(code)) {
        ++m_setbit_count;
        setBit(code);
    }
  }

#else
  //alternate version; is this any faster?
  //! insert a kmer into the hash table
  void insert(uint64_t &code)  {
    if (!find(code)) {
        ++m_setbit_count;
        setBit(code);
    }
  }
#endif


  //! returns the number of kmers by traversing the table
  void getCount(uint64_t &is_set, uint64_t &not_set) const;

  //! returns the number of kmers by accessing a private variable
  uint64_t getSetCount() {
    return m_setbit_count;
  }

  //! returns the number of kmers that are not in the table
  //!by accessing a private variable
  uint64_t getUnsetCount() {
    return m_maxnumofkmers - m_setbit_count;
  }

  //! get possible number of kmers
  uint64_t getMaxPossibleKmerCount() {
    return m_maxnumofkmers;
  }

  //! returns the percentage of possible kmers that are
  //! actually in the DB
  double getPercentPresent() {
    return 100.0 * (double)m_setbit_count/(double)m_maxnumofkmers;
  }

  //! prints ascii representation of kmers to cout
  void printTextKmers() const;


  //! returns true if the kmer is in the hash table
  bool find(uint64_t kmer) const {
    return m_tab[DIVWORDSIZE(kmer)] & ITHBIT(MODWORDSIZE(kmer));
  }


  //======================================================================
  //public methods below this line are legacy; some will go away,
  //don't use them for now

  //void printTextKmers(const char *fn, bool print_kmers = false) const;
  void printTable(const char *fn) const;

  //void printFastaFormattedKmers(const char *fn) const;
  void printBinaryKmers(const char *fn, bool print_kmers = false) const;
  void loadTable(const char *fn);


  void flipBits();


  void resetCounters() {
    m_setbit_count = 0;
  }

  void getSetBitCounts(uint64_t &fwd, uint64_t &rc) {
    fwd = m_set_fwd_count;
    rc = m_set_rc_count;
  }

  uint64_t m_maxnumofkmers;

  void save() {
    m_hash_table_2 = m_hash_table;
  }

  void restore() {
    m_hash_table = m_hash_table_2;
  }

  uint64_t size() {
    return m_hash_table.size();
  }





//private:

  uint64_t m_setbit_count;

  uint64_t m_set_fwd_count;
  uint64_t m_set_rc_count;

  std::string m_characters;

  std::vector<uint64_t> m_hash_table;
  std::vector<uint64_t> m_hash_table_2;
  uint64_t *m_tab;
  int m_k;
  int m_kk;



  void setBit(uint64_t code) {
    m_tab[DIVWORDSIZE(code)] |= ITHBIT(MODWORDSIZE(code));
  }

  char m_kmer[33];
  char m_kmer_rc[33];

  int m_rank;
  int m_size;
};


} //namespace unwords

#endif
