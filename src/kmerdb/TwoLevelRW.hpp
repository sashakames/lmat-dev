#ifndef __SORTED_DB_HPP__
#define __SORTED_DB_HPP__



#include <iostream>
#include <set>
#include <ext/hash_map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <cstdlib>
#include <string.h>


#include <metag_const.h>
#include "metag_typedefs.hpp"

#if (WITH_PJMALLOC == 0)
#define JEMALLOC_P(x) x
#endif


#undef PAGE_SIZE 
#undef MAX_PAGE



// 18-mer

#define TT_BLOCK_COUNT_18 134217728
#define BITS_PER_2ND_18 9
#define MASK_2ND_18 0x00000000000001ff
#define LENGTH_MAX_2ND_18 513

// 20-mer

#define TT_BLOCK_COUNT 536870912
#define BITS_PER_2ND 11

#define LENGTH_MAX_2ND 2049




#include <perm.h>



namespace metag {

class MyPair
{
public:
  MyPair(unsigned int f, uint32_t s) : first(f), second(s) {}

  const bool operator < (const MyPair &mp ) const {
    return (first < mp.first);
  }

  unsigned int first;
  uint32_t second;
};




  int kmer_rec_comp(const void *a, const void *b);


 
  class SortedDbRW {

public:



    SortedDb(size_t n_kmers, size_t space_size, const int n_threads)

  {

    top_tier_block = new (JEMALLOC_P(malloc)(sizeof(uint64_t)*TT_BLOCK_COUNT)) uint64_t[TT_BLOCK_COUNT];


    //    allocator_str = "\n\ndatabase built with perm-je\n\n";
    idx_config = IDX_CONFIG;

    m_n_kmers = 0;

    bzero((void*)top_tier_block, sizeof(uint64_t) *TT_BLOCK_COUNT);
    

  }



  void set_kmer_length(char c)
  {
    m_kmer_length = c;
    idx_config = IDX_CONFIG;
  }

  char get_kmer_length() {
    return m_kmer_length;
    
  }
         
  size_t size() {
    return m_n_kmers;
  }

private:


  int idx_config;

  size_t m_n_kmers;

  uint8_t m_kmer_length;

  uint64_t *top_tier_block;


  // used for taxid lists


  size_t  m_list_offset ;
  uint16_t m_cur_page;
  uint32_t m_cur_offset;



  };

}
#endif
