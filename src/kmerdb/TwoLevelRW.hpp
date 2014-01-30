#ifndef __TWOLEVEL_DB_HPP__
#define __TWOLEVEL_DB_HPP__



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


 
  class TwoLevelRW {

public:
    

int size_classes = { 4, 8, 12, 16, 24, 32 , 48, 64, 96, 128, 192, 256, 384, 512, 786, 1024, 2048 };

int get_class(int key)
{
  int idx=0;

  while (key > size_classes[idx])
    idx++;

  return size_classes[idx];
  
}




    void addTaxid(uint64_t & kmer_rec)
    {
      int tid_count = get_tid_count(kmer_rec);

      if (get_class(tid_count) > get_class(tid_count+1) {
	  
	  //	  reallocate

	} else {
	  //  add to current structure

	  
	}

    }





    void addKmer(uint64_t kmer, uint32_t taxid) {
      
        size_t top_index =  (kmer >> BITS_PER_2ND);  // & 0x0000000007ffffff;
  

    }

    void addData(char *filename, uint32_t taxid) ;


  private:


  size_t m_n_kmers;

  uint8_t m_kmer_length;

  uint64_t *top_tier_block;


  // used for taxid lists



  };

}
#endif
