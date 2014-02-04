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
    

int SIZE_CLASSES = { 4, 8, 12, 16, 24, 32 , 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 2048 };

int get_class(int key)
{
  int idx=0;

  while (key > SIZE_CLASSES[idx])
    idx++;

  return SIZE_CLASSES[idx];
  
}


      char *get_addr(uint64_t kmer_rec) {
          
          
          
      }

      
                        template <int, int >
      
      void reallocate_cell<OLD,NEW> {
                  TidRecCell<OLD> *cur_cell = reinterpret_cast<TidRecCell<OLD>*>(get_addr  );
                  TidRecCell<NEW> *new_cell = new JEMALLOC_P(malloc)(sizeof(unit16_t)*new_size_class);
                  
                  char * src;
                  char * dest;
                  
                  int pos = cur_cell->find_item(in_tid, tid_count);
                  
                if (pos > -1)
                    break;
                  
                  abspos = (-pos) -1;
                  src = cur_cell->get();
                  dest = new_cell->get();
                  
                  if (abspos > 0)
                      memcpy(dest, src, (sizeof(uint16_t) * abspos));
                  
                  delete(cur_cell);
          ￼
      // store the offset
      }

    void addTaxid(uint64_t & kmer_rec, uint16_t in_tid) {
        
          int tid_count = get_tid_count(kmer_rec);

         int old_size_class = get_class(tid_count)
         int new_size_class = get_class(tid_count+1)
      if ( old_size_class > new_size_class) {


          
          switch (old_size_class) {
              case <#constant#>4:
                    reallocate_cell(in_tid, tid_count, kmer_rec)

                  
                  break;
                  
              default:
                  break;
          }
          

          
	  //	  reallocate

      } else {
	  //  add to current structure
          
          bool rc = FALSE;
          
          switch (old_size_class) {
                  
              case 4:
          TidRecCell<4> &cur_cell = reinterpret_cast<TidRecCell<4>>(get_addr(kmer))
                rc = cur_cell.insert_item();                  ￼
                  break;
                                case 8:
          TidRecCell<8> &cur_cell = reinterpret_cast<TidRecCell<8>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                                                  case 12:
          TidRecCell<12> &cur_cell = reinterpret_cast<TidRecCell<12>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                                                  case 8:
          TidRecCell<16> &cur_cell = reinterpret_cast<TidRecCell<16>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                                                  case 24:
          TidRecCell<24> &cur_cell = reinterpret_cast<TidRecCell<24>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                                                  case 32:
          TidRecCell<32> &cur_cell = reinterpret_cast<TidRecCell<32>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                                                  case 48:
          TidRecCell<48> &cur_cell = reinterpret_cast<TidRecCell<48>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                                                  case 64:
          TidRecCell<64> &cur_cell = reinterpret_cast<TidRecCell<64>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
            case 96:
          TidRecCell<96> &cur_cell = reinterpret_cast<TidRecCell<96>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 128:
                  TidRecCell<128> &cur_cell = reinterpret_cast<TidRecCell<128>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 192:
                  TidRecCell<192> &cur_cell = reinterpret_cast<TidRecCell<192>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 256:
                  TidRecCell<256> &cur_cell = reinterpret_cast<TidRecCell<256>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 384:
                  TidRecCell<384> &cur_cell = reinterpret_cast<TidRecCell<384>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 512:
                  TidRecCell<512> &cur_cell = reinterpret_cast<TidRecCell<512>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 768:
                  TidRecCell<768> &cur_cell = reinterpret_cast<TidRecCell<768>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 1024:
                  TidRecCell<1024> &cur_cell = reinterpret_cast<TidRecCell<1024>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
              case 2048:
                  TidRecCell<2048> &cur_cell = reinterpret_cast<TidRecCell<2048>>(get_addr(kmer))
                  rc = cur_cell.insert_item();                  ￼
                  break;
                  
              default:
                  break;
          }
          

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
