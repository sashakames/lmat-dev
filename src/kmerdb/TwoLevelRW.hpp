#ifndef __TWOLEVEL_DB_HPP__
#define __TWOLEVEL_DB_HPP__



#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <cstdlib>
#include <string.h>
#include <stdio.h>

#include <perm.h>

#include <metag_const.h>
#include "metag_typedefs.hpp"

#include "ListCellBase.hpp"

#include "kencode.hpp"



using namespace kencode_ns;

#if (WITH_PJMALLOC == 0)
#define JEMALLOC_P(x) x
#endif


#undef PAGE_SIZE 
#undef MAX_PAGE


typedef    std::map<uint64_t, uint64_t, std::less<uint64_t>, PERM_NS::allocator<std::pair<const uint64_t, uint64_t> > > kmer_map;

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


  int  SIZE_CLASSES[] = { 4, 8, 12, 16, 24, 32 , 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 2048 };

 
  class TwoLevelRW {

  public:
    
    TwoLevelRW(int in_kmer_len) {
          
      m_kmer_length = in_kmer_len;
          
      }
      
      

int get_class(int key)
{
  int idx=0;

  while (key > SIZE_CLASSES[idx])
    idx++;

  return SIZE_CLASSES[idx];
  
}

      char *mmap_base_addr;


    int get_tid_count(uint64_t kmer_rec)
    {

      return ( 0x0003FF0000000000 & kmer_rec )>> 40;
    }
      
      void set_base_addr(char *in_addr)
      {
          
          
          mmap_base_addr = in_addr;
      }
      

      char *get_addr(uint64_t kmer_rec) {
          
          size_t offset = kmer_rec & 0x00000FFFFFFFFFF ;
          
          return mmap_base_addr + offset;
      }

      
      void set_offset(char *addr, uint64_t &kmer_rec) {
          
          size_t offset = addr - mmap_base_addr;
          
          kmer_rec = kmer_rec | (offset & 0x000000FFFFFFFFFF) ;
          

      }

    template <int SZCLASS> 

    int  grow_cell(uint64_t & kmer_rec, uint16_t tid, int count)
    {
      
      ListCellBase<uint16_t, SZCLASS> *cur_cell = reinterpret_cast<ListCellBase<uint16_t, SZCLASS>*>(get_addr(kmer_rec));
      kmer_rec = (kmer_rec & (0x000000FFFFFFFFFF) ) | ( ((uint64_t)(count + 1 )) << 40);
      return  cur_cell->insert_item(tid, count );
  }

      
    template <int OLD, int NEW>
      void reallocate_cell(uint16_t in_tid, int tid_count  , uint64_t & kmer_rec ) {
      ListCellBase<uint16_t, OLD> *cur_cell = reinterpret_cast<ListCellBase<uint16_t, OLD>*>(get_addr(kmer_rec));
      ListCellBase<uint16_t, NEW> *new_cell = new(JEMALLOC_P(malloc)(sizeof(uint16_t)*NEW)) ListCellBase<uint16_t, NEW>  ;
      
      char * src;
      char * dest;
      
      int pos = cur_cell->find_item(in_tid, tid_count);
                  
      if (pos > -1)
	return;
      
      int abspos = (-pos) -1;
      src = cur_cell->get();
      dest = new_cell->get();
                  
      if (abspos > 0)
	memcpy(dest, src, (sizeof(uint16_t) * abspos));
      
      new_cell->m_items[abspos] = in_tid;
                
      if ( tid_count > abspos)
	delete(cur_cell);
      set_offset(new_cell->get(), kmer_rec);
    }

  void addTaxid(uint64_t & kmer_rec, uint16_t in_tid) {
    
    int tid_count = get_tid_count(kmer_rec);
    
    int old_size_class = get_class(tid_count);
    int new_size_class = get_class(tid_count+1);

    bool rc = false;
    
      if ( old_size_class > new_size_class) {
	
	switch (old_size_class) {
	case 4:
	  reallocate_cell<4,8>(in_tid, tid_count, kmer_rec);
	  
	  break;
	case 8:
	  reallocate_cell<8,12>(in_tid, tid_count, kmer_rec);
	  break;
	case 12:
	  reallocate_cell<12, 16>(in_tid, tid_count, kmer_rec);
	case 16:
	  reallocate_cell<16, 24>(in_tid, tid_count, kmer_rec);
	case 24:
	  reallocate_cell<24, 32>(in_tid, tid_count, kmer_rec);
	case 32:
	  reallocate_cell<32, 48>(in_tid, tid_count, kmer_rec);
	case 48:
	  reallocate_cell<48, 64>(in_tid, tid_count, kmer_rec);
	case 64:
	  reallocate_cell<64, 96>(in_tid, tid_count, kmer_rec);
	case 96:
	  reallocate_cell<96, 128>(in_tid, tid_count, kmer_rec);
	case 128:
	  reallocate_cell<128, 192>(in_tid, tid_count, kmer_rec);
	case 192:
	  reallocate_cell<192, 256>(in_tid, tid_count, kmer_rec);
	case 256:
	  reallocate_cell<256, 384>(in_tid, tid_count, kmer_rec);
	case 384:
	  reallocate_cell<384, 512>(in_tid, tid_count, kmer_rec);
	case 512:
	  reallocate_cell<512, 768>(in_tid, tid_count, kmer_rec);
	case 768:
	  reallocate_cell<768, 1024>(in_tid, tid_count, kmer_rec);
	case 1024:
	  reallocate_cell<1024, 2048>(in_tid, tid_count, kmer_rec);
	  break;
          
	default:
	  assert(0);
	  break;
	}
	

          
	  //	  reallocate

      } else {
	  //  add to current structure
          

          
          switch (old_size_class) {
                  
	  case 4:
	    rc = grow_cell<4>(kmer_rec, in_tid, tid_count );
	    break;
	  case 8:
	    rc = grow_cell<8>(kmer_rec, in_tid, tid_count );
	    break;
	  case 12:
	    rc = grow_cell<12>(kmer_rec, in_tid, tid_count );
	    break;
	  case 16:
	    rc = grow_cell<16>(kmer_rec, in_tid, tid_count );
	    break;
	  case 24:
	    rc = grow_cell<24>(kmer_rec, in_tid, tid_count );
	    break;
	  case 32:
	    rc = grow_cell<32>(kmer_rec, in_tid, tid_count );
	    break;
	  case 48:
	    rc = grow_cell<48>(kmer_rec, in_tid, tid_count );
	    break;
	  case 64:
	    rc = grow_cell<64>(kmer_rec, in_tid, tid_count );
	    break;
	  case 96:
	    rc = grow_cell<96>(kmer_rec, in_tid, tid_count );
	    break;
	  case 128:
	    rc = grow_cell<128>(kmer_rec, in_tid, tid_count );
	    break;
	  case 192:
	    rc = grow_cell<192>(kmer_rec, in_tid, tid_count );
	    break;
	  case 256:
	    rc = grow_cell<256>(kmer_rec, in_tid, tid_count );
	    break;
	  case 384:
	    rc = grow_cell<384>(kmer_rec, in_tid, tid_count );
	    break;
	  case 512:
	    rc = grow_cell<512>(kmer_rec, in_tid, tid_count );
	    break;
	  case 768:
	    rc = grow_cell<768>(kmer_rec, in_tid, tid_count );
	    break;
	  case 1024:
	    rc = grow_cell<1024>(kmer_rec, in_tid, tid_count );
	    break;
	  case 2048:
	    rc = grow_cell<2048>(kmer_rec, in_tid, tid_count );
	    break;
              default:
                  break;
                  

	  }
          
      }

      assert(rc);

  }





  void addKmer(uint64_t kmer, uint16_t taxid) {
      
      //  size_t top_index =  (kmer >> BITS_PER_2ND);  // & 0x0000000007ffffff;
  
    // Use whatever index we got for now....
        
        addTaxid(m_kmer_map[kmer], taxid);
        

    }

      

uint64_t read_encode(FILE *f, kencode_c &ken) {

  char buf[33];


  fscanf(f,"%s", buf);

  if (strlen(buf) > 0) {
    uint64_t kmer = ken.kencode(buf);
    //    cout << kmer << "\n";
  
    return kmer;
  }
  else
    return ~0;

}

      
      void addData(const char *filename, uint16_t taxid) {
	
	kencode_c ken(m_kmer_length);
	
	FILE *f = fopen(filename, "r");
	
	  uint64_t kmer = 0;
	  
	  while (~kmer) {
	    
	    kmer = read_encode(f, ken); 
	    addKmer(kmer, taxid);

	  }
	  
	  fclose(f);
          
      }


  private:


  size_t m_n_kmers;

  uint8_t m_kmer_length;

//  uint64_t *top_tier_block;
      

    kmer_map m_kmer_map;

  // used for taxid lists



  };

}
#endif
