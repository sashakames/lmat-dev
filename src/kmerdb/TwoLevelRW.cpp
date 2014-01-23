#include "SortedDb.hpp"

#include <cassert>
#include <fstream>
#include <sstream>
#include <set>
#include <ext/hash_map>
#include <queue>

#include <metag_const.h>
#include "KmerFileMetaData.hpp"


using namespace std;
using namespace metag;

int metag::kmer_rec_comp(const void *a, const void *b)
{
  const kmer_record *kr_a = reinterpret_cast<const kmer_record*>(a);
  const kmer_record *kr_b = reinterpret_cast<const kmer_record*>(b);

  return kr_a->kmer_lsb - kr_b->kmer_lsb;

}

size_t ext_taxids = 0 ;
size_t singletons = 0;

size_t reduced_kmers = 0;
size_t cut_kmers = 0;

int size_classes = { 4, 8, 12, 16, 24, 32 , 48, 64, 96, 128, 192, 256, 385, 512, 786, 1024, 2048 };


void SortedDbRW::addKmer(uint64t kmer)
{

  // look up k-mer

    size_t top_index =  (kmer >> BITS_PER_2ND);  // & 0x0000000007ffffff;

    

    // check the slot if it has been written to yet

    if (top_tier_block[top_index] == 0) {
	

      // if empty allocate a new slot
      top_tier_block[top_index] = (uint16_t*)malloc(size_classes * sizeof(uint16_t));
      top_tier_block[top_index][0] = 1;
      

    } else


    
    kmer_table[m_list_offset].kmer_lsb = kmer_lsb_in;
    
    assert (top_tier_block[top_index] >> 48 == start_count);

    uint16_t tmp_tid_count = tid_count;

    set<uint32_t> write_set;
    priority_queue<MyPair> taxid_q;
     
    if (tid_cutoff > 0 && tid_count > tid_cutoff) {

      if (species_map.size() == 0) {
	tmp_tid_count = 0;

	for (uint16_t k=0; k<tid_count; k++) {
	  // scan through and dump
	  assert(fread(&tid, 4, 1, in) == 1);        
	}
      } else {

	if (strainspecies) {

	  cout << "functionality disabled!\n";
	  exit(1);
	  for (uint16_t k=0; k<tid_count; k++) {

	    
	    assert(fread(&tid, 4, 1, in) == 1);        
	    
	    if (species_map.find(tid) == species_map.end()) {
	      // alternate version - try to not write the tid because doing so might be keeping the accuracy down 
	      write_set.insert(tid);
	      
	    } else {

	    uint32_t mapped_tid = species_map[tid];
	    write_set.insert(mapped_tid);
	    
	    }
	  }
	
	  tmp_tid_count = write_set.size();
      
	// if not enough reduction , then we revert to the old method
	  if (tmp_tid_count > tid_cutoff)
	    tmp_tid_count = 0;
	} else {
	  	  
	  for (uint16_t k=0; k<tid_count; k++) {

	    assert(fread(&tid, 4, 1, in) == 1);        

	    const MyPair  pp(species_map[tid], tid);
	    taxid_q.push(pp);
	    
	    
	  }

	  // iterate on taxid priority queue in batches of taxons with          
	  // the same rank                                                  
	
	  while(!taxid_q.empty()) {

	    // get current priorty                                            
	    int cur_priority = taxid_q.top().first;
	      
	    // pull out elements that match top priority                        
	    while(taxid_q.top().first == cur_priority) {
	      
	      const MyPair res = taxid_q.top();
	      taxid_q.pop();
		
	      //	      write_set.insert(res.second);

	      if (taxid_q.empty())
		break;
	    }
	    // if we are under the cut, then we can stop                        
	    if (taxid_q.size() <= tid_cutoff) {
	      //            cout << "Cut to rank: " << cur_priority << " org \
	      // cout  " << m_taxid_count << " new count" <<  m_filtered_list.size()  << "\n";  
	      tmp_tid_count = write_set.size();

	      break;
	      
	    }
	    //else {
                // prepare for next batch of element poping                     
	    //write_set.clear();
	    //}
	    
	  }
	  if (taxid_q.size() == 0) {
	    tmp_tid_count = 1;
	    const MyPair pp(1, 1);
	    taxid_q.push(pp);
	    cut_kmers++;
	  }
	  
	}
      }
      
    }

    if (tmp_tid_count > 1) {
      
      if (16+m_cur_offset+tmp_tid_count*4 > PAGE_SIZE) {
	  
	m_cur_page ++;
	m_cur_offset = 0;
	
      }
	
      kmer_table[m_list_offset].page_id = m_cur_page;
      kmer_table[m_list_offset].page_offset = m_cur_offset;
    }
    else if (tid_count == 1) {
      
      singletons++;
	
      assert(fread(&tid, 4, 1, in) == 1);        
	//	assert (tid <= MAX_TID && tid != INVALID_TID_2 );
      

      kmer_table[m_list_offset].page_id = MAX_PAGE;

      if (p_br_map) {
	uint16_t tid_16 = (*p_br_map)[tid];
	if (!(tid_16 < p_br_map->size()+3))  // pad for starting the count at 2
	{
	  cout << "bad single: " << tid << " " << tid_16 << "\n";
	  assert(0);
	}
	kmer_table[m_list_offset].page_offset = tid_16;
      }
      else {

	kmer_table[m_list_offset].page_offset = tid;
      }

    } else if (tmp_tid_count == 1) {

           
      tid = taxid_q.top().second;
      kmer_table[m_list_offset].page_id = MAX_PAGE;

      if (p_br_map) {
	uint16_t tid_16 = (*p_br_map)[tid];
	if (!(tid_16 < p_br_map->size()+1))
	  {
	    cout << "bad single: " << tid << " " << tid_16 << "\n";
	    assert(0);
	  }
	kmer_table[m_list_offset].page_offset = tid_16;
      } else {
	kmer_table[m_list_offset].page_offset = tid;

      }
      reduced_kmers++;
    } else {
      assert(tmp_tid_count == 0);      
      kmer_table[m_list_offset].page_id = MAX_PAGE;
      kmer_table[m_list_offset].page_offset = 1;
      cut_kmers++;
    }
    
    
    m_list_offset++;
    start_count++;

    assert(start_count <= LENGTH_MAX_2ND );

	//write the kmer
      //      memcpy(m_data[m_cur_page]+m_cur_offset, &kmer, 8);
    if (taxid_q.size() > 0 && tmp_tid_count > 0)	{

	// check for no reduction

      reduced_kmers++;
      if (kmer % 4096 == 0) {
	mcpyinsdb(kmer, 8);
	m_cur_offset += 8;
	  
      }	
	
      mcpyinsdb(tmp_tid_count, 2);
      m_cur_offset += 2;

      for (int i=0; i < taxid_q.size(); i++) {
	
	tid = taxid_q.top().second;
	taxid_q.pop();
	
	if (p_br_map) {

	  uint16_t tid_16 = (*p_br_map)[tid];
	  if (!(tid_16 < p_br_map->size()+1)) {
	    cout << "bad set: " << tid << " " << tid_16 << "\n";
	    assert(0);
	  }
	  assert (tid_16 > 0);
	      mcpyinsdb(tid_16,2);
	      m_cur_offset += 2;
	}  else {

	      mcpyinsdb(tid,4);
	      m_cur_offset += 4;
	    }

	ext_taxids++;
      } 
      
      
    }
    // no attempt to reduce list; copy in the taxid list
    else if (write_set.size() == 0 && tmp_tid_count > 1) {
	  
      if (kmer % 4096 == 0) {
	mcpyinsdb(kmer, 8);
	m_cur_offset += 8;
      }
      
      mcpyinsdb(tid_count, 2);
      m_cur_offset += 2;

      uint16_t tmpcount;

      if (tmpcount != count_marker) {

	cout << "changed to " << tmpcount << " at " << kmer << "\n";
	count_marker = tmpcount;

      }
      
      //write the tuples
      for (uint16_t k=0; k<tid_count; k++) {
	
	  ext_taxids++;
	  assert(fread(&tid, 4, 1, in) == 1);        

	if (p_br_map) {


	  uint16_t tid_16 = (*p_br_map)[tid];
	  if (!(tid_16 < p_br_map->size()+1)) {
	    cout << "bad read: " << tid << " " << tid_16 << "\n";
	    assert(0);
	  }
	  
	  assert (tid_16 > 0);

	  mcpyinsdb(tid_16, 2);
	  m_cur_offset += 2;
	} else {

	  mcpyinsdb(tid, 4);
	  m_cur_offset += 4;
	}

      }
    } else {
      assert (tid_count == 1 || tmp_tid_count < 2);
      
    }

    
    if (use_tax_histo_format) {
	
	if ((i+1) % TAX_HISTO_SANITY_COUNT == 0) {
	  assert(fread(&test, sizeof(uint64_t), 1, in) == 1);
	  assert(test == sanity);
	}  
      } else {
	if ((i+1) % KMER_SANITY_COUNT == 0) {
	  assert(fread(&test, sizeof(uint64_t), 1, in) == 1);
	  assert(test == sanity);
	}
      }
    last_kmer = kmer;
  }

  m_n_kmers ++;
 

}


