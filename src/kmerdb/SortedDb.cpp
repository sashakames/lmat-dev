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

template <class tid_T>
void SortedDb<tid_T>::add_data(const char *filename, size_t stopper = 0, bool use_tax_histo_format = true, bitreduce_map_t *p_br_map = NULL   ,  my_map &species_map = NULL, int tid_cutoff = 0, bool strainspecies = false)
{
  

  static uint64_t last_kmer = 0;

  FILE *in = fopen(filename, "r");
  assert(in != NULL);
  //assert(fread(&kmer_count, 8, 1, fp) == 1);
  //for (uint64_t j=0; j<kmer_count; j++) 
  fseek(in, 0, SEEK_END);
  long f = ftell(in);
  fseek(in, 0, SEEK_SET);

   KmerFileMetaData metadata;
   metadata.read(in);
   if (use_tax_histo_format)
     assert(metadata.version() == TAX_HISTO_VERSION);

   uint64_t kmer_ct = metadata.size();
   
   uint64_t sanity = ~0, test;

   if (stopper == 0)
     stopper = ~0;

  uint64_t kmer;
  uint32_t tid;
  uint16_t tid_count; //, genome_count, p_count, tuple_count;

  uint32_t tid_count_32;

  ////uint16_t genome_count;

  static long long int start_count;
  static long long int start_offset;

  static uint16_t count_marker = 0;

  //  cout << "stopper set to: " << stopper << "\n";

  for  ( uint64_t i=0; i<kmer_ct; i++)  {

    if (i > stopper)
      break;

    //loop exit condition
    if (ftell(in) == f) break;

    //read kmer, taxid count, and tuple count
    assert(fread(&kmer, 8, 1, in) == 1);        

    if ((last_kmer > 0) &&  (kmer <= last_kmer)) {
      cout << "Kmers arriving out of order.  New: " << kmer << " last: " << last_kmer << "\n";
      exit(1);
    }


    if (use_tax_histo_format) {
      assert(fread(&tid_count, 2, 1, in) == 1);
    } else {
      assert(fread(&tid_count_32, 4, 1, in) == 1);

      tid_count = (uint16_t)tid_count_32;

    }

    // First get the mapping - the msb for the kmer - assuming 27 for now
    // TODO: make this configurable at runtime 

    size_t top_index =  (kmer >> BITS_PER_2ND);  // & 0x0000000007ffffff;


    // check the slot if it has been written to yet

    if (top_tier_block[top_index] == 0) {
	


      // if empty write the current offset - this will be the start of the short list within the second tier
      start_offset = m_list_offset; 
      start_count = 1;  
    }

    uint16_t kmer_lsb_in = MASK_2ND & kmer;

    top_tier_block[top_index] = ((uint64_t) start_count << 48) | start_offset;
    
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
    if (taxid_q.size() > 1 && tmp_tid_count > 1)	{
      
	// check for no reduction

      reduced_kmers++;
      if (kmer % 4096 == 0) {
	mcpyinsdb(kmer, 8);
	m_cur_offset += 8;
	  
      }	
	
      mcpyinsdb(tmp_tid_count, 2);
      m_cur_offset += 2;
      

      for (int i=0; i < tmp_tid_count; i++) {
	
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

  m_n_kmers += kmer_ct;
  
  fclose(in);

  cout << "storage used for tax ids, counts, etc.."  << m_cur_page << " - "  << m_cur_offset << "\n";
  cout << "kmer count: " << m_n_kmers << "\n";
  cout << "offset at: " << m_list_offset << "\n";
  cout << "singletons: " << singletons << "\n";
  cout << "taxids in storage: " << ext_taxids << "\n";
  cout << "kmers reduced: " << reduced_kmers << "\n";
  cout << "kmers cut to 1: " << cut_kmers << "\n";

}


template <class tid_T>
void SortedDb<tid_T>::add_data(const char *filename, size_t stopper, const int n_threads, const int thread_no)
{
  
  assert(0);
  cout << "adding file " << filename << " thread number " << thread_no << endl;

  FILE *in = fopen(filename, "r");
  assert(in != NULL);
  //assert(fread(&kmer_count, 8, 1, fp) == 1);
  //for (uint64_t j=0; j<kmer_count; j++) 
  fseek(in, 0, SEEK_END);
  long f = ftell(in);
  fseek(in, 0, SEEK_SET);


  uint64_t kmer;
  uint32_t tid;
  uint16_t tid_count, genome_count, p_count, tuple_count;

  ////uint16_t genome_count;

  size_t i =0;
  


  long long int start_count = -1;
  long long int start_offset;


  while (true) {

    i++;
    
    if (i == stopper)
      break;

    //loop exit condition
    if (ftell(in) == f) break;

    //read kmer, taxid count, and tuple count
    assert(fread(&kmer, 8, 1, in) == 1);        
    assert(fread(&genome_count, 2, 1, in) == 1);        
    assert(fread(&tid_count, 2, 1, in) == 1);        
    assert(fread(&tuple_count, 2, 1, in) == 1);        

    

    /*
    if(genome_count == 1 && tid_count == 1 && tuple_count == 1) {
        
      assert(fread(&tid, 4, 1, in) == 1);        
      assert (tid <= MAX_TID && tid != INVALID_TID_2 );
        
      //        memcpy(m_data[m_cur_page]+m_cur_offset, &tid, 4);


      uint16_t present;
      assert(fread(&present, 2, 1, in) == 1);        
      assert(present == 1); 
    */
        
      //      (*this)[kmer] =  pair<uint32_t, uint8_t>(tid , MAX_PAGE);

	//set entry in hash: kmer -> <offset, page>

      // 	(*this)[kmer] = pair<uint32_t, uint8_t>(m_cur_offset, m_cur_page);

      // First get the mapping - the msb for the kmer - assuming 27 for now
      // TODO: make this configurable at runtime 



    uint32_t top_index =  (kmer >> BITS_PER_2ND); // & 0x0000000007ffffff;

      // check the slot if it has been written to yet

      if (top_tier_block[top_index] == 0) {
	

	// if empty write the current offset - this will be the start of the short list within the second tier
	start_offset = list_offset_arr[thread_no]; 
	start_count = 1;  
      } else if (start_count == -1) {
	start_count = top_tier_block[top_index] >> 48;
	assert (start_count < LENGTH_MAX_2ND && start_count > 0);
      }
	


      uint16_t kmer_lsb_in = 0x00000000000001ff & kmer;
	
      top_tier_block[top_index] = ((uint64_t) start_count << 48) | start_offset;

      kmer_table[list_offset_arr[thread_no]].kmer_lsb = kmer_lsb_in;
      
      assert (top_tier_block[top_index] >> 48 == start_count);



      
      if (tuple_count > 1) {
	if (16+m_cur_offset_arr[thread_no]+tuple_count*6 > PAGE_SIZE) {
	
	  m_cur_page_arr[thread_no] += n_threads;
	  m_cur_offset_arr[thread_no] = 0;
	  cout << "Page Added: " << m_cur_page_arr[thread_no] << "-  ";
	  cout << "Kmers ingested: " << i << "\n"; 

	}

	kmer_table[list_offset_arr[thread_no]].page_id = m_cur_page_arr[thread_no];
	kmer_table[list_offset_arr[thread_no]].page_offset = m_cur_offset_arr[thread_no];
      }
      else {

	assert(fread(&tid, 4, 1, in) == 1);        
	//	assert (tid <= MAX_TID && tid != INVALID_TID_2 );

	uint16_t present;
	assert(fread(&present, 2, 1, in) == 1);        

	kmer_table[list_offset_arr[thread_no]].page_id = MAX_PAGE;
	kmer_table[list_offset_arr[thread_no]].page_offset = tid;
      }


      list_offset_arr[thread_no]++;
      start_count++;

      assert(start_count <= LENGTH_MAX_2ND);


	//write the kmer
      //      memcpy(m_data[m_cur_page]+m_cur_offset, &kmer, 8);

      if (tuple_count > 1) {

	if (kmer % 4096 == 0) {
	  mcpyinsdbt(kmer, 8, thread_no);
	  m_cur_offset_arr[thread_no] += 8;

	  
	}
	
      //write taxid count, genome count, and tuple count
      // memcpy(m_data[m_cur_page]+m_cur_offset, &tid_count, 2);
	//	mcpyinsdb(tid_count, 2);
	//	m_cur_offset += 2;
      //      memcpy(m_data[m_cur_page]+m_cur_offset, &genome_count, 2);
	//	mcpyinsdb(genome_count, 2);
	//	m_cur_offset += 2;
      //      memcpy(m_data[m_cur_page]+m_cur_offset, &tuple_count, 2);

	mcpyinsdbt(tuple_count, 2, thread_no);
	m_cur_offset_arr[thread_no] += 2;

      //write the tuples
	for (uint16_t k=0; k<tuple_count; k++) {

	  assert(fread(&tid, 4, 1, in) == 1);        
	  assert (tid <= MAX_TID && tid != INVALID_TID_2 );
	
	//        memcpy(m_data[m_cur_page]+m_cur_offset, &tid, 4);
	  mcpyinsdbt(tid, 4, thread_no);
	  m_cur_offset_arr[thread_no] += 4;

	  uint16_t present;
	  assert(fread(&present, 2, 1, in) == 1);        

	//        memcpy(m_data[m_cur_page]+m_cur_offset, &present, 2);
	  //	  mcpyinsdb(present, 2);

	  //	  m_cur_offset += 2;
	  //cout << tid << " ";
	  /*
	    assert(fread(&p_count, 2, 1, in) == 1);        
	    memcpy(m_data[m_cur_page]+m_cur_offset, &p_count, 2);
	    m_cur_offset += 2;
	  */
	}
	//cout << endl;
      }
  }

#pragma omp critical
  {
    m_n_kmers += (i-1);

  }

  cout << "kmer count: " << m_n_kmers << "\n";
  
  fclose(in);


}

template class SortedDb<DBTID_T>;

