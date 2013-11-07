#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <omp.h>
#include "kencode.hpp"

#include "all_headers.hpp"

#include "TaxNodeStat.hpp"


#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;

size_t perm_bytes_allocd;

static char getComp(char base) {
   char c='X';
    switch (tolower(base)) {
      case 'a' : c = 't'; break;
      case 'c' : c = 'g'; break;
      case 't' : c = 'a'; break;
      case 'g' : c = 'c'; break;
      default  : c = 'X';
    }
  return c;
}

static bool isBase(char base) {
   bool res=false;
    switch (tolower(base)) {
      case 'a' : 
      case 'c' : 
      case 't' : 
      case 'g' : res=true; break;
      default  : res=false; 
    }
  return res;
}



/*
** step through each k-mer in the read in a single direction and update the "presence/absence" vector (pa_vec) when a k-mer is found
**
** h => k-mer database
** read => query read
** mer_len => k-mer length
** pa_vec => a 'bit' vector, with one element for each position in the read. Records whether the kmer for that position  is found in
**           the database.
**
** Note that the read is search in each direction separately, however, since the k-mers are stored in canonical order
** strand differentiation is not made, thus presence absence from either strand is stored in the single direction
** 
*/

void retrieve_kmer_labels(INDEXDB<uint32_t>* table, const string& read, const size_t mer_len, bool revStrnd, ofstream &ofs)
{
  bool verbose=false;
  kencode_c ken(mer_len);
  size_t read_len = read.length();
  uint64_t kmer_id = 0x0;
  bool bad_last_kmer=true;
  unsigned next_valid_pos = 0;
  
  assert(read.size() >= mer_len);  
  for(size_t j = 0; j < read.size(); ++j) {
    if(!isBase(read[j])) {
        bad_last_kmer=true;
        next_valid_pos = (j + mer_len);
        continue;
    }
    if( bad_last_kmer && (j >= mer_len - 1) && j >= next_valid_pos ) {
        const string nkmer=read.substr((j-mer_len)+1,mer_len);
	kmer_id = ken.kencode(nkmer);
	//	kmer_id = parse_dna::mer_string_to_binary(nkmer.c_str(), mer_len);

        bad_last_kmer=false;
    } else if( !bad_last_kmer) {
      kmer_id = ken.kencode(read[j]);
	     //kmer_id =parse_dna::mer_string_to_binary(read.c_str()+(j-mer_len), mer_len);
      bad_last_kmer=false;
    }
    size_t pos = j - mer_len + 1;
    if (revStrnd) {
       pos = read_len - j - 1;
    }
    //cout<<"debug kmer "<<kmer_id<<" "<<j<<" "<<label_vec[pos].first<<" "<<read[j]<<" "<<(bad_last_kmer?"true":"false")<<endl;  
    if(!bad_last_kmer) {
      TaxNodeStat<uint32_t> *h = new TaxNodeStat<uint32_t>(*table);
       h->begin(kmer_id);
       
       const uint16_t ng = h->taxidCount();
       if ( ng > 0 ) {
	 ofs<< ng ;
	         while( h->next() ) {  
	            const uint32_t tid = h->taxid();
	            //collect the unique set of tax ids
	           
	            ofs<<" "<<tid ;
	         }
       } else {
	      ofs<<" none";
       }
	    ofs<<endl;
	    delete h;
    }
  }
}

void proc_line(int ri_len, string &line, int k_size, INDEXDB<uint32_t> *table, ofstream &ofs ) {

  assert(ri_len >= 0 && ri_len <= line.length());


  for(unsigned chk = 0; (chk+ri_len) <= line.length(); chk += (ri_len+1)) {
    string rval;
    if( ri_len < line.length() ) {
      rval = line.substr(chk,ri_len);
    } else {
      rval = line;
    }

    const unsigned read_len = rval.length();
    if( static_cast<signed>(read_len) < ri_len ) {
      break;
    }
    if( read_len < k_size ) {
      if( chk == 0 ) /// read too short 
	ofs << "none!!" << endl;
      break;
    }
    retrieve_kmer_labels(table, rval, k_size, false, ofs);

       string rev_cmp(rval.size(),'\0');
       assert(rval.size()>0);
       for (int j=(rval.size()-1); j >= 0; --j) {
         const int bidx = (rval.size()- 1) - j;
         rev_cmp[j] = getComp(rval[bidx]);
       }
       retrieve_kmer_labels(table, rev_cmp, k_size, true, ofs);

  }
   
}


void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <jellylist kmerdb file> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  
}

size_t *split_file(int n_threads, ifstream &file)
{
  size_t * arr = new size_t[n_threads+1];

  arr[0] = 0;

  file.seekg(0, ios::end);
  size_t end = file.tellg();

  for (size_t i = 1; i<n_threads; i++)  {

    file.seekg( i * (end / n_threads ));
    
    string junk;
    getline(file, junk);
    
    arr[i] = file.tellg();
  }

  arr[n_threads] = end;

  return arr;

}

int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1, ri_len = -1;
   int n_threads = 1;


   float threshold = 0.0;
   string kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options;


   unsigned long read_count = 0,  max_reads = 0;
   


   while ((c = getopt(argc, argv, "o:t:k:i:d:l:q:")) != -1) {

      switch(c) {
      case 'l':
         ri_len = atoi(optarg);
         break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 'd':
         kmer_db_fn = optarg;
         break;
      case 'q':

	max_reads = atoi(optarg);
	break; 
      case 't':
        n_threads = atoi(optarg);
        break;
      case 'o':
	ofbase = optarg;
	break;	

      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }

   if (kmer_db_fn == "" || query_fn == "")  {
     usage(argv[0]);
     return -1;

   }


   cout << "Start kmer DB load\n";
   //KmerDB *table;
   INDEXDB<uint32_t> *taxtable;
   TaxNodeStat<uint32_t>* table = NULL;
   


#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDB<uint32_t>*, std::size_t> ret = mfile.find<INDEXDB<uint32_t > >("KmerDB");

     taxtable = ret.first;
     k_size = taxtable->get_kmer_length();
     cout << "k size:  " << k_size  <<   endl ;
#else

#if WITH_PJMALLOC == 1


     perm(&taxtable, sizeof(taxtable));
     mopen(kmer_db_fn.c_str(), "r", MMAP_SIZE);

     k_size = taxtable->get_kmer_length();


#endif


#endif


   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;



   assert(k_size > 0 );
   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   size_t *arr = split_file(n_threads, ifs);
   ifs.close();

   bool finished;

   size_t pos = 0 ;


     //unsigned cnt = 0;
   finished = false;

   ofstream ofs;
  
   StopWatch clock;
   clock.start();
   
   omp_set_num_threads(n_threads);
	
#pragma omp parallel shared(arr, k_size, query_fn, ofbase, taxtable,  max_reads)  private(ifs, ri_len, finished, pos, ofs, ofname, line, read_count)

   {

     read_count = 0;
     ri_len = -1;
     


     finished = false;
     int use_len = -1; 

   StopWatch clock;
   clock.start();
  
   size_t read_count = 0; 

     //     cout<<"Read query file: "<<query_fn<<endl;
     ifs.open(query_fn.c_str());
     if(!ifs) {
         cerr<<"did not open for reading: "<<query_fn<<endl;
         exit(-1);
     }
     
     ofname = ofbase;

     std::stringstream outs;

     outs << omp_get_thread_num();
     ofname += outs.str();

     ofname += ".out" ;

     ofs.open(ofname.c_str());

     ifs.seekg(arr[omp_get_thread_num()]);

     string read_buff;

     while (!finished)   {


       getline(ifs, line);

       pos = ifs.tellg();

       if ((pos >= arr[1+omp_get_thread_num()]) || (pos == -1)) {
	      finished = true;
       } 

       //#pragma omp critical    
       //cout << omp_get_thread_num()  << " - pos:" << ifs.tellg() << "\n"; 

       if (line[0] != '>' && line.length() > 1) {
          read_buff += line;
       } 
       if( (line[0] == '>' || finished) && read_buff.length() > 0 ) {

          if( ri_len == -1 ) {
            // if not specified use first read as the fixed length
            use_len = read_buff.length();
          } else {
            use_len = ri_len;
          } 
	       //#pragma omp critical
	       ofs<<"read: "<<read_buff<< " pos: " << ifs.tellg() <<  endl;
       //split read into multiples of "ri_len" (read interval length)

	       proc_line(use_len, read_buff, k_size, taxtable, ofs);
          read_buff="";
	 read_count ++;
	 if (read_count == max_reads)
	   finished = true;

       }

     }
     ofs.close();
   }

   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
