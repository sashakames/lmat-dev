#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <omp.h>
#include "kencode.hpp"
#include "all_headers.hpp"






using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;

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

void print_mer_counts(const KmerDB* h, const string& read, const size_t mer_len, vector<bool>& pa_vec, bool revStrnd)
{
//dah debug:
  kencode_c ken(mer_len);
  size_t read_len = read.length();
  uint64_t kmer_id; 
  bool bad_last_kmer=true;
  unsigned next_valid_pos = 0;
  for(size_t j = 0; j < (read.size()-mer_len+1); ++j) {
    if(!isBase(read[j])) {
        bad_last_kmer=true;
        next_valid_pos = (j + mer_len);
        continue; 
    }
    if( bad_last_kmer && (j >= mer_len) && j >= next_valid_pos ) { 
        const string nkmer=read.substr((j-mer_len)+1,mer_len);
        kmer_id = ken.kencode(nkmer);
        bad_last_kmer=false;
    } else if( !bad_last_kmer) {
        kmer_id = ken.kencode(read[j]);
        bad_last_kmer=false;
    }
    if(!bad_last_kmer) {
       KmerDB* tfind = const_cast<KmerDB*>(h);
       if (tfind->lookup(kmer_id)) {
         size_t pos = j;
         if (revStrnd) {
            pos = read_len - j;
         }
         pa_vec[pos] = true;
       }
    }
  }
}

void proc_line(int ri_len, string &line, int k_size, KmerDB *table, ofstream &ofs) {

    for(unsigned chk = 0; (chk+ri_len) <= line.length(); chk += (ri_len+1)) {
       string rval = line.substr(chk,ri_len);
       const unsigned read_len = rval.length();
       if( static_cast<signed>(read_len) < ri_len ) {
         break;
       }  
       vector<bool> pa_vec(read_len);

       print_mer_counts(table, rval, k_size, pa_vec, false);

       string rev_cmp(rval.size(),'\0');
       assert(rval.size()>0);
       for (int j=(rval.size()-1); j >= 0; --j) {
         const int bidx = (rval.size()- 1) - j;
         rev_cmp[j] = getComp(rval[bidx]);
       }
       print_mer_counts(table, rev_cmp, k_size, pa_vec, true);
       unsigned kmer_cnt = 0;
       for(unsigned i = 0; i < read_len; ++i) {
          if(pa_vec[i]) {
            ++kmer_cnt;
          }
       }
       const float fval = static_cast<float>(kmer_cnt) / static_cast<float>(read_len);
       ofs<<fval<<" "<<kmer_cnt<<" "<<read_len<< endl;
    }
}


int *split_file(int n_threads, ifstream &file)
{
  int * arr = new int[n_threads+1];

  arr[0] = 0;

  file.seekg(0, ios::end);
  unsigned long int end = file.tellg();

  for (int i = 1; i<n_threads; i++)  {

    file.seekg( i * (end / n_threads ));
    
    string junk;
    file>>junk;
    
    arr[i] = file.tellg();
  }

  arr[n_threads] = end;

  return arr;

}

void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <jellylist kmerdb file> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1, ri_len = -1;
   int n_threads = 0;
#if WITH_PJMALLOC == 1
   bool restore = false;
#endif 

   double refresh = 0;

   string kmer_db_fn, query_fn, ofname, ofbase;


   size_t mmap_size = 0;


   while ((c = getopt(argc, argv, "k:i:d:l:t:r o:z:s:")) != -1) {

      switch(c) {
#if WITH_PJMALLOC == 1
      case 'r':
        restore = true;
        break;
#endif 
      case 's':
	mmap_size = strtoull(optarg, NULL, 10);
	cout << "mmap size: " << mmap_size << endl;
	break;

      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'z':
	refresh = atof(optarg);
	break;
      case 'l':
	ri_len = atoi(optarg);
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 'd':
        kmer_db_fn = optarg;
        break;
      case 'o':
        ofbase = optarg;
        break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }

   if (ofbase == "" || n_threads == 0 || kmer_db_fn == "" || query_fn == "")  {
     usage(argv[0]);
     return -1;

   }

   cout << "Start kmer DB load\n";
   KmerDB *table;

#if WITH_PJMALLOC == 1
   if (restore) {
     perm(&table, sizeof(table));
     mopen(kmer_db_fn.c_str(), "r", mmap_size);
   } else
#endif
   {
     table = new KmerDB;
     //bool use_offsets = true;
     table->ingest(kmer_db_fn.c_str() /*, use_offsets*/ );
   }

   cout << "End kmer DB load\n";
   cout << "DB size is " << table->size() << endl;


   k_size = table->get_kmer_length();

   assert(k_size > 1);

   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   int *arr = split_file(n_threads, ifs);
   ifs.close();

   bool finished;

   int pos = 0 ;

   ofstream ofs;

   StopWatch clock;
   clock.start();
   

   double next_stop = refresh;



#pragma omp parallel shared(next_stop, refresh, arr, k_size, table, query_fn, ofbase)  private(ifs, ri_len, finished, pos, ofs, ofname, line)

   {

     finished = false;
     ri_len = -1;

     ifs.open(query_fn.c_str());
     
     ofname = ofbase;

     std::stringstream outs;

     outs << omp_get_thread_num();
     ofname += outs.str();

     ofname += ".out" ;

     ofs.open(ofname.c_str());

     ifs.seekg(arr[omp_get_thread_num()]);

     ifs >> line;

     while (!finished)   {

       //#pragma omp critical    
       //cout << omp_get_thread_num()  << " - pos:" << ifs.tellg() << "\n"; 

       if (line[0] != '>') {
  	 
	 if( ri_len == -1) {
	   // if not specified use first read as the fixed length
	   ri_len = line.length();
	 } 

	 //#pragma omp critical
	 ofs<<"read: "<<line << " pos: " << ifs.tellg() <<  endl;

       //split read into multiples of "ri_len" (read interval length)
	 proc_line(ri_len, line, k_size, table, ofs);


       }
       ifs >> line;

       pos = ifs.tellg();

       if ((pos >= arr[1+omp_get_thread_num()]) || (pos == -1))   
	 finished = true;

#ifdef HAVE_ZAP_REGION       
       if (refresh && omp_get_thread_num() == 0 && clock.queryElapsedTime() >=next_stop) {
	 next_stop += refresh;
	 zap_region();

       }
#endif

     }
     ofs.close();

   }

   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
