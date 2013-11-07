#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <omp.h>
#include "kencode.hpp"


#include "SortedDb.hpp"

#define USE_SORTED_DB 1


#include "TaxNodeStat.hpp"
#include "StopWatch.hpp"

#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;


#define KMER_IN_MAX 134217728

size_t perm_bytes_allocd;


void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <jellylist kmerdb file> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  
}

my_map tid_rank_map;

void proc_kmer(SortedDb<uint32_t> *table, uint64_t kmer_id)
{
  TaxNodeStat<uint32_t> *h = new TaxNodeStat<uint32_t>(*table);

    h->begin(kmer_id, tid_rank_map, 1000, false);
  //  h->begin(kmer_id); //, tid_rank_map, 500, false);
   //  cout<<"debug kmer "<<kmer_id<<" "<<pos;
    cout<<h->taxidCount() << " - " ;

   while( h->next() ) {  
     const uint32_t tid = h->taxid();
     //collect the unique set of tax ids
     
     cout<<" "<<tid;
   }

       /*       else {
	      cout<<" none";
       }
       */
       cout<<endl; 
   delete h;
}

int compar(const void *x, const void *y)
{
  return ((*(size_t*)x) - (*(size_t*)y));


}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1, ri_len = -1;
   int n_threads = 0;


   bool full_lookup = true;

   float threshold = 0.0;
   string kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options;

   unsigned long read_count = 0, read_max = 0;
      
   size_t mmap_size;
   bool sort = false;

   int profile = 0;

   string rank_table_file;

   while ((c = getopt(argc, argv, "p:o i:d:t:s:q:h m:")) != -1) {
   
      switch(c) {
      case 'm':
	rank_table_file = optarg;
	break;
      case 'p':
	profile = atoi(optarg);
	break;
      case 'o':
        sort = true;
        break;
      case 'h':
	full_lookup = false;
	break;
	case 'q':
	  read_max = atoi(optarg);
	  break;
      case 's':

	mmap_size = atoi(optarg);
	mmap_size *=(1<<30);
        cout << "Input heap size: " << mmap_size << endl;
        break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'i':
         query_fn = optarg;
         break;
      case 'd':
         kmer_db_fn = optarg;
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
   SortedDb<uint32_t> *taxtable;

   if (rank_table_file.length() > 0) {

       FILE * rmfp = fopen(rank_table_file.c_str(), "r");

       uint32_t src, dest;

       while (fscanf(rmfp,"%d%d", &src, &dest) > 0) {
	 tid_rank_map[src] = dest;
       }

       fclose(rmfp);

   }







#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDB*, std::size_t> ret = mfile.find<INDEXDB>("KmerDB");

     taxtable = ret.first;
     taxtable->conv_ptrs();

#else
   perm(&taxtable, sizeof(taxtable));
   mopen(kmer_db_fn.c_str(), "r", mmap_size);
#endif
   
   k_size = taxtable->get_kmer_length();



   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;







   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */



   size_t count = 0;
   
   if (read_max == 0)
     read_max=KMER_IN_MAX;

   size_t *buf = (size_t*)malloc(sizeof(size_t)*read_max);

   string line;

   ifstream ifs(query_fn.c_str());

   StopWatch clock;

   clock.start();


   int file_num = 0;

   while (true)   {



       getline(ifs, line);

       if (ifs.tellg() == -1)
	 break;
       
       ifstream   ifs2(line.c_str());

       

       while (true) {
	 
	 getline(ifs2, line);

	 sscanf(line.c_str(), "%lld", &buf[count]);
	 count++;

       
	 if (ifs.tellg() == -1)
	   break;
	 if (count > (read_max / n_threads) * file_num)
	   break;
       }
       
       file_num++;

       if (count > read_max)
	 break;
       
   }
   
   count--;
   
   size_t i;
 
   cout << "Load time, kmers: " << clock.stop() << " " << count << endl;
   clock.reset();

   if (sort) {
     clock.start();   
  
     qsort(buf,count,sizeof(size_t), compar);
     cout << "Sort time: " << clock.stop() << endl;
     clock.reset();
     
   }


   clock.start();

   if (profile > 0) {

     uint32_t *dest_arr = new uint32_t[profile];

     taxtable->get_values(dest_arr, profile);

     int i;

     for (i=0; i < profile; i++) {

       cout << dest_arr[i] << "\n";

     }


     delete dest_arr;

   } 

   else {
#pragma omp parallel shared(taxtable, buf, count, n_threads) private(i)
   {

     for (i=omp_get_thread_num()*(count/n_threads); i < (1+omp_get_thread_num())*(count/n_threads); i++)  {
	 


       proc_kmer(taxtable, buf[i]);

	 
       
     }


   }
   }

   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
