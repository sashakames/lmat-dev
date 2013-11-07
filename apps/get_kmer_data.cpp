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


#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;


size_t perm_bytes_allocd;

#define KMER_IN_MAX 134217728

size_t perm_bytes_allocd;


void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <jellylist kmerdb file> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  
}


void proc_kmer(TaxTable *table, uint64_t kmer_id)
{
   TaxNodeStat *h = new TaxNodeStat(*table);

   h->begin(kmer_id);
   //  cout<<"debug kmer "<<kmer_id<<" "<<pos;
   const uint16_t ng = h->genomeCount();
       if ( ng > 0 ) {
	 //	         cout<<" "<<ng<<" "<<h->taxidCount();
	         while( h->next() ) {  
	            const uint32_t tid = h->taxid();
	            //collect the unique set of tax ids
	            const uint16_t pr_cnt = h->present();
		    //         cout<<" "<<tid<<" "<<pr_cnt;
	         }
       } 
       /*       else {
	      cout<<" none";
       }
       cout<<endl; */
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

   while ((c = getopt(argc, argv, "i:d:t:s:q:h o ")) != -1) {
   
      switch(c) {
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
   TaxTable *taxtable;







#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<TaxTable*, std::size_t> ret = mfile.find<TaxTable>("KmerDB");

     taxtable = ret.first;
#else
   perm(&taxtable, sizeof(taxtable));
   mopen(kmer_db_fn.c_str(), "r", mmap_size);
#endif
   
   k_size = taxtable->get_kmer_length();

   cout << taxtable->size() << " " << taxtable->bucket_count() << endl;



   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;




   ifstream ifs(query_fn.c_str());


   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */



   size_t count = 0;
   
   if (read_max == 0)
     read_max=KMER_IN_MAX;

   size_t *buf = (size_t*)malloc(sizeof(size_t)*read_max);

   string line;


   while (true)   {
       getline(ifs, line);
       sscanf(line.c_str(), "%lld", &buf[count]);
       count++;
       if (count > read_max)
	 break;
       
       if (ifs.tellg() == -1)
	 break;


   }

   StopWatch clock;


   if (sort) {
     clock.start();   
  
     qsort(buf,count,sizeof(size_t), compar);
     cout << "Sort time: " << clock.stop() << endl;
     clock.reset();
     
   }

   
   size_t i;
 

   cout << count << "\n" ;


   clock.start();   



#pragma omp parallel shared(taxtable, buf, count, n_threads) private(i)
   {

     for (i=omp_get_thread_num()*(count/n_threads); i < (1+omp_get_thread_num())*(count/n_threads); i++)  {
	 

       if (full_lookup)
	 proc_kmer(taxtable, buf[i]);
       else 
	 bool foo = taxtable->find(buf[i]) != taxtable->end();
	 
       
     }


   }

   cout << "throughput: " << double(read_max) / clock.stop() << endl;

   return 0; 
}
