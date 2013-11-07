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
#include <tr1/unordered_map>


#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;

size_t perm_bytes_allocd;
id_convback_map_t conv_map;

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

int proc_mer(INDEXDB<DBTID_T>* table, const string& read, const size_t mer_len, FILE *ofs, ofstream &log)
{
  bool verbose=false;
  kencode_c ken(mer_len);

  kmer_t kmer_id = 0x0;
  kmer_id = ken.kencode(read);


  TaxNodeStat<DBTID_T> *h = new TaxNodeStat<DBTID_T>(*table);
  h->begin(kmer_id, & conv_map);

  
  const uint16_t ng = h->taxidCount();

  if ( ng > 0 ) {

    assert(fwrite(&kmer_id,8, 1, ofs) == 1);
    assert(fwrite(&ng,2, 1, ofs) == 1);
    while( h->next() ) {  
      const uint32_t tid = h->taxid();
      //collect the unique set of tax ids

      assert(fwrite(&tid,4, 1, ofs) == 1);
		    

    }


    delete h;
    return 1;

  } else {
    log << kmer_id << endl;
    delete h;
    return 0;
  }
  

}
 

void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <tax_histo file> -i <query fasta file> -t <number of threads> -o <output filename base> -k <kmer_len> [-r]\n";
  cout << "specify -r if using mmap DB\n";
  cout << "-k <kmer_len> does not need to be specified if using mmap DB\n";
  //cout << execname << " -d <jellylist kmerdb file> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  
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
   int k_size=-1;
   int n_threads = 0;
#if WITH_PJMALLOC == 1
   bool restore = false;
#endif 

   bool ascii = false;

   string kmer_db_fn, query_fn, ofname, ofbase, id_bit_conv_fn;


   unsigned long read_count = 0,  max_reads = 0;
   

   while ((c = getopt(argc, argv, "o:t:i:d:q:f:")) != -1) {
      switch(c) {
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
        omp_set_num_threads(n_threads);
        break;
      case 'o':
	ofbase = optarg;
	break;	
      case 'f':
        id_bit_conv_fn = optarg;
        break;

      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }

   if (kmer_db_fn == "" || query_fn == "" || n_threads == 0)  {
     usage(argv[0]);
     return -1;

   }


   cout << "Start kmer DB load\n";
   //KmerDB *table;
   INDEXDB<DBTID_T> *taxtable;


   perm(&taxtable, sizeof(taxtable));
   mopen(kmer_db_fn.c_str(), "r", MMAP_SIZE);
   
   k_size = taxtable->get_kmer_length();

   cout << "DB size is " << taxtable->size() << endl;

   assert(k_size > 0 );

   if (id_bit_conv_fn.length() > 0) {
     cout << "Loading map file,\n";
     FILE * tfp = fopen(id_bit_conv_fn.c_str(), "r");

     uint32_t src;
     uint16_t dest;

     while (fscanf(tfp,"%d%hd", &src, &dest) > 0) {

       conv_map[dest] = src;
     }
     fclose(tfp);
   }



   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   size_t *arr = split_file(n_threads, ifs);
   ifs.close();

   bool finished;

   size_t pos = 0 ;


     //unsigned cnt = 0;
   finished = false;

   FILE * out_file;
  
   StopWatch clock;
   clock.start();

   KmerFileMetaData metadata;  

#pragma omp parallel shared(arr, k_size, query_fn, ofbase, taxtable,  max_reads, clock)  private(ifs, finished, pos, out_file, ofname, line, read_count, metadata)

   {

     read_count = 0;

     finished = false;

     //     cout<<"Read query file: "<<query_fn<<endl;
     ifs.open(query_fn.c_str());
     if(!ifs) {
         cerr<<"did not open for reading: "<<query_fn<<endl;
         exit(-1);
     }
     
     ofname = ofbase;

     std::stringstream outs;

     ofstream logfs;

     outs << omp_get_thread_num();
     ofname += outs.str();

     ofname += ".out" ;


     cerr << "opening for writing: " << ofname.c_str() << endl;
     FILE *out_file = fopen(ofname.c_str(), "wb");
     assert(out_file);



     //write metadata
     metadata.setVersion(TAX_HISTO_VERSION);
     metadata.setSize(0); // try for now?
     metadata.setKmerLength(k_size);
     metadata.setDefaultDataStart();

     metadata.write(out_file);


     ofname += ".log" ;

     logfs.open(ofname.c_str());

     ifs.seekg(arr[omp_get_thread_num()]);

     string read_buff;

     uint64_t sanity = ~0;


     while (!finished)   {

       getline(ifs, line);


       pos = ifs.tellg();

       if ((pos >= arr[1+omp_get_thread_num()]) || (pos == -1)) {
	      finished = true;
       } 

       //#pragma omp critical    
       //cout << omp_get_thread_num()  << " - pos:" << ifs.tellg() << "\n"; 

       //if (line[0] != '>' && line.length() > 1) {



            // if not specified use first read as the fixed length

	       //#pragma omp critical

       //split read into multiples of "ri_len" (read interval length)
       int written = 0;

       if (line.length() > 0) {

	 written = proc_mer(taxtable, line, k_size, out_file, logfs);
	 read_count+=written;
       }
	 //	 if (read_count == max_reads)
	 //  finished = true;

       if (written &&  (read_count % TAX_HISTO_SANITY_COUNT == 0)) {
	 assert(fwrite(&sanity, 8, 1, out_file) == 1);
	 
       }

       
     }

     fseek(out_file, 0, SEEK_SET);

     metadata.setSize(read_count);
     metadata.write(out_file);
     fclose(out_file);

   }

  
   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
