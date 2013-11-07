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

void retrieve_kmer_labels(const string& read, const size_t mer_len, bool revStrnd)
{
  bool verbose=false;
  kencode_c ken(mer_len);
  size_t read_len = read.length();
  uint64_t kmer_id = 0x0;
  bool bad_last_kmer=true;
  unsigned next_valid_pos = 0;
  
  if (read.size() < mer_len)
    return;

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

      cout << kmer_id << endl;
    }
  }
}

void proc_line(string &rval, int k_size) {



       if (rval.size() < k_size)
	 return;

       retrieve_kmer_labels(rval, k_size, false);

       string rev_cmp(rval.size(),'\0');
 


       for (int j=(rval.size()-1); j >= 0; --j) {
         const int bidx = (rval.size()- 1) - j;
         rev_cmp[j] = getComp(rval[bidx]);
       }
       retrieve_kmer_labels(rev_cmp, k_size, true);

    }
   
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

   bool ascii = false;

   float threshold = 0.0;
   string kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options;


   unsigned long read_count = 0, read_max = 0;
   

   while ((c = getopt(argc, argv, "k:i:l:q:")) != -1) {
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
      case 'q':
	read_max = atoi(optarg);
	break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }

   if (k_size < 1 || query_fn == "")  {
     usage(argv[0]);
     return -1;

   }





   assert(k_size > 0 );
   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */


   bool finished;

   int pos;


     //unsigned cnt = 0;
   finished = false;
     
   getline(ifs, line);


   while (!finished)   {


       if (line[0] != '>') {
  	 
	       //          cout<<"read: "<<line << " pos: " << ifs.tellg() <<  endl;

	       proc_line(line, k_size);
	       read_count++;
	       if (read_count == read_max)
	         break;
       }

       getline(ifs, line);


       if( line[0] =='\0') {
         finished=true;
       }
   }
     

   return 0; 
}
