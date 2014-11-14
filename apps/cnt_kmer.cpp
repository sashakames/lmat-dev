#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <list>
#include <stdlib.h>
#include <omp.h>
#include <queue>
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"
#include "tid_checks.hpp"
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <gzstream.h>

#include <version.h>

#define MMAP_SIZE 0
#define TID_T uint32_t

#define HUMAN_TAXID 9606


#define QUEUE_SIZE_MAX 2000000000


using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;
using namespace metag;

static bool verbose=false;
size_t perm_bytes_allocd;
static map<TID_T,string> gRank_table;

typedef map<TID_T,unsigned> cmap_t;
typedef pair<TID_T,float> ufpair_t;
typedef pair<TID_T,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<int16_t,tax_data_t> label_info_t;
typedef map<TID_T,TID_T> hmap_t;
typedef std::tr1::unordered_map<TID_T,float> ufmap_t;
typedef std::tr1::unordered_map<TID_T,vector<float> > uvfmap_t;
typedef std::tr1::unordered_map<TID_T,string> usmap_t;
typedef std::tr1::unordered_map<uint16_t, uvfmap_t> u_ufmap_t;
typedef std::tr1::unordered_map<uint16_t, usmap_t> u_usmap_t;
typedef pair<string,string> read_pair;

#define _USE_KPATH_IDS 0



struct ScoreOptions {
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(false), _strict_kmer_match(false), _prn_all(false), _imap(imap), _diff_thresh(1.0), _diff_thresh2(3.0) {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh, _diff_thresh2;
   u_ufmap_t _rand_hits; 
   u_usmap_t _rand_class; 
   bool _comp_rand_hits;
};


#define ENCODE(t, c, k) \
switch (c) { \
case 'a': case 'A': t = 0; break; \
case 'c': case 'C': t = 1; break; \
case 'g': case 'G': t =2; break; \
case 't': case 'T': t = 3; break; \
default: k = 0; continue; \
}

/*
   Given a sequence of nucleotide characters,
   break it into canonical k-mers in one pass.
   Nucleotides are encoded with two bits in
   the k-mer. Any k-mers with ambiguous characters
   are skipped.
   str:  DNA sequence (read)
   slen: DNA sequence length in nucleotides
   klen: k-mer length in nucleotides
*/

static
int retrieve_kmer_labels(const char* str, const int slen, const kmer_t klen) {
    int j; /* position of last nucleotide in sequence */
    int k = 0; /* count of contiguous valid characters */
    int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
    kmer_t mask = ((kmer_t)1 << klen*2)-1; /* bits covering encoded k-mer */
    kmer_t forward = 0; /* forward k-mer */
    kmer_t reverse = 0; /* reverse k-mer */
    kmer_t kmer_id; /* canonical k-mer */
    set<kmer_t> no_dups;
    unsigned valid_kmers=0; 
    for (j = 0; j < slen; j++) {
        register int t;
        const char base=str[j];
        ENCODE(t, base, k);
        forward = ((forward << 2) | t) & mask;
        reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
        if (++k >= (signed)klen) {
           kmer_id = (forward < reverse) ? forward : reverse;
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
            /* zero based position of forward k-mer is (j-klen+1) */
            /* zero based position of reverse k-mer is (slen-j-1) */
           //const int pos = j-klen+1;
           /* kmer_lookup(kmer); do k-mer lookup here... */
           no_dups.insert(kmer_id);
           valid_kmers++;
       }
   }
   return valid_kmers;
}

int proc_line(int ri_len, string &line, int k_size, ofstream &ofs, const ScoreOptions& sopt) {
      int valid_kmers=-1;
     if(ri_len < 0 || ri_len > (signed)line.length()) {
         cout<<"unexpected ri_len value: "<<ri_len<<endl;
     } else if( ri_len < k_size ) {
     } else {
        valid_kmers = retrieve_kmer_labels(line.c_str(), ri_len, k_size);
     }
     return valid_kmers;
}



void usage(char *execname)
{
  cout << "LMAT version " << LMAT_VERSION  << "\n";
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l <human bias>]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:rank table]\n";
  cout << "[-m <rank/tid-map-file>] [-g <tid-cutoff>] [-w <with-strain-species-map> (affects -m option)]\n";
  cout << "[-h:turn phiX screening off]\n";
 
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=20;
   int n_threads = 0;

   string query_fn,ofbase,ofname;
   //string rank_map_file,rank_ids, kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, depth_file, rand_hits_file, rank_table_file, id_bit_conv_fn;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   bool fastq=false; 

   while ((c = getopt(argc, argv, "u:ahn:j:b:ye:w:pk:c:v:k:i:d:l:t:r:sm:o:x:f:g:z:qV")) != -1) {

      switch(c) {
      case 'y':
         verbose = true;
         break;
      case 'q':
         fastq=true;
         break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 'o':
         ofbase = optarg;
         break;
    case 'V':
      cout << "kmer cnt version " << LMAT_VERSION  << "\n";
	exit(0);
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }
   if (ofbase == "") cout << "ofbase\n";
   if (n_threads == 0) cout << "n_threads\n";
   if (query_fn == "") cout << "query_fn\n";

   string line;

   omp_lock_t buffer_lock;
   std::queue < read_pair > read_buffer_q;
   omp_init_lock(&buffer_lock);

   bool finished;

   ofstream ofs;

   StopWatch clock;
   clock.start();
   size_t read_count_in =0;
   size_t read_count_out = 0;
   string read_buff, hdr_buff, save_hdr;

   bool in_finished = false;

   ifstream tmpstream;

   istream ifs(cin.rdbuf());
   if (query_fn != "-") {
     tmpstream.open(query_fn.c_str());

     if(!tmpstream) {
	  cerr<<"did not open for reading: "<<query_fn<<endl;
	  
	  exit(-1);
	  
     }  
     ifs.rdbuf(tmpstream.rdbuf());
   }

   




#pragma omp parallel shared(k_size, query_fn, ofbase, sopt,  in_finished, read_count_in, read_count_out, ifs)  private(finished, ofs, ofname, line, read_buff, hdr_buff, save_hdr)
  {

    finished = false;

    ofname = ofbase;
    std::stringstream outs;
    outs << omp_get_thread_num();
    ofname += outs.str();
    ofname += ".out" ;
    ofs.open(ofname.c_str());
    bool eof = false;
    while (!finished)   {
      if ((in_finished == false) && (omp_get_thread_num() == 0)) {
	      int j = 0 ;
	      int queue_size = 0;
	      omp_set_lock(&buffer_lock);
	      queue_size = read_buffer_q.size();
	      omp_unset_lock(&buffer_lock);
	      string last_hdr_buff;
	      while (queue_size < QUEUE_SIZE_MAX && j< 2* n_threads && (!in_finished)) {
            eof = !getline(ifs, line);
   	      if (eof) {
	            in_finished = true;
	    
		    if(verbose) cout << line.size() << " line length\n";
		    line = "";
	  }


	  if (line[0] == '>' || (fastq && line[0] == '@') ) {

	    last_hdr_buff = hdr_buff;
	    // skip the ">"                                                        
	    hdr_buff=line.substr(1,line.length()-1);
	  }

	  if (line[0] != '>' && line.length() > 1 && !fastq) {
	    read_buff += line;
	    line = "";
	  }

	  if( fastq && line[0] != '@' && line[0] != '+' && line[0] != '-' ) {
	    read_buff += line;
	    line = "";
	  }
	  if( ((line[0] == '>' || in_finished) || (fastq && (line[0] == '+' ||
	     line[0] == '-'))) && read_buff.length() > 0 ) {

	    omp_set_lock(&buffer_lock);
	    
	    if (in_finished)
		read_buffer_q.push(read_pair(read_buff, hdr_buff));
	    else
		read_buffer_q.push(read_pair(read_buff, last_hdr_buff));

	    read_count_in++;
	    omp_unset_lock(&buffer_lock);

	    read_buff="";

	    j ++;
	    
	    if(fastq) eof = !getline(ifs, line); // skip quality values for now       

	  }

	  if (in_finished) {
	    
	    cout << read_count_in << " reads in\n";
	    break;

	  }
	}
	
      }
      read_buff = "";
      save_hdr = "";
      omp_set_lock(&buffer_lock);
      if (!read_buffer_q.empty()) {
         read_pair in_pair = read_buffer_q.front();
         read_buff = in_pair.first;
         save_hdr = in_pair.second;
         read_buffer_q.pop();
         read_count_out++;

      }

      omp_unset_lock(&buffer_lock);
      if (read_buff.length() > 0) {
        if(save_hdr[0] == '\0') {
           ostringstream ostrm;
           ostrm<<"unknown_hdr:"<<read_count_out;
           save_hdr=ostrm.str();
        }
        int result = proc_line(read_buff.length(), read_buff, k_size, ofs, sopt); 
        ofs<<save_hdr<<"\t"<<result<<endl; 
	     read_buff="";
	
      }
      if ((read_count_in == read_count_out) && in_finished)
	   finished = true;
    }
    ofs.close();
  }
   return 0; 
}
