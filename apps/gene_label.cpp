#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <omp.h>
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"
#include <gzstream.h>
#include <cmath>


#define MMAP_SIZE 0
#define TID_T uint32_t


#define TAXNODE TaxNode<uint32_t>
#define TAXNODESTAT TaxNodeStat<uint32_t>
#define TAXTREE TaxTree<uint32_t>
#define INDEXDBSZ INDEXDB<uint32_t>

//#define TAXNODE TaxNode
//define TAXNODESTAT TaxNodeStat
//define TAXTREE TaxTree
//define INDEXDBSZ INDEXDB


using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

using namespace metag;

static bool verbose=false;

bool add_root_on_kmer_drop = true;

size_t perm_bytes_allocd;

typedef uint32_t taxid_t;
typedef uint32_t geneid_t;
typedef pair<uint32_t,uint32_t> upair_t;
typedef pair<uint32_t,float> ufpair_t;
typedef pair<uint32_t,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<int16_t,tax_data_t> label_info_t;
typedef map<uint32_t,uint32_t> map_t;
typedef map<uint32_t,float> ufmap_t;
typedef map<uint16_t, ufmap_t> u_ufmap_t;
typedef map<uint32_t,uint32_t> hmap_t;

vector <int> read_len_vec;
vector <int> read_len_avgs;

#define USE_KPATH_IDS 0


bool badGenomes(uint32_t tid) {
   bool isBad=false;
   switch(tid) {
#if USE_KPATH_IDS == 1
      case 373758:
      case 452244:
      case 383216:
      case 814203:
      case 346687:
      //case 44004: //neanderthal (?)
      case 435171:
         isBad=true;
         break;
#endif
   }
   return isBad;
}


typedef std::pair<TAXNODE*,uint16_t> tncpair_t;

struct CmpDepth {
   CmpDepth(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<uint32_t,float>& a, const pair<uint32_t,float>& b) const {
      const int adepth = (*_imap.find(a.first)).second;
      const int bdepth = (*_imap.find(b.first)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

struct Cmp{
   bool operator()(const upair_t& a, const upair_t b) const {
      return a.second > b.second;
   }
};


struct CmpDepth1 {
   CmpDepth1(const hmap_t& imap) : _imap(imap) {}
   bool operator()(uint32_t a, uint32_t b) const {
      const int adepth = (*_imap.find(a)).second;
      const int bdepth = (*_imap.find(b)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

static bool isHuman(uint32_t taxid)  {
   bool res = false;
   switch (taxid) {
#if _USE_KPATH_IDS == 0
   case 10000348:
   case 10000349:
   case 10000350:
   case 10000351:
   case 10000352:
   case 10000353:
   case 10000354:
   case 10000355:
   case 10000356:
   case 10000357:
   case 10000358:
   case 10000359:
   case 10000360:
   case 10000361:
   case 10000362:
   case 10000363:
   case 10000364:
   case 10000365:
   case 10000366:
   case 10000367:
   case 10000368:
   case 10000369:
   case 10000370:
   case 10000371:
   case 10000372:
   case 10000373:
   case 9606:
   case 63221: //neanderthal:
#else
   case 8122:
   case 44004:
   case 8121:
   case 177579:
   case 8120:
#endif
      res=true;
      break;
   default:
      res=false;
      break;
   }
   return res;
}

void
fill_in_labels(const TAXTREE& taxtree, vector<uint32_t>& row, const tax_data_t& taxids, const hmap_t& idx2taxid, const hmap_t& taxid2idx) {
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const uint32_t col_tax_val = (*mtch).second;
      const TAXTREE::const_iterator col_tax_val_it = taxtree.find(col_tax_val); 
      if(col_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
         continue;
      }
      TAXNODE* col_tax_val_node = (*col_tax_val_it).second;
      // see if this id's genome count needs to be update
      //  based on its children   - can skip if its a leaf
      if( !col_tax_val_node->isLeaf() && row[tax_idx] == 0 ) {
         tax_data_t::const_iterator it = taxids.begin();
         const tax_data_t::const_iterator is = taxids.end();
         unsigned has_cnt_sum = 0;
         for(; it != is; ++it) {   
            const uint32_t tax_id = (*it).first;
            const TAXTREE::const_iterator tax_val_it = taxtree.find(tax_id); 
            if(tax_val_it == taxtree.end()) {
               cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
               continue;
            }
            TAXNODE* tax_val_node = (*tax_val_it).second;
            // test if tax_id is a descendant of tax_val
            // save max value
            const hmap_t::const_iterator idx_mtch = taxid2idx.find(tax_id);
            assert(idx_mtch != taxid2idx.end());
            const uint32_t tax_val_idx = (*idx_mtch).second;
            
            if( col_tax_val != tax_id && row[tax_val_idx] > 0 && tax_val_node->isAncestor(col_tax_val)) {
               // want to capture the descendant tax nodes closest to col_tax_val_node 
               has_cnt_sum = 1;
               break;
            }
         }
         if( has_cnt_sum > 0 ) {
            row[tax_idx] = has_cnt_sum;
            if(verbose) cout<<"Final Tax Node Val "<<col_tax_val<<" final transfer cnt for taxid="<<row[tax_idx]<<endl;
         }
      }
   }
   if(verbose) cout<<"fill_row: "; 
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const uint32_t tax_val = (*mtch).second;
      if(verbose) cout<<"["<<tax_val<<","<<row[tax_idx]<<"]";
   }
   if(verbose) cout<<endl;
}

struct TCmp { 
   TCmp(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<uint32_t,float>& a, const pair<uint32_t,float>& b) const {
	   if( fabs(a.second- b.second) < 0.001 ) {
		   const int adepth = (*_imap.find(a.first)).second;
		   const int bdepth = (*_imap.find(b.first)).second;
		   return adepth < bdepth;
      } 
	   return a.second < b.second; 
   }
   const hmap_t& _imap;
};

struct ScoreOptions {
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(false), _strict_kmer_match(false), _prn_all(false), _imap(imap), _diff_thresh(1.0), _diff_thresh2(3.0) {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh, _diff_thresh2;
   u_ufmap_t _rand_hits; 
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
unsigned retrieve_kmer_labels(INDEXDBSZ* table, const char* str, const int slen, const kmer_t klen, 
                          list<uint32_t>& geneid_lst, 
                          uint16_t max_count, hmap_t& gene_track) 
   {
    unsigned valid_cnt = 0;
    int j; /* position of last nucleotide in sequence */
    int k = 0; /* count of contiguous valid characters */
    int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
    kmer_t mask = ((kmer_t)1 << klen*2)-1; /* bits covering encoded k-mer */
    kmer_t forward = 0; /* forward k-mer */
    kmer_t reverse = 0; /* reverse k-mer */
    kmer_t kmer_id; /* canonical k-mer */

    set<kmer_t> no_dups;
    for (j = 0; j < slen; j++) {
        register int t;
        ENCODE(t, str[j], k);
        forward = ((forward << 2) | t) & mask;
        reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
        if (++k >= (signed)klen) {
            /* zero based position of forward k-mer is (j-klen+1) */
            /* zero based position of reverse k-mer is (slen-j-1) */
           // const int pos = j-klen+1;
           kmer_id = (forward < reverse) ? forward : reverse;
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
           no_dups.insert(kmer_id);
           ++valid_cnt; 
           //if(verbose) cout<<"debug valid: "<<j<<" "<<pos<<endl;
           TAXNODESTAT *h = new TAXNODESTAT(*table);
           //if(verbose) cout<<"lookup kmer at posit:"<<j<<endl;
           h->begin(kmer_id,NULL);
           while( h->next() ) {
              uint32_t gid = h->taxid();
              if( gene_track.find(gid) == gene_track.end() ) {
                  gene_track.insert(make_pair(gid,1));
                  geneid_lst.push_back(gid);
              } else {
                  gene_track[gid] += 1; 
              }
              //if(verbose) cout<<"found this: gene_tid="<<tid<<" label="<<assigned_taxid<<" gid="<<gid<<endl;
           }
           //if( verbose ) cout<<"num taxids: "<<geneid_lst.size();
           //if(dcnt > 0 && verbose ) cout<<" end k-mer lookup"<<endl;
           //if(mtch == 0 && verbose ) cout<<" no k-mer matches "<<endl;
           delete h;
       }
   }
   return valid_cnt;
}

void proc_line(const string &line, int k_size, INDEXDBSZ *table, ofstream &ofs, 
               const ScoreOptions& sopt, uint16_t max_count, hmap_t& gene_update, float min_score, int min_kmer) {

     const int ri_len = line.length();
     if(ri_len < 0 ) {
         cout<<"unexpected ri_len value: "<<ri_len<<endl;
         return;
     } else if( ri_len < k_size ) {
         ofs<<"\t"<<-1<<" "<<-1<<"\tNone\t"<<"\t-1 -1 ReadTooShort"<<endl;
     } else {
        vector<label_info_t> label_vec(ri_len-k_size+1,make_pair(-1,tax_data_t()));
        list<geneid_t> geneid_lst; 
        hmap_t gene_track;
        const unsigned cnt = retrieve_kmer_labels(table, line.c_str(), ri_len, k_size,geneid_lst,max_count,gene_track);
        if( !geneid_lst.empty() ) {  
           list<uint32_t>::const_iterator it = geneid_lst.begin();
           const list<uint32_t>::const_iterator is = geneid_lst.end();
           vector< upair_t > gsort(geneid_lst.size()); 
           for(unsigned i = 0; it != is; ++it, ++i) {
             upair_t st = *(gene_track.find(*it));
             gsort[i] = st;
             if(verbose) cout<<"gene scores: "<<gsort[i].first<<" "<<gsort[i].second<<endl;
           }
           sort(gsort.begin(),gsort.end(),Cmp()); 
           const float gscore = (float)gsort[0].second/(float)cnt;
           const uint32_t gl = gsort[0].first;
           ofs<<"\t"<<-1<<" "<<gsort[0].second<<" "<<cnt<<"\t"<<gl<<" "<<gscore<<" GL"<<endl;
           if( gscore > min_score && (signed)cnt > min_kmer) {
               ++gene_update[gl];
           } 
        } else {
           ofs<<"\t"<<-1<<" "<<-1<<"\tNone\t"<<"\t-1 -1 NoDbHits"<<endl;
        }
     }
  }


   size_t *split_file(int n_threads, ifstream &file)
   {
     size_t * arr = new size_t[n_threads+1];

     arr[0] = 0;

     file.seekg(0, ios::end);
  size_t end = file.tellg();

  for (size_t i = 1; (signed)i<n_threads; i++)  {

    file.seekg( i * (end / n_threads ));
    
    string junk;
    getline(file, junk);
    
    arr[i] = file.tellg();
  }

  arr[n_threads] = end;

  return arr;

}

string *split_file_names(int& n_threads, const string& file_lst)
   {
     ifstream ifs(file_lst.c_str());
     list<string> lst;
     string filename;
     while(ifs>>filename) {
         lst.push_back(filename);
     } 
     if( n_threads != 0 && n_threads != (signed)lst.size() ) {
         cout<<"warning, thread count overwritten (for now assume when a list of LMAT taxonomy classification files are given, a thread is created for each file)"<<endl;
         n_threads=lst.size();
     } else if( n_threads == 0 ) { 
         n_threads=lst.size();
     }
     string * arr = new string[n_threads+1];

     list<string>::const_iterator it = lst.begin(); 
     const list<string>::const_iterator is = lst.end(); 
      for(unsigned cnt = 0; it != is; ++it, ++cnt) {
         arr [cnt] = *it;
      }
  return arr;

}


void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:restore persistent db] [-a:ascii input format][-h:max_tid_count\n";
  cout << "note: -r makes -k unnecessary\n"; 
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int n_threads = 0, k_size = -1;
   bool restore = false, ascii = false;
   float min_score = 0.0;
   int min_kmer = 0;

   string genefile, kmer_db_fn, query_fn, query_fn_lst, ofname, ofbase;
   hmap_t imap;	
   ScoreOptions sopt(imap);
   size_t mmap_size = 0;
   uint16_t max_count = ~0;

   while ((c = getopt(argc, argv, "h:n:jb:ye:wmpk:c:v:k:i:d:l:t:s:r o:x:f:g:z:q:a ")) != -1) {
      switch(c) {
      case 'g' :
         genefile = optarg;
         break;
      case 'h' :
        max_count = atoi(optarg);
        if (max_count < 0) {
          max_count *= -1;
          add_root_on_kmer_drop = true;
        }
        break;
      case 's':
         mmap_size = atoi(optarg);
         mmap_size = mmap_size * (1<<30);
         cout << "Input heap size: " << mmap_size << endl;
         break;
      case 'j':
         verbose = true;
         break;
      case 'b':
         sopt._diff_thresh = atof(optarg);
         break;
      case 'l':
         query_fn_lst= optarg;
         break;
      case 'y':
         verbose = true;
         break;
      case 'p':
         sopt._prn_all = true;
         break;   
      case 'r':
        restore = true;
        break;
      case 'a':
        ascii = true;
        break;
      case 't':
        n_threads = atoi(optarg);
        break;
      case 'x':
         min_score = atof(optarg);
         break;
      case 'q':
         min_kmer = atoi(optarg);
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
      case 'o':
         ofbase = optarg;
         break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }
   if (ofbase == "" || kmer_db_fn == "" )  {
     cout<<"essential arguments missing: ["<<ofbase<<"] ["<<n_threads<<"] ["<<kmer_db_fn<<"] ["<<query_fn<<"] "<<endl; 
     usage(argv[0]);
     return -1;
   }
   if (!restore && k_size == -1)  {
     cout << "missing kmer size!\n";
     usage(argv[0]);
     return -1;

   }
   if(verbose) cout<<"verbose on"<<endl;
   cout << "Start kmer DB load\n";
   INDEXDBSZ *taxtable;
#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDBSZ*, std::size_t> ret = mfile.find<INDEXDBSZ>("KmerDB");

     taxtable = ret.first;
     k_size = taxtable->get_kmer_length();
     cout << "k size:  " << k_size  <<   endl ;
#else


#if WITH_PJMALLOC == 1

   if (restore) {

     perm(&taxtable, sizeof(taxtable));
     if(mopen(kmer_db_fn.c_str(), "r", mmap_size) != 0 ) {
         cout<<"Error opening db file, must exit:"<<kmer_db_fn<<endl;
         return -1;
     }

     if (k_size < 1)
       k_size = taxtable->get_kmer_length();

     cout << "num kmers: " << taxtable->size() << " - " << k_size  <<   endl ;

   } else 
    
#endif
{
#if (USE_SORTED_DB == 0)
     taxtable = new INDEXDBSZ;
     ifstream qifs(kmer_db_fn.c_str());
     if( !qifs ) {
       cerr<<"Unable to open: "<<kmer_db_fn<<endl;
       return -1;
     }
     string fname;
     
     while(qifs>>fname) {
       cout<<"register file: "<<fname<<endl;
       taxtable->registerFile(fname.c_str());
     }
     taxtable->ingest(ascii);
#endif
   }

#endif
   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;

   assert(k_size > 0 );

   size_t * arr= NULL;
   string * file_lst=NULL;
   ifstream ifs;
   if( query_fn.length() > 0 ) {
      /* ASSUME READ FITS ON ONE LINE FOR NOW!! */
      ifs.open(query_fn.c_str());
      arr = split_file(n_threads, ifs);
      ifs.close();
   } else if( query_fn_lst.length() > 0 ) {
      file_lst = split_file_names(n_threads,query_fn_lst);
      cout<<"set threads="<<n_threads<<endl;
   }
   omp_set_num_threads(n_threads);
   vector< map<TID_T,hmap_t> > gtrackall(n_threads);
   string line;
   bool finished;
   size_t pos = 0 ;
   ofstream ofs;

   StopWatch clock;
   clock.start();
   size_t read_count = 0; 

#pragma omp parallel shared(arr, file_lst, k_size, query_fn, ofbase, taxtable, sopt, gtrackall, min_score, min_kmer)  private(ifs,finished, pos, ofs, ofname, line, read_count)

   {
     read_count = 0;
     finished = false;
     const char* fn = NULL;
     if( arr ) {
         fn = query_fn.c_str();
     } else if( file_lst ) {
         fn = file_lst[ omp_get_thread_num() ].c_str();
     }  
     //cout<<"did not open for reading: ["<<fn<<"] tid: ["<<omp_get_thread_num()<<"]"<<endl;
     ifs.open(fn);
     if(!ifs) {
         cerr<<"did not open for reading: ["<<fn<<"] tid: ["<<omp_get_thread_num()<<"]"<<endl;
         exit(-1);
     }
         
        
     ofname = ofbase;

     std::stringstream outs;

     outs << omp_get_thread_num();
     ofname += outs.str();

     ofname += ".out" ;

     ofs.open(ofname.c_str());

      
     if( arr ) {
      ifs.seekg(arr[omp_get_thread_num()]);
     }

     while (!finished)   {


       getline(ifs, line);

       pos = ifs.tellg();
       if ((signed)pos == -1) {
          finished = true;   

       }
       if( arr ) {
          if ((pos >= arr[1+omp_get_thread_num()]) ) {
            finished = true;
          } 
       }

       size_t p1 = line.find('\t');
       string hdr = line.substr(0,p1);

       size_t p2 = line.find('\t',p1+1);
       const string read_buff = line.substr(p1+1,p2-p1-1);

       size_t p3 = line.find('\t',p2+1);
       string ignore_stats = line.substr(p2+1,p3-p2-1);
       size_t p4 = line.find('\t',p3+1);
       string ignore_alt_scores = line.substr(p3+1,p4-p3-1);
       size_t p5 = line.find('\t',p4+1);
       string taxid_w_scores = line.substr(p4+1,p5-p4);

       istringstream istrm(taxid_w_scores.c_str());
       float score;
       taxid_t taxid;
       string match_type;
       
       istrm >>taxid>>score>>match_type; 

       //skip human genes for now
       if(isHuman(taxid)) continue;
       //if(score < tax_threshold) continue;
       ofs<<hdr<<"\t"<<read_buff<<"\t"<<taxid<<"\t";
      
       map<TID_T,hmap_t>& track = gtrackall[omp_get_thread_num()];
       hmap_t& gtrack = track[taxid];
	    proc_line(read_buff, k_size, taxtable, ofs, sopt, max_count, gtrack, min_score, min_kmer);
	    read_count ++;
     }
     ofs.close();
   }
   typedef map<TID_T,uint32_t> tmap_t; 
   map<uint32_t,tmap_t>  merge_cnt;
   for(unsigned thread = 0; thread < gtrackall.size();  ++thread) {
      map<TID_T,hmap_t>& gt = gtrackall[thread];
      map<TID_T,hmap_t>::const_iterator it = gt.begin();
      const map<TID_T,hmap_t>::const_iterator is = gt.end();
      for(; it != is; ++it) {
         TID_T tid = (*it).first;
         const hmap_t& hm = (*it).second;
         hmap_t::const_iterator it1 = hm.begin();
         const hmap_t::const_iterator is1 = hm.end();
         for(; it1 != is1; ++it1) {
            uint32_t gid = (*it1).first;
            uint32_t cnt = (*it1).second;
            if( merge_cnt.find(gid) == merge_cnt.end() ) {
               merge_cnt.insert( make_pair(gid,tmap_t()));
            }
            if( merge_cnt[gid].find(tid) == merge_cnt[gid].end() ) {
               merge_cnt[gid][tid] = cnt;
            } else {
               merge_cnt[gid][tid] += cnt;
            }
         }
      }
   }
   igzstream zipifs(genefile.c_str());
   if( !zipifs ) {
      cout<<"no way"<<genefile<<endl;
      return -1;
   }
   ostringstream output;
   output<<ofbase<<"."<<min_score<<"."<<min_kmer<<".genesummary";
   ofstream sum_ofs(output.str().c_str());
   if( !sum_ofs ) {
      cerr<<"Can't write to "<<output.str()<<endl;
      return -1;
   }
   const unsigned buff_size = 20000;
   char buff[buff_size]; 
   while(zipifs.getline(buff,buff_size)) {
      istringstream istrm(buff);
      TID_T tid;
      uint32_t gid;
      istrm>>tid>>gid;
      if( merge_cnt.find(gid) != merge_cnt.end() ) {
         tmap_t::const_iterator ti = merge_cnt[gid].begin();
         const tmap_t::const_iterator ts = merge_cnt[gid].end();
         for(; ti != ts; ++ti) {
            TID_T label = (*ti).first;
            uint32_t cnt = (*ti).second;
            sum_ofs<<cnt<<"\t"<<label<<"\t"<<buff<<endl;
         }
      } 
   }
   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
