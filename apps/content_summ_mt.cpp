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
#include "tid_checks.hpp"
#include "all_headers.hpp"
#include <gzstream.h>

#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

using namespace metag;

static bool verbose=false;

bool add_root_on_kmer_drop = true;

bool VERBOSE=false;

size_t perm_bytes_allocd;

typedef pair<TID_T,signed> tid_call_t;
typedef pair<TID_T,float> ufpair_t;
typedef map<TID_T,TID_T> hmap_t;
typedef map<TID_T,float> ufmap_t;
typedef map<kmer_t,uint64_t> kmer_cnt_t;


#define _USE_KPATH_IDS 0

static float compReadCnt(const vector<map<TID_T,float> >& track, TID_T tid);
static float compKmerCov(const vector<map<TID_T,kmer_cnt_t> >&, TID_T, map<TID_T,float>&, kmer_cnt_t&) ;

// as long cont_lst is not too big (which it shouldn't be) let's just do brute force for now
static bool checkContLst( const string& buff_str, const list<string>& cont_lst ) {
   bool isCont=false;
   list<string>::const_iterator it = cont_lst.begin();
   const list<string>::const_iterator is = cont_lst.end();
   for(; it != is; ++it) {
      const string& cont_id = *it;
      //cout<<"needle="<<cont_id<<" haystack="<<buff_str<<endl;
      if( buff_str.find( cont_id ) != string::npos  ) {
         if(VERBOSE) cout<<"identified as genome with possible contaminants "<<cont_id<<endl;
         isCont=true;  
         break;
      }
   }
   return isCont;
}

struct QCmp {
   bool operator()(const tid_call_t& a, const tid_call_t& b) {
      return (a.second > b.second);
   }
};

//my_map tid_rank_map;

// this determines which style of tid reduction we want, default is
// the "comphrehensive reduction; requires the tid to rank mapping
// table.  -w option enables strain to species only reduction.
bool tid_map_is_strain_species = false;

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
int retrieve_kmer_labels(/*INDEXDB<TID_T>* table,*/ const char* str, const int slen, const kmer_t klen, const list<TID_T>& cand_lst, map< TID_T, 
                         kmer_cnt_t >& kmer_track) {
    int j; /* position of last nucleotide in sequence */
    int k = 0; /* count of contiguous valid characters */
    int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
    kmer_t mask = ((kmer_t)1 << klen*2)-1; /* bits covering encoded k-mer */
    kmer_t forward = 0; /* forward k-mer */
    kmer_t reverse = 0; /* reverse k-mer */
    kmer_t kmer_id; /* canonical k-mer */
    set<kmer_t> no_dups;
    int valid_kmers=0;
    for (j = 0; j < slen; j++) {
        register int t;
        ENCODE(t, str[j], k);
        forward = ((forward << 2) | t) & mask;
        reverse = ((kmer_t)(t^3) << highbits) | (reverse >> 2);
        if (++k >= (signed)klen) {
           valid_kmers++;
           kmer_id = (forward < reverse) ? forward : reverse;
            /* zero based position of forward k-mer is (j-klen+1) */
            /* zero based position of reverse k-mer is (slen-j-1) */
           //const int pos = j-klen+1;
           /* kmer_lookup(kmer); do k-mer lookup here... */
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
           no_dups.insert(kmer_id);
           const TID_T& taxid = cand_lst.front();
           if( kmer_track[taxid].find(kmer_id) == kmer_track[taxid].end() ) {
               kmer_track[taxid].insert(make_pair(kmer_id,1));
           } else { 
               kmer_track[taxid][kmer_id] += 1;
           } 
         }
    }
    return valid_kmers;
}

pair<TID_T,float> proc_line(const TaxTree<TID_T> & tax_tree, const string &read_str, const string& seq_buff, int k_size, const vector< set<TID_T> >& cand_call_order, const map<TID_T,TID_T>&  call_lin,
               map< TID_T, kmer_cnt_t >& kmer_track, TID_T call, float call_score) {
   
   list<TID_T> cand_lst;
   TID_T used_tid = call;
   istringstream istrm(seq_buff.c_str());
   TID_T cand_tid;
   float match_score, fnd_score = call_score;
   bool fnd = false;
   vector<pair<TID_T,float> > save;
   bool debug=false;
   list<TID_T> plas_lst;
   while( istrm>>cand_tid>>match_score ) {
      if( fnd && match_score < fnd_score ) {
         break;
      }
      if( call_lin.find(cand_tid) != call_lin.end() ) {
         TID_T ctid = (*call_lin.find(cand_tid)).second;
         ostringstream os1, os2;
         os1<<" "<<ctid<<" ";
         string spos1 = seq_buff.find( os.str().c_str() );
         os2<<ctid<<" ";
         string spos2 = seq_buff.find( os.str().c_str() );
         // check that the actually org call is somewhere in the list of candidates
         if( spos1 != string::npos || (spos1 != string::npos && spos1 == 0 ) ) {
         fnd=true;
         fnd_score = match_score;
         bool skipThisChr=false;
         if( ctid >= 10000000 ) {
            plas_lst.push_back(ctid);
         } else if( !plas_lst.empty() ) {
         // if we have plasmids, let the tie score go to the plasmid
         // don't store the parent. By default the parent will always
         // have more reads and be chosen instead
         // more work may eventually be needed to validate reads that are reported 
         // there is a worst case now where a read could be from a plasmid and a chromosome
         // there's no good way to differentiate these now, they're just assigned to the plasmid.
            list<TID_T>::const_iterator ti = plas_lst.begin(); 
            const list<TID_T>::const_iterator is = plas_lst.end(); 
            for(; ti != is; ++ti) {
               const TID_T pid = *ti; 
               vector<TID_T> ptor;
               TaxTree<TID_T>& tref = const_cast<TaxTree<TID_T>&>(tax_tree); 
               tref.getPathToRoot(pid,ptor);
               for(unsigned i = 0; i < ptor.size(); ++i) {
                  if( ptor[i] == ctid ) {
                     skipThisChr=true;
                     break;
                  }
               }
               if( skipThisChr ) break;
            }
         }
         if( !skipThisChr ) {
            if(debug) cout<<"debug entry1 "<<ctid<<" "<<match_score<<" "<<cand_tid<<endl;
            save.push_back(make_pair(ctid,match_score));
         } 
         }
      }
   }
   if( !fnd || save.size() == 0 ) {
      //fnd_score set to 0 for now
      if(debug) cout<<"debug entry0 "<<call<<endl;
      cand_lst.push_back(call);
   } else {
      for(unsigned j = 0; j < cand_call_order.size(); ++j) {
         for(unsigned i = 0; i < save.size(); ++i) {
            if( cand_call_order[j].find( save[i].first ) != cand_call_order[j].end() ) {
               used_tid = save[i].first;
               fnd_score = save[i].second;
               if(debug) cout<<"debug entry2 "<<used_tid<<" "<<fnd_score<<endl;
               cand_lst.push_back(used_tid);
               break;
            } 
         }
         if( !cand_lst.empty() ) {
            break;
         }
      }
      assert(!cand_lst.empty());
   }
   const int ri_len = read_str.length();
   retrieve_kmer_labels(read_str.c_str(), ri_len, k_size, cand_lst, kmer_track);
   if( debug) cout<<"debug entry3 "<<used_tid<<" "<<fnd_score<<endl;
   return make_pair(used_tid,fnd_score);
}

static string *split_file_names(int& n_threads, const string& file_lst)
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


#if 0
static void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:restore persistent db] [-h:max_tid_count\n";
  cout << "[-r <rank/tid-map-file>] [-h <tid-cutoff>] [-w <with-strain-species-map> (affects -r option)]\n";
  cout << "note: -r makes -k unnecessary\n"; 
}
#endif

struct SaveRes {
   SaveRes(TID_T tid, float gc, float rat, float wrdc) : _tid(tid), _gen_copy(gc), _ratio(rat), _wrdc(wrdc) { }
   TID_T _tid;
   float _gen_copy;
   float _ratio;
   float _wrdc; // read cnt weighted by each read's individual score
};


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=17;

   float threshold = 0.0;
   map<TID_T, uint64_t> ref_kmer_cnt;
   string sep_plas_file,query_fn_lst,lmat_sum,kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, kmer_cnt_file, depth_file, rank_table_file, cont_genomes_lst;
   hmap_t imap;	
   bool skipHuman=false;
   string thresh_str;
   while ((c = getopt(argc, argv, "m:f:ah:n:jb:ye:wp:k:c:v:k:i:d:l:t:sr:o:x:f:g:z:q:")) != -1) {
      switch(c) {
      case 's':
         skipHuman = true;
         break;
      case 'p':
         sep_plas_file = optarg;
         break;
      case 'b':
         cont_genomes_lst = optarg;
         break;
      case 'r':
         rank_table_file = optarg;
         break;
      case 'y':
         verbose = true;
         break;
      case 'm':
	     kmer_cnt_file = optarg;
        break;
      case 'l':
	     lmat_sum = optarg;
        break;
      case 'v':
         threshold = atof(optarg);
         break;
      case 'c':
         tax_tree_fn = optarg;
	      break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'f':
         query_fn_lst = optarg;
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 't':
         thresh_str = optarg;
         break;
      case 'o':
         ofbase = optarg;
         break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }
   StopWatch clock;
   const unsigned buff_size = 2024;
   char buff[buff_size];
   set<TID_T> separate_plasmid;
   if( sep_plas_file.size() > 0 ) {
      ifstream ifs(sep_plas_file.c_str());
      if( !ifs ) {
         cout<<"Error reading "<<sep_plas_file<<endl;
         return -1;
      }
      while(ifs.getline(buff,buff_size)) {
         istringstream istrm(buff);
         TID_T tid;
         int kmer_cnt;
         istrm>>tid>>kmer_cnt;
         separate_plasmid.insert(tid);
         
      }
   }
   list<string> cont_lst;
   if( cont_genomes_lst.size() > 0 ) {
      ifstream ifs(cont_genomes_lst.c_str());
      if( !ifs ) {
         cout<<"Error reading "<<cont_genomes_lst<<endl;
         return -1; 
      }
      while(ifs.getline(buff,buff_size)) {
         cont_lst.push_back(buff);
      }
   }
   float min_euk, min_prok,min_vir, min_contam, min_plas, min_wrdcnt, same_strain_thresh,min_avg_wght;
   istringstream thrsh_istrm(thresh_str.c_str());
   thrsh_istrm>>min_euk>>min_prok>>min_vir>>min_contam>>min_plas>>min_wrdcnt>>same_strain_thresh>>min_avg_wght;
   cout<<"Thresh: "<<thresh_str<<endl;
   if( !kmer_cnt_file.c_str() ) {
      cout<<"Require an expected k-mer count file"<<endl;
      return 0;
   }
   cout<<"expected kmer cnt:"<<kmer_cnt_file<<endl;
   ifstream ifs2(kmer_cnt_file.c_str());
   TID_T tid;
   uint64_t cnt;
   while(ifs2.getline(buff,buff_size)) {
      istringstream istrm(buff);
      istrm>>tid>>cnt;
      ref_kmer_cnt.insert(make_pair(tid,cnt));
   }
   map<TID_T,string> rank_table;
   if( rank_table_file.size() > 0) {
      ifstream ifs1(rank_table_file.c_str());
      string rank;
      while(ifs1>>tid>>rank) {
         rank_table.insert(make_pair(tid,rank));
      }
   }
   if( k_size <= 0 ) {
      cerr<<"kmer size="<<k_size<<" must be non-zero"<<endl;
      return -1;
   }
   int n_threads=0;
   string * file_lst = split_file_names(n_threads,query_fn_lst);
   cout<<"set threads="<<n_threads<<endl;
   assert(n_threads>=1);
   omp_set_num_threads(n_threads);
   cout<<"Read taxonomy tree: "<<tax_tree_fn<<endl;
   TaxTree<TID_T> tax_tree(tax_tree_fn.c_str());
   cout<<"Done Read taxonomy tree: "<<tax_tree_fn<<endl;
   set<TID_T> call_set;
   map<TID_T,int> merge_read_cnt;
   for(unsigned sample_iter = 0; sample_iter < 2; ++sample_iter) {
      set<TID_T> curr_call_set;
      ifstream call_ifs(lmat_sum.c_str());
      if( !call_ifs ) {
         cerr<<"Failed to open "<<lmat_sum<<" must exit now"<<endl;
         return -1;
      }
      map<TID_T,signed> save_call;
      map<TID_T,string> hold_sum;
      set<TID_T> non_leaf_call;
      set<TID_T> cont_set;
      while(call_ifs.getline(buff,buff_size)) {
         string buff_str=buff;
         if( buff_str.find("\tNULL\t") == string::npos) {
            istringstream istrm(buff);
            TID_T tid;
            unsigned read_cnt = 0;
            string descrip;
            float wght_rc = 0;
            istrm>>wght_rc>>read_cnt>>tid>>descrip;
            // assumes sorted order
            if( sample_iter == 0 ) {
               if( wght_rc < min_wrdcnt ) {
                  break;
               }
               if( (wght_rc/(float)read_cnt) < min_avg_wght ) {
                  continue;
               }
               if( checkContLst( buff_str, cont_lst ) ) {
                  cont_set.insert(tid);
               }
            } else if( sample_iter == 1 ) { 
               // Second time around we start with the calls that were determined in previous round
               // technically, re-reading from file is not needed here
               // HUMAN by default is now ignored in reporting, but must still be considered
               // on the second iteration to prevent human reads from being re-assigned.
               if( call_set.find(tid) == call_set.end() && tid != HUMAN_TID ) {
                  continue;
               }
            }
            vector<TID_T> ptor;
            tax_tree.getPathToRoot(tid,ptor);
            if( ptor.size() > 0 && tid < 10000000 ) { //tid >= 10000000 are plasmids treat their parent's as leaves also
               assert (ptor[0] != tid ); // confirm tid is included here
               for(unsigned i = 0; i < ptor.size(); ++i) {
                  non_leaf_call.insert(ptor[i]);
               }
            } else if( ptor.size() <= 0 ) {
               non_leaf_call.insert(tid);
            }
            save_call.insert(make_pair(tid,read_cnt));
            if( VERBOSE ) cout<<"check candidate: "<<buff_str<<endl;
            hold_sum.insert(make_pair(tid,buff_str));
         }
      }
      list<tid_call_t> cand_lst;
      map<TID_T,signed>::const_iterator it1 = save_call.begin();
      const map<TID_T,signed>::const_iterator is1 = save_call.end();
      for(;  it1 != is1; ++it1) {
         TID_T tid = (*it1).first;
         unsigned read_cnt_sum= (*it1).second;
         if( sample_iter == 0 ) {
            if( non_leaf_call.find(tid) == non_leaf_call.end() ) {
               // without this, plasmids will get selected first
               if( tid < 10000000 ) {
                  vector<TID_T> ptor;
                  tax_tree.getPathToRoot(tid,ptor);
                  for(unsigned i = 0; i < ptor.size(); ++i) {
                     if(VERBOSE) cout<<"add to: "<<tid<<" "<<ptor[i]<<" "<<save_call[ptor[i]]<<" "<<read_cnt_sum<<endl;
                     read_cnt_sum+=save_call[ptor[i]];
                  }
               }
               if(VERBOSE) cout<<"final read cnt: "<<tid<<" "<<read_cnt_sum<<endl;
               cand_lst.push_back( make_pair(tid, read_cnt_sum) );
            }
         } else {
            assert( merge_read_cnt.find(tid) != merge_read_cnt.end());
            unsigned read_cnt_sum = merge_read_cnt[tid];
            if(VERBOSE) cout<<"step read cnt: "<<tid<<" "<<read_cnt_sum<<endl;
            cand_lst.push_back( make_pair(tid, read_cnt_sum) );
         }
      }
      vector< tid_call_t > call_vec(cand_lst.size());
      list<tid_call_t>::const_iterator cit = cand_lst.begin();
      const list<tid_call_t>::const_iterator cis = cand_lst.end();
      for(unsigned i = 0; cit != cis; ++cit, ++i) {
         const tid_call_t& val = *cit;
         assert( i < call_vec.size());
         call_vec[i] = val;
      }
      sort( call_vec.begin(), call_vec.end(), QCmp() );
      set<TID_T> leaf_calls;
      vector< set<TID_T> > cand_call_order(call_vec.size());
      map<TID_T,TID_T> save_call_lin;
      for(unsigned cvi =0; cvi < call_vec.size(); ++cvi) {
         const TID_T tid = call_vec[cvi].first;
         vector<TID_T> localptor;
         tax_tree.getPathToRoot(tid,localptor);
         if(VERBOSE) cout<<"save call lineage: "<<cvi<<" "<<tid<<endl;
         cand_call_order[cvi].insert( tid );
         save_call_lin.insert( make_pair( tid,tid) );
         leaf_calls.insert(tid);
         for(unsigned li = 0; li < localptor.size(); ++li) {
            if(save_call_lin.find(localptor[li]) == save_call_lin.end()) {
               if(VERBOSE) cout<<"Associate "<<li<<" "<<localptor[li]<<" with "<<tid<<endl;
               save_call_lin.insert( make_pair(localptor[li],tid) );
               cand_call_order[cvi].insert( localptor[li] );
            } 
         }
      }
      string line;
      bool finished = false;
      clock.start();
      size_t read_count = 0; 
      vector<map< TID_T, kmer_cnt_t> >  kmer_track(n_threads);
      vector<map<TID_T,float> > weighted_readcnt(n_threads);
      ifstream ifs;
      vector<map<TID_T,int> > read_cnts(n_threads); 
#pragma omp parallel shared(file_lst, k_size, tax_tree, cand_call_order,save_call_lin, leaf_calls, kmer_track,weighted_readcnt,read_cnts)  private(ifs,finished, ofname, line, read_count)
      {
      read_count = 0;
      finished = false;
      const signed thread = omp_get_thread_num();
      const char* fn = file_lst[ thread ].c_str();
      string in_fn, out_fn;
      ostringstream to;
      if( sample_iter > 0 ) {
         to<<fn<<".ras."<<(sample_iter-1)<<"\0";
         in_fn = to.str();
         ostringstream to2;
         to2<<fn<<".ras."<<sample_iter<<"\0";
         out_fn = to2.str();
      } else {
         in_fn = fn;
         to<<fn<<".ras."<<sample_iter<<"\0";
         out_fn = to.str();
      }
      ifs.open(in_fn.c_str());
      if(!ifs) {
         cerr<<"did not open for reading: ["<<in_fn<<"] tid: ["<<omp_get_thread_num()<<"]"<<endl;
         exit(-1);
      }
      ofstream rofs(out_fn.c_str());
      if(!rofs) {
         cerr<<"did not open for writing: ["<<out_fn<<"] tid: ["<<omp_get_thread_num()<<"]"<<endl;
         exit(-1);
      }
      while (!finished)   {
       if(!getline(ifs, line)) break;
       size_t pos = ifs.tellg();
       if ((signed)pos == -1) {
          finished = true;
       }
       const size_t p1 = line.find('\t');
       const size_t p2 = line.find('\t',p1+1);
       const size_t p3 = line.find('\t',p2+1);
       const size_t p4 = line.find('\t',p3+1);
       const size_t p5 = line.find('\t',p4+1);
       const string hdr = line.substr(0,p1-1);
       const string read_buff = line.substr(p1+1,p2-p1-1);
       const string ignore_stats = line.substr(p2+1,p3-p2-1);
       const string alt_scores = line.substr(p3+1,p4-p3-1);
       const string taxid_w_scores = line.substr(p4+1,p5-p4-1);
       // would be faster to check for a -1 not a string but read_label *may* need to be modified
       if( taxid_w_scores.find( "NoDbHits" ) != string::npos || taxid_w_scores.find( "ReadTooShort") != string::npos ) {
         continue;
       }
       istringstream istrm(taxid_w_scores.c_str());
       float score;
       TID_T taxid;
       string match_type;
       istrm >>taxid>>score>>match_type;
       if(isHuman(taxid) && skipHuman) continue;
       if(score < threshold) continue;
       pair<TID_T,float> used = proc_line(tax_tree,read_buff, alt_scores, k_size, cand_call_order, save_call_lin, kmer_track[thread], taxid, score);
       map<TID_T,int>& read_cnt = read_cnts[thread];
       if( weighted_readcnt[thread].find(used.first) == weighted_readcnt[thread].end() ) { 
         weighted_readcnt[thread].insert(make_pair(used.first,used.second));
         assert( read_cnt.find(used.first) == read_cnt.end());
         read_cnt[used.first] = 1;
       } else {
         weighted_readcnt[thread][used.first]+=used.second;
         assert( read_cnt.find(used.first) != read_cnt.end());
         read_cnt[used.first] += 1;
       }
       string  new_match_type;
       if( used.first != taxid ) {
           new_match_type="Reassigned";
       } else {
           new_match_type=match_type;
       }
       if( sample_iter == 0 ) {
         rofs<<hdr<<"\t"<<read_buff<<"\t"<<ignore_stats<<"\t"<<alt_scores<<"\t"<<used.first<<" "<<used.second<<" "<<new_match_type<<endl;
       } else {
         // on last iteration, don't print out alternative scores to save space
         rofs<<hdr<<"\t"<<read_buff<<"\t"<<ignore_stats<<"\t"<<"-1 -1"<<"\t"<<used.first<<" "<<used.second<<" "<<new_match_type<<endl;
       }
       read_count ++;
     }
   }
      map<TID_T,int> curr_merge_read_cnt;
      for(signed ti = 0; ti < n_threads; ++ti) {
         map<TID_T,int>& read_cnt = read_cnts[ti];
         map<TID_T,int>::const_iterator mit = read_cnt.begin();
         const map<TID_T,int>::const_iterator mis = read_cnt.end();
         for( ; mit != mis; ++mit) {
            const TID_T tid = (*mit).first;
            const int rc = (*mit).second;
            if( curr_merge_read_cnt.find(tid) == curr_merge_read_cnt.end() ) {
               curr_merge_read_cnt[tid] = rc;
            } else {
               curr_merge_read_cnt[tid] += rc;
            }
         }
      }
   kmer_cnt_t merge_kmer_cnt;
   list<SaveRes> save_out; 
   float tot_gen_copy=0;
   set<TID_T> seen;
   map<TID_T,float> save_cov, merge_kmer_track, save_last_ratio;
   for(unsigned i =0; i < call_vec.size(); ++i) {
         const TID_T tid = call_vec[i].first;
         if(isHuman(tid) && skipHuman) continue;
         assert(ref_kmer_cnt.find(tid) != ref_kmer_cnt.end());
         const uint64_t expect_cov = ref_kmer_cnt[tid];
         if(VERBOSE) cout<<"candidate organism "<<i<<" "<<call_vec[i].first<<" "<<call_vec[i].second<<" "<<expect_cov<<endl;
         const float targ_gen_copy = compKmerCov(kmer_track,tid,merge_kmer_track,merge_kmer_cnt);
         const float wrdc = compReadCnt(weighted_readcnt,tid);
         uint64_t k_mer_cnt = merge_kmer_cnt[tid];
         bool isEuk = false, isVir = false, isProk = false, isFrancisella = false, isBrucella = false, isEphage=false;
         if(VERBOSE) cout<<sample_iter<<" Target coverage for "<<tid<<" is "<<targ_gen_copy<<" add "<<k_mer_cnt<<" wrdc="<<wrdc<<endl;
         // for the moment only use plasmid specific reads to measure its k-mer coverage
         // There is a problem now with it pulling in k-mer coverage from higher order reads, 
         // that may belong to a chromosome, and the ratio appears inflated
         set<string> seen_rank;
         float last_sibling_ratio = -1;
         if( tid < 10000000 ) { // not a plasmid
            vector<TID_T> ptor;
            tax_tree.getPathToRoot(tid,ptor);
            for(unsigned li = 0; li < ptor.size(); ++li) {
               const TID_T ptid = ptor[li];
               
               if( seen.find( ptid ) != seen.end() ) {
                  const string rank = rank_table[ptid];
                  if(VERBOSE) cout<<"seen? "<<ptid<<" "<<rank<<endl;
                  seen_rank.insert(rank);
                  if( li == 1 ) { 
                     last_sibling_ratio = save_last_ratio[ptid];
                  }
               }
            /*
            ** Area for improvement is knowing when to merge calls 
            ** The default is to merge strain calls (with low coverage)
            ** but for the list below, there are numerous near neighbor 
            ** species that likely should be treated as neighbor strains
            */
         
               if( ptid == 262 ) isFrancisella = true;
               if( ptid == 234 ) isBrucella = true;
               if( ptid == 10663 ) isEphage= true;
               if(isEukId(ptid) ) isEuk = true;
               if(isVirId(ptid) ) isVir = true;
               if(isProkId(ptid) ) isProk = true;
               if( ptid == 2 || ptid == 131567 ) continue;
               float gen_copy = 0;
               if( save_cov.find(ptid) == save_cov.end() ) {
                  gen_copy = compKmerCov(kmer_track,ptid,merge_kmer_track,merge_kmer_cnt);
                  save_cov.insert(make_pair(ptid,gen_copy));
               } else {
                  gen_copy = save_cov[ptid];
               }
               if( gen_copy > 0 ) {
                  const int gen_left_over=targ_gen_copy - gen_copy;
                  if(VERBOSE) cout<<"genome copy estimate "<<ptid<<" is "<<gen_copy<<" "<<targ_gen_copy<<" "<<gen_left_over<<endl;
                  signed add_kmer_cnt = merge_kmer_cnt[ptid];
                  // This means potentially more than one genome
                  if( merge_kmer_cnt[ptid] > expect_cov ) {
                     // only add k-mers up to targ_k_kmer_cnt amount
                     add_kmer_cnt=expect_cov;
                  }
                  // how many genomes left?
                  if ( gen_left_over <= 0 ) {
                     assert(add_kmer_cnt >= 0);
                     merge_kmer_cnt[ptid] -= add_kmer_cnt;
                  } else {
                     if(VERBOSE) cout<<sample_iter<<" Warning multiple genome copies may make this complicated"<<" "<<gen_left_over<<" "<<targ_gen_copy <<" "<<gen_copy<<endl;
                  }
                  k_mer_cnt += add_kmer_cnt;
                  if(VERBOSE) cout<<sample_iter<<" Adding to k-mer count "<<ptid<<" "<<add_kmer_cnt<<" remaining: "<<merge_kmer_cnt[ptid]<<" running summ "<<k_mer_cnt<< endl;
               } else {
                  if(VERBOSE) cout<<"Should be No reads for this tax id "<<ptid<<endl;
               }
            }
         } else {
            // this now assumes the only custom ids are plasmids!!!
            isProk=true;
         }
         float pass_thresh = 0;
         if( (isProk && tid >= 10000000) || (separate_plasmid.find(tid) != separate_plasmid.end()) ) {
            //up the min_threshold for plasmid coverage
            if(VERBOSE) cout<<"select plasmid threshold "<<" "<<min_plas<<endl;
            pass_thresh = min_plas;
         } else if( cont_set.find(tid) != cont_set.end() ) {
            if(VERBOSE) cout<<"select contaminant threshold "<<min_contam<<endl;
            pass_thresh = min_contam;
         } else if( isEuk ) {
            if(VERBOSE) cout<<"select euk threshold "<<min_euk<<endl;
            pass_thresh = min_euk;
         } else if( isProk ) {
            if(VERBOSE) cout<<"select prok threshold "<<min_prok<<endl;
            pass_thresh = min_prok;
         } else if( isVir ) {
            if(VERBOSE) cout<<"select vir threshold "<<min_vir<<endl;
            pass_thresh = min_vir;
         }
         if (seen_rank.find("species") != seen_rank.end() || (seen_rank.find("genus") != seen_rank.end() && (isBrucella || isFrancisella || isEphage )) ) {
            if(last_sibling_ratio < 0) {
               if(VERBOSE) cout<<"Warning less than 0 "<<last_sibling_ratio<<endl;
            }
            if( last_sibling_ratio < same_strain_thresh ) {
               if(VERBOSE) cout<<"iter= "<<sample_iter<<" Not enough coverage of the first strain to infer the possibility of additional strains: "<<last_sibling_ratio<<" "<<same_strain_thresh<<endl;
               continue;
            } else {
               if(VERBOSE) cout<<"Wow, we've got enough coverage for the first strain, should we keep inferring additional strains? "<<last_sibling_ratio<<" "<<same_strain_thresh<<endl;
            }
         }
         const float ratio = (float)k_mer_cnt / (float)expect_cov;
         if( ratio >= pass_thresh ) {
            if(VERBOSE) cout<<sample_iter<<" Serious result: "<<tid<<" "<<targ_gen_copy<<" "<<ratio<<" "<<k_mer_cnt<<" "<<expect_cov<<" read_cnt="<<call_vec[i].second<<" thr="<<pass_thresh<<" "<<wrdc<<" wrdc="<<min_wrdcnt<<endl;
            if( expect_cov < 100000 && isEuk ) {
               if(VERBOSE) cout<<"Warning, euk genome has less than 100000 k-mers "<<tid<<" "<<expect_cov<<endl;
               continue;
            }
#if 0
            if( (targ_gen_copy > 1 && isEuk && ratio < 0.1  ) {
               // should be a clear case of mobile element of some sort that probably should not be used
               continue;
            }
#endif
            if( isProk && targ_gen_copy > 1 && ratio < 0.01 ) {
               // starting look like repetitive bacterial elements may need to be ignored.
               if(VERBOSE) cout<<"Warning, possible repetitive element will not include in content summary "<<tid<<" "<<targ_gen_copy<<" "<<ratio<<endl;
               continue;
               
            } else {
               if( ratio >= 1 ) {
                  tot_gen_copy += targ_gen_copy;
               } else {
                  tot_gen_copy += ratio;
               }
               save_out.push_back ( SaveRes(tid,targ_gen_copy,ratio,wrdc) );
               curr_call_set.insert(tid);
               vector<TID_T> ptor;
               tax_tree.getPathToRoot(tid,ptor);
               for(unsigned li = 0; li < ptor.size(); ++li) {
                  const TID_T ptid = ptor[li];
                  if( seen.find( ptid ) == seen.end() ) {
                     if(VERBOSE) cout<<"Set this: "<<ptid<<endl;
                     seen.insert(ptid);
                     save_last_ratio.insert(make_pair(ptid,ratio));
                  }
               }
            }
         } else {
            if(VERBOSE) cout<<"iter="<<sample_iter<<" Fail Serious result: "<<tid<<" "<<targ_gen_copy<<" "<<ratio<<" "<<k_mer_cnt<<" "<<expect_cov<<" read_cnt="<<call_vec[i].second<<" thr="<<pass_thresh<<" "<<wrdc<<" wrdc="<<min_wrdcnt<<endl;
         }
      }

      ofstream ofs(ofbase.c_str());
      typedef pair<TID_T,pair<float,pair<int,int> > > res_t;
      list<SaveRes>::const_iterator ito = save_out.begin(); 
      const list<SaveRes>::const_iterator iso = save_out.end(); 
      ofs<<"Abundance\tGenome Copies\tGenome Covered\tReads\tWReads\tWReads1\tReads1\tTaxID\tName"<<endl;
      for(; ito != iso; ++ito) {
          const int gen_copy = (*ito)._gen_copy;
          float ratio = (*ito)._ratio;
          const TID_T tid = (*ito)._tid;
          const float wrdc = (*ito)._wrdc;
          if( curr_merge_read_cnt.find(tid) == curr_merge_read_cnt.end()) {
              cerr<<"Warning this tid: "<<tid<<" was not found among merged read count, ignore this tid and continue"<<endl;
              continue;
          }
          const int tot_read_cnt = curr_merge_read_cnt[tid];
          float abund = 0;
          if( ratio >= 1 ) {
            abund = gen_copy / tot_gen_copy;
          } else if( ratio >= 0 ) {
            abund = ratio / tot_gen_copy;
          } else {
            assert(tid==HUMAN_TID);
            abund=-1;
            ratio=-1;
          }
          assert(hold_sum.find(tid) != hold_sum.end());
          const string in_str = hold_sum[tid];
          ofs<<abund<<"\t"<<gen_copy<<"\t"<<ratio<<"\t"<<tot_read_cnt<<"\t"<<wrdc<<"\t"<<in_str<<endl; 
      }
      string to = ofbase + ".leftover";
      ofstream left_ofs(to.c_str());
      map<TID_T,int>::const_iterator mit = curr_merge_read_cnt.begin();
      const map<TID_T,int>::const_iterator mis = curr_merge_read_cnt.end();
      for( ; mit != mis; ++mit) {
         const TID_T tid = (*mit).first;
         const int rc = (*mit).second;
         if( call_set.find(tid) == call_set.end() ) {
            left_ofs<<tid<<" "<<rc<<endl;
         }
      }
      call_set = curr_call_set;
      merge_read_cnt = curr_merge_read_cnt;
   }
   cout << "query time: " << clock.stop() << endl;

   return 0; 
}

   static float compReadCnt(const vector<map<TID_T,float> >& track, TID_T tid) {
      float sum = 0;
      for(unsigned thread = 0; thread < track.size(); ++thread) {
         if( track[thread].find(tid) != track[thread].end()) {
            sum += (*track[thread].find(tid)).second;
         }
      }
      return sum;
   }
static float compKmerCov(const vector<map<TID_T,kmer_cnt_t> >& track, TID_T tid, map<TID_T,float>& merge_track, kmer_cnt_t& merge_kmer_cnt) {
   float res = -1.0; 
   if( merge_track.find(tid) == merge_track.end() ) {
      kmer_cnt_t  local_merge;
      uint64_t kmer_cnt = 0;
      for(unsigned thread = 0; thread < track.size(); ++thread) {
         if( track[thread].find(tid) != track[thread].end() ) {
            const kmer_cnt_t& local_track = (*track[thread].find(tid)).second;
            kmer_cnt_t::const_iterator it = local_track.begin();
            const kmer_cnt_t::const_iterator is = local_track.end();
            for(; it != is; ++it ) {
               const kmer_t kmer_val = (*it).first;
               const uint64_t cnt = (*it).second;
               if( local_merge.find(kmer_val) == local_merge.end() ) {
                  local_merge[kmer_val] = cnt;
                  ++kmer_cnt; 
               } else {
                  local_merge[kmer_val] += cnt;
               }
            }
         }
      }
      if(VERBOSE) cout<<"What is the merge count: "<<tid<<" "<<kmer_cnt<<endl;
      merge_kmer_cnt[tid] = kmer_cnt;
      kmer_cnt_t::const_iterator it1 = local_merge.begin();
      const kmer_cnt_t::const_iterator is1 = local_merge.end();
      vector<int> vec(local_merge.size());
      float sum = 0;
      for(unsigned i = 0; it1 != is1; ++it1, ++i ) {
         const int cnt = (*it1).second;
         assert(i < vec.size());
         vec[i] = cnt;
         sum += cnt;
      }
      if( vec.size() > 0 ) {
         sort(vec.begin(),vec.end());
         const int mid = (float)vec.size() / 2.0 ;
         assert(mid >= 0 && mid < (signed)vec.size());
         res = vec[mid];
         sum /= (float)vec.size();
         if(VERBOSE) cout<<"Debug: "<<tid<<" median = "<<res<<" avg = "<<sum<<" "<<vec[0]<<" "<<vec[vec.size()-1]<<endl;
      } else {
         if(VERBOSE) cout<<"Really no k-mers for "<<tid<<endl;
      }
      merge_track.insert(make_pair(tid,res));
   } else {
      res = merge_track[tid];
   } 
   return res;
}
