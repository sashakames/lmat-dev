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
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <gzstream.h>

#define MMAP_SIZE 0

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

using namespace metag;


bool add_root_on_kmer_drop = true;

//static bool verbose=false;
static bool VERBOSE=false;

size_t perm_bytes_allocd;

typedef map<TID_T,list<TID_T> > cand_lin_t;
typedef std::tr1::unordered_map<TID_T,string> rank_map_t;

typedef pair<TID_T,signed> tid_call_t;
typedef pair<TID_T,float> ufpair_t;
typedef map<TID_T,TID_T> hmap_t;
typedef map<TID_T,float> ufmap_t;
typedef map<kmer_t,uint64_t> kmer_cnt_t;
typedef std::tr1::unordered_map<TID_T,TID_T> tid_map_t;

static std::tr1::unordered_set<int> gLowNumPlasmid;

struct TObj {
   bool operator()(const ufpair_t& a, const ufpair_t& b) {
      return (a.second > b.second);
   }
};

#define isPlasmid(tid) ((tid >=10000000 || (gLowNumPlasmid.find(tid) != gLowNumPlasmid.end())) ? true : false)

static void loadLowNumPlasmids(const string& file) {
   ifstream ifs_lst(file.c_str());
   if( !ifs_lst) {
      cerr<<"Unexpected reading error: "<<file<<endl;
      return;
   }
   TID_T pid;
   while(ifs_lst>>pid) {
      gLowNumPlasmid.insert(pid);
   }
}


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

pair<TID_T,float> proc_line(const TaxTree<TID_T> & tax_tree, const string &read_str, const string& seq_buff, int k_size, const cand_lin_t&  call_lin,
                            map< TID_T, kmer_cnt_t >& kmer_track, TID_T call, float call_score, const rank_map_t& rank_table, const tid_map_t& top_strain, const tid_map_t& strain2spec) {

   list<TID_T> cand_lst;
   istringstream istrm(seq_buff.c_str());
   TID_T cand_tid = call;
   float match_score, fnd_score = call_score;
   bool fnd = false;
   vector<ufpair_t> save;
   bool debug=false;
   list<TID_T> plas_lst;
   std::tr1::unordered_set<string> rank_track;
   std::tr1::unordered_set<TID_T> no_dups;
   unsigned rank_cnt = 0;
   // check if tax id is in the lineage of one of the called organisms
   // only check if top score tax id isn't among the final organism call taxids
   if( call_lin.find(cand_tid) != call_lin.end() ) {
         if(debug) cout<<"what the hell: "<<cand_tid<<" "<<match_score<<endl;
         // retrieve original child call and check to see if it was called in this read
         std::tr1::unordered_set<TID_T> once;
         bool fnd_call=false;
         const list<TID_T>& lst_ctid = (*call_lin.find(cand_tid)).second;
         list<TID_T>::const_iterator ib = lst_ctid.begin();
         const list<TID_T>::const_iterator ie = lst_ctid.end();
         for(; ib != ie; ++ib) {
               size_t spos = string::npos;
               TID_T ctid = *ib;
               assert(ctid != cand_tid);
               for(unsigned chk=0; chk < 2; ++chk) {
                  ostringstream os1, os2;
                  //cout<<"look for: "<<ctid<<endl;
                  os1<<" "<<ctid<<" ";
                  const size_t spos1 = seq_buff.find( os1.str().c_str() );
                  os2<<ctid<<" ";
                  const size_t spos2 = seq_buff.find( os2.str().c_str() );
                  /////////////////////////////////
                  // check that the actually org call is somewhere in the list of candidates
                  if( spos1 != string::npos  ) {
                     spos=spos1;
                  }
                  if( spos2 != string::npos && spos2 == 0 ) {
                     spos=spos2;
                  }
                  if( spos != string::npos ) {
                     break;
                  }
                  // for plasmids, don't worry about checking at the species level
                  if( isPlasmid(ctid) ) break;
                  // if a strain taxid is not found in score list, check to see if species is found
                  if( spos == string::npos && strain2spec.find(ctid) !=  strain2spec.end() ) {
                     // note by definition if a strain is present in call list, the species call should not be present in call list
                     ctid = (*strain2spec.find(ctid)).second;
                     // only do this once per strain
                     if( once.find(ctid) != once.end() ) {
                        break;
                     } 
                     once.insert(ctid);
                  }
               }
               if( spos != string::npos ) {
                  const size_t sposa = seq_buff.find( " ", spos+1 );
                  const size_t sposb = seq_buff.find( " ", sposa+1 );
                  float mscore=0;
                  const string debug_str= seq_buff.substr(sposa,sposb-sposa);
                  istringstream is(seq_buff.substr(sposa,sposb-sposa));
                  is>>mscore;
                  fnd=true;
                  fnd_call=true;
                  save.push_back(make_pair(ctid,mscore));
                  if( ctid < 10000000 ) {
                     const string& rval=(*rank_table.find(ctid)).second;
                     if(rank_track.find( rval ) == rank_track.end() ) {
                        rank_track.insert( rval );
                        ++rank_cnt;
                    }
                  }
                  if(debug) cout<<"now check the match call: "<<ctid<<" "<<mscore<<" "<<rank_cnt<<endl;
               }
            }
            // if the calls for this top hit weren't found still register the rank to avoid over aggressive reassignment
            if( !fnd_call ) {
               if( cand_tid < 10000000 ) {
                  const string& rval=(*rank_table.find(cand_tid)).second;
                  if(rank_track.find( rval ) == rank_track.end() ) {
                     rank_track.insert( rval );
                     ++rank_cnt;
                  }
               }

            }
      }
   if( !fnd || save.size() == 0 ) {
         //fnd_score set to 0 for now
      if(debug) cout<<"stick with orig1: "<<call<<endl;

      cand_lst.push_back(call);
   } else if( save.size() == 1 ) {
      // we've got one top best hit to one of our target orgs, so we're done, do a straight reassignment
      if(debug) cout<<"wow one top: "<<save[0].first<<endl;
      fnd_score = save[0].second;
      cand_lst.push_back(save[0].first);
   } else if( save.size() > 1 ) {
      if( rank_cnt == 1 || (rank_cnt == 2 && rank_track.find("species") != rank_track.end() && rank_track.find("strain") != rank_track.end()) ) {
         sort(save.begin(),save.end(),TObj());
         assert(save[0].second >= save[1].second);
         if( save[0].second > save[1].second ) {
            // found a top hit we're done
            if(debug) cout<<"wow one top2: "<<save[0].first<<endl;
            fnd_score = save[0].second;
            cand_lst.push_back(save[0].first);
         } else {
            // if just a tie with a plasmid, you can potentially greedily select a plasmid (what if multiple plasmids?)
            // address that here

            if(debug) cout<<"ties can't resolve revert to "<<call<<endl;
            cand_lst.push_back(call);
         }
      } else {
      // things get complicated with variable ranks, a read may hit a higher order rank just because the genome is more complete
      // so we can't make an assignment in this case without more work so keep as is
         if(debug) cout<<"variable ranks don't try: "<<call<<endl;
         cand_lst.push_back(call);
      }
   }
#if 0
   // IF after all this effort at reassignment, we're still stuck on a species call
   // go ahead and dump it off to whatever the top strain hit is
   if( cand_lst.front() == call && rank_table.find(call) != rank_table.end()) {
      const string& sval = (*rank_table.find(call)).second;
      if( sval == "species" && top_strain.find(call) != top_strain.end()  ) {
         const TID_T strain_id = (*top_strain.find(call)).second;
         cand_lst.clear();
         cand_lst.push_back(strain_id);
      }
   }
#endif
   assert(!cand_lst.empty());
   const TID_T used_tid=cand_lst.front();
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
   string low_num_plasmid_file;
   hmap_t imap;	
   bool skipHuman=false;
   string thresh_str;
   while ((c = getopt(argc, argv, "m:f:ah:n:jb:ye:wp:k:c:v:k:i:d:l:t:sr:o:x:f:g:z:q:")) != -1) {
      switch(c) {
      case 'p':
         low_num_plasmid_file = optarg;
         break;
      case 's':
         skipHuman = true;
         break;
      case 'b':
         cont_genomes_lst = optarg;
         break;
      case 'r':
         rank_table_file = optarg;
         break;
      case 'y':
         VERBOSE = true;
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
#if 0
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
#endif
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
   if( low_num_plasmid_file.length() > 0 ) {
      loadLowNumPlasmids(low_num_plasmid_file);
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
   rank_map_t rank_table;
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
   for(unsigned sample_iter = 0; sample_iter < 1; ++sample_iter) {
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
                  // still want to save any higher level species calls to report back their support
                  hold_sum.insert(make_pair(tid,buff_str));
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
            if( merge_read_cnt.find(tid) == merge_read_cnt.end()) {
               cerr<<"Warning tid not found among merged counts "<<tid<<" will continue "<<endl;
            } else {
               unsigned read_cnt_sum = merge_read_cnt[tid];
               if(VERBOSE) cout<<"step read cnt: "<<tid<<" "<<read_cnt_sum<<endl;
               cand_lst.push_back( make_pair(tid, read_cnt_sum) );
            }
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
      tid_map_t top_strain, strain2spec;
      set<TID_T> leaf_calls;
      map<TID_T,TID_T> get_species;
      cand_lin_t save_call_lin;
      for(unsigned cvi =0; cvi < call_vec.size(); ++cvi) {
         const TID_T tid = call_vec[cvi].first;
         vector<TID_T> localptor;
         tax_tree.getPathToRoot(tid,localptor);
         if(VERBOSE) cout<<"save call lineage: "<<cvi<<" "<<tid<<endl;
         list<TID_T> lst; // intentionally empty
         save_call_lin.insert( make_pair( tid,lst) );
         leaf_calls.insert(tid);
         for(unsigned li = 0; li < localptor.size(); ++li) {
            if( rank_table[localptor[li]] == "species" ) {
               if( top_strain.find( localptor[li] ) == top_strain.end() ) {
                  assert( tid != localptor[li] );
                  top_strain.insert( make_pair(localptor[li],tid) );
               }
               strain2spec.insert( make_pair(tid,localptor[li]) );
               if(VERBOSE) cout<<"save species mapping: "<<tid<<" to "<<localptor[li]<<endl;
            
               get_species.insert( make_pair(tid, localptor[li]) );
            }
            save_call_lin[ localptor[li] ].push_back(tid);
            if(VERBOSE) cout<<"why so slow: "<<tid<<" "<<save_call_lin[localptor[li]].size()<<endl;
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
#pragma omp parallel shared(file_lst, k_size, tax_tree, strain2spec, save_call_lin, top_strain,leaf_calls, kmer_track,weighted_readcnt,read_cnts)  private(ifs,finished, ofname, line, read_count)
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
       pair<TID_T,float> used = proc_line(tax_tree,read_buff, alt_scores, k_size, save_call_lin, kmer_track[thread], taxid, score, rank_table,top_strain, strain2spec);
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
   set<TID_T> seen;
   map<TID_T,float> merge_kmer_track;
   ofstream ofs(ofbase.c_str());
      map<TID_T,list<TID_T> > child;
      typedef pair<TID_T,pair<float,pair<int,int> > > res_t;
      //list<SaveRes>::const_iterator ito = save_out.begin(); 
      //const list<SaveRes>::const_iterator iso = save_out.end(); 
      for(unsigned i =0; i < call_vec.size(); ++i) {
         const TID_T tid = call_vec[i].first;
         vector<TID_T> ptor;
         tax_tree.getPathToRoot(tid,ptor);
         TID_T child_node=tid;
         for(unsigned li = 0; li  < ptor.size(); ++li) {
            const TID_T ptid = ptor[li];
            if( seen.find(child_node) == seen.end() ) {
               seen.insert( child_node );
               if( child.find(ptid) == child.end() ) {
                  list<TID_T> nodes;
                  nodes.push_back(child_node);
                  child.insert( make_pair(ptid,nodes) );
               } else {
                  child[ptid].push_back(child_node);
               }
            }
            child_node= ptid;
         }
      }
      ofs<<"Name\tTaxID\tAbundance\tGenome Copies\tGenome Covered\tReads\tWReads"<<endl;
      map<TID_T,list<char> > tab_lst;
      TID_T cnode = 1; // should always be root
      list<TID_T> open;
      open.push_back(cnode);
      while( !open.empty() ) {
         const TID_T tid = open.front();
         open.pop_front();
         const list<TID_T>& lst = child[tid];
         list<TID_T>::const_iterator ib = lst.begin();
         const list<TID_T>::const_iterator ie = lst.end();
         list<char> chk = tab_lst[tid];
         chk.push_back('\t');
         for(; ib != ie; ++ib) {
            tab_lst[*ib] = chk;
            open.push_front(*ib);
         }
         const unsigned tot_read_cnt = curr_merge_read_cnt[tid];
         float wrdc=0,ratio=0;
         unsigned median_kmer_cnt=0;
         uint64_t expect_cov=0,k_mer_cnt=0;
         if( tot_read_cnt > 0 ) {
            wrdc = compReadCnt(weighted_readcnt,tid);
            median_kmer_cnt = compKmerCov(kmer_track,tid,merge_kmer_track,merge_kmer_cnt);
            expect_cov = ref_kmer_cnt[tid];
            k_mer_cnt = merge_kmer_cnt[tid];
            ratio = (float)k_mer_cnt / (float)expect_cov;
         }
         const string call_name = tax_tree.getName(tid);
         list<char>::const_iterator tab_ib = tab_lst[tid].begin();
         list<char>::const_iterator tab_ie = tab_lst[tid].end();
         for( ; tab_ib != tab_ie; ++tab_ib) {
            ofs<<*tab_ib;
         }
         ofs<<call_name<<"\t"<<tid<<"\t"<<median_kmer_cnt<<"\t"<<ratio<<"\t"<<k_mer_cnt<<"\t"<<expect_cov<<"\t"<<tot_read_cnt<<"\t"<<wrdc<<endl;
      }
      string to = ofbase + ".leftover";
      ofstream left_ofs(to.c_str());
      map<TID_T,int>::const_iterator mit = curr_merge_read_cnt.begin();
      const map<TID_T,int>::const_iterator mis = curr_merge_read_cnt.end();
      for( ; mit != mis; ++mit) {
         const TID_T tid = (*mit).first;
         const int rc = (*mit).second;
         if( curr_call_set.find(tid) == curr_call_set.end() ) {
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
