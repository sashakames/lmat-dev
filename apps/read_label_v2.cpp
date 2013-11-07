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
#include "kencode.hpp"
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"
#include <gzstream.h>

#define MMAP_SIZE 0
#define TID_T uint32_t
#define DBTID_T uint32_t


using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

using namespace kencode_ns;
using namespace metag;

static bool verbose=false;

bool add_root_on_kmer_drop = true;

size_t perm_bytes_allocd;

typedef pair<TID_T,float> ufpair_t;
typedef pair<TID_T,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<int16_t,tax_data_t> label_info_t;
typedef map<TID_T,TID_T> hmap_t;
typedef map<TID_T,float> ufmap_t;
typedef map<uint16_t, ufmap_t> u_ufmap_t;

vector <int> read_len_vec;
vector <int> read_len_avgs;

#define _USE_KPATH_IDS 0


my_map tid_rank_map;
id_convback_map_t conv_map;


// this determines which style of tid reduction we want, default is
// the "comphrehensive reduction; requires the tid to rank mapping
// table.  -w option enables strain to species only reduction.
bool tid_map_is_strain_species = false;

bool badGenomes(TID_T tid) {
   bool isBad=false;
#if _USE_KPATH_IDS == 1
   switch(tid) {
      case 1154764:
      case 1218173:
         isBad=true;
         break;
      }
#else 
   switch(tid) {
      // comment in NCBI is that these genomes are likely HIV-1 but labeled only HIV
      // Thus their distinct lineage on another branch away from HIV-1 confounds HIV-1 labeling!
      // for now ignore sequences with this taxid.
      case 12721:
      case 693660:
         isBad=true;
         break;
   }
   
#endif
      return isBad;
}

// vec is assumed to be sorted
int closest(int value)
{
  unsigned i;

  for (i =0; i<read_len_avgs.size(); i++) {

    if (value <= read_len_avgs[i]) 
      return read_len_vec[i];

  }

  return (read_len_vec[i]);
  
}



static int getReadLen(int rl)
{
  int len = closest(rl);

  if (len > 0)
    return len;

  return 80;

}


typedef std::pair<TaxNode<TID_T>*,TID_T> tncpair_t;

static bool isAncestor(const TaxTree<TID_T>& tax_tree, TID_T prev_taxid /*ancestor */, TID_T curr_taxid /* descendant */ ) {
   bool isOk=false;
   vector<TID_T> path;
   TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
   tax_tree_tmp.getPathToRoot(curr_taxid,path);
   for(unsigned p = 0; p < path.size(); ++p) {
      if(path[p] == prev_taxid) {
         isOk=true;
         break;
      }
   }
   return isOk;
}


struct SimpleCmp {
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
      return a.second > b.second;
   }
};

struct CmpDepth {
   CmpDepth(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
      const int adepth = (*_imap.find(a.first)).second;
      const int bdepth = (*_imap.find(b.first)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

struct CmpDepth1 {
   CmpDepth1(const hmap_t& imap) : _imap(imap) {}
   bool operator()(TID_T a, TID_T b) const {
      const int adepth = (*_imap.find(a)).second;
      const int bdepth = (*_imap.find(b)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};



enum nomatch_t { eReadTooShort=0, eNoDbHits };
enum match_t { eDirectMatch, eMultiMatch, ePartialMultiMatch, eNoMatch, eNoLCAError };
static string match2str(match_t match) {
   string str = "error";
   switch(match) {
   case eDirectMatch : 
      str="DirectMatch";
      break;
   case eMultiMatch : 
      str="MultiMatch";
      break;
   case ePartialMultiMatch : 
      str="PartialMultiMatch";
      break;
   case eNoMatch : 
      str="NoMatch";
      break;
   case eNoLCAError : 
      str="LCA_ERROR";
      break;
   }
   return str;
}

static bool addToCandLineage(const ufpair_t cand, list<ufpair_t>& lineage, const hmap_t& dmap, const TaxTree<TID_T>& tax_tree) {
   bool addNode = false;
   if( lineage.empty() ) {
      addNode=true;
   } else {
      const hmap_t::const_iterator mtch0 = dmap.find(cand.first);
      unsigned cand_depth = 0;
      if( mtch0 != dmap.end()) {
         cand_depth = (*mtch0).second;
      }
      addNode = true;
      list<ufpair_t>::iterator it = lineage.begin();
      const list<ufpair_t>::iterator is = lineage.end();
      for(; it != is; ++it) {
         const TID_T taxid = (*it).first;
         const hmap_t::const_iterator mtch = dmap.find(taxid);
         unsigned chk_depth = 0;
         if( mtch != dmap.end() ) {
            chk_depth = (*mtch).second;
         }
         if( chk_depth > cand_depth && !isAncestor(tax_tree,cand.first,taxid) ) {
            addNode=false;
            break;
         } else if( chk_depth < cand_depth && !isAncestor(tax_tree,taxid,cand.first)) {
            addNode=false;
            break;
         } else if( chk_depth == cand_depth ) {
            addNode=false;
            break;
         }
      }
   }
   if( addNode ) {
      if(verbose) cout<<"Add to Lineage: "<<cand.first<<" "<<cand.second<<endl;
      lineage.push_back(cand);
   }
   return addNode; 
}

static bool isHuman(TID_T taxid)  {
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
   //case 44004: //neanderthal
#endif
      res = true;
      break;
   default:
      res=false;
      break;
   }
   return res;
}

static bool cmpCompLineage(ufpair_t cand, const vector<ufpair_t>& lineage, set<TID_T>& no_good, float diff_thresh, const TaxTree<TID_T>& tax_tree) {
   const float undef = -10000;
   bool keep_going=true;
   for(unsigned i = 0; i < lineage.size(); ++i) {
      if( isAncestor(tax_tree, lineage[i].first, cand.first )) {
         break;
      } 
      if( lineage[i].second != undef &&  (lineage[i].second - cand.second) > diff_thresh  ) {
         if(verbose) cout<<"competing lineage is too far to care: "<<cand.first<<" "<<cand.second<<" "<<lineage[i].first<<" "<<lineage[i].second<<" "<<diff_thresh<<endl;
         keep_going=false;
         break;
      }
      if( (lineage[i].second - cand.second) <= diff_thresh ) {
         if(verbose) cout<<"competing lineage is too close: "<<cand.first<<" "<<cand.second<<" "<<lineage[i].first<<" "<<lineage[i].second<<" "<<diff_thresh<<endl;
         no_good.insert(lineage[i].first);
      }
   }
   return keep_going;
}

static pair<ufpair_t,match_t> 
findReadLabelVer2(const vector<ufpair_t>& rank_label, float diff_thresh, const TaxTree<TID_T>& tax_tree, const hmap_t& taxid2idx, list<ufpair_t>& cand_lin, const hmap_t& dmap, const ufmap_t& all_cand_set) {
   match_t match = eNoMatch;
   unsigned lowest_depth = 0, highest_depth = 0;
   ufpair_t lowest=make_pair(0,0), highest = make_pair(0,0);
   signed lidx = -1;
   assert( cand_lin.empty() );
   for(signed i = rank_label.size()-1; i >= 0; --i) {
      if( !addToCandLineage(rank_label[i],cand_lin,dmap, tax_tree) )  {
         lidx = i;
         break;
      } else {
         const hmap_t::const_iterator mtch = dmap.find(rank_label[i].first);
         if( (*mtch).second> lowest_depth || (i == (signed)rank_label.size()-1) ) {
            lowest = rank_label[i];
            lowest_depth = (*mtch).second;
         }
         if( (*mtch).second < highest_depth || i == (signed)rank_label.size()-1 ) {
            highest = rank_label[i];
            highest_depth = (*mtch).second;
         }
      }
   }
   set<TID_T> add_set;
   if( highest_depth != 0 ) {
      vector<TID_T> path;
      TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
      tax_tree_tmp.getPathToRoot(highest.first,path);
      for(unsigned i = 0; i < path.size(); ++i) {
         add_set.insert( path[i] );
         ufmap_t::const_iterator mtch = all_cand_set.find( path[i] );
         if( mtch != all_cand_set.end() ) {
            ufpair_t val = make_pair( path[i], (*mtch).second);
            cand_lin.push_back(val);
         } else {
            //compute log_odds
            const float undef = -10000;
            cand_lin.push_back( make_pair(path[i], undef) );
         }
      }
   } 
   vector<ufpair_t> cand_lin_vec(cand_lin.size());
   list<ufpair_t>::const_iterator it = cand_lin.begin();
   const list<ufpair_t>::const_iterator is = cand_lin.end();
   for(unsigned i= 0; i < cand_lin_vec.size(); ++it, ++i) {
      cand_lin_vec[i] = *it;
   }
   CmpDepth cd(dmap);
   sort(cand_lin_vec.begin(),cand_lin_vec.end(),cd);

   set<TID_T> no_good;

   for(signed i = lidx; i >= 0; --i) {
      if( add_set.find( rank_label[i].first ) == add_set.end() ) {
         if( !cmpCompLineage(rank_label[i],cand_lin_vec,no_good,diff_thresh,tax_tree) )  {
            if(verbose) cout<<" competing taxid too low scoring, quit search "<<rank_label[i].first<<" "<<rank_label[i].second<<endl;
            break;
         }
      }
   }
   ufpair_t taxid_call;
   if( cand_lin.empty() && no_good.empty() ) {
      match = eNoMatch;
   } else if( !cand_lin.empty() && no_good.empty() ) {
      taxid_call = lowest; 
      match = eDirectMatch;
   } else {
      vector<ufpair_t> cand_vec(cand_lin.size());
      list<ufpair_t>::const_iterator it = cand_lin.begin(); 
      const list<ufpair_t>::const_iterator is = cand_lin.end(); 
      unsigned cnt= 0;
      for(cnt= 0; it != is; ++it, ++cnt) {
         if(verbose) cout<<"merging cand_lst "<<cnt<<" "<<(*it).first<<" "<<(*it).second<<endl;
         cand_vec[cnt] = *it;
      }  
      CmpDepth cd(dmap);
      sort(cand_vec.begin(),cand_vec.end(),cd);
      float sval = 0, max_val = -10000;
      pair<TID_T,bool> res = make_pair(0,false);
      int root_idx = -1;
      for(unsigned i = 0; i < cand_vec.size(); ++i) {
         const TID_T tax_i = cand_vec[i].first;
         max_val = std::max(cand_vec[i].second,max_val);
         sval += cand_vec[i].second;
         ++cnt; 
         if( no_good.find(tax_i) == no_good.end() ) {
            res = make_pair( cand_vec[i].first, true );
            root_idx = i;
            break;
         }
      }
      if(!res.second) {
         taxid_call = make_pair(0,-1);
         match = eNoLCAError; 
      }  else {
         const TID_T lca_tid = res.first;
         match = eMultiMatch; 
         if( all_cand_set.find(lca_tid) != all_cand_set.end() ) {
            assert(root_idx != -1);
            if( max_val < cand_vec[root_idx].second ) {
               match = ePartialMultiMatch;
               max_val = cand_vec[root_idx].second;
            }
         }
         taxid_call = make_pair(lca_tid, max_val);
      }
   }
   if(verbose) cout<<"I'm confused: "<<taxid_call.first<<" "<<taxid_call.second<<endl;
   return make_pair(taxid_call,match);
}

void
fill_in_labels(const TaxTree<TID_T>& taxtree, vector<TID_T>& row, const tax_data_t& taxids, const hmap_t& idx2taxid, const hmap_t& taxid2idx) {
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const TID_T col_tax_val = (*mtch).second;
      const TaxTree<TID_T>::const_iterator col_tax_val_it = taxtree.find(col_tax_val); 
      if(col_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
         continue;
      }
      TaxNode<TID_T>* col_tax_val_node = (*col_tax_val_it).second;
      // see if this id's genome count needs to be update
      //  based on its children   - can skip if its a leaf
      if( !col_tax_val_node->isLeaf() && row[tax_idx] == 0 ) {
         tax_data_t::const_iterator it = taxids.begin();
         const tax_data_t::const_iterator is = taxids.end();
         unsigned has_cnt_sum = 0;
         for(; it != is; ++it) {   
            const TID_T tax_id = (*it).first;
            const TaxTree<TID_T>::const_iterator tax_val_it = taxtree.find(tax_id); 
            if(tax_val_it == taxtree.end()) {
               cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
               continue;
            }
            TaxNode<TID_T>* tax_val_node = (*tax_val_it).second;
            // test if tax_id is a descendant of tax_val
            // save max value
            const hmap_t::const_iterator idx_mtch = taxid2idx.find(tax_id);
            assert(idx_mtch != taxid2idx.end());
            const TID_T tax_val_idx = (*idx_mtch).second;
            
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
      const TID_T tax_val = (*mtch).second;
      if(verbose) cout<<"["<<tax_val<<","<<row[tax_idx]<<"]";
   }
   if(verbose) cout<<endl;
}

struct TCmp { 
   TCmp(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<TID_T,float>& a, const pair<TID_T,float>& b) const {
	if( fabs(a.second- b.second) < 0.001 ) {
		const int adepth = (*_imap.find(a.first)).second;
		const int bdepth = (*_imap.find(b.first)).second;
		return adepth < bdepth;
        } 
	return a.second < b.second; }
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

static void loadRandHits(const string& file_lst, u_ufmap_t& rand_hits_all) {
   ifstream ifs_lst(file_lst.c_str());
   if( !ifs_lst) {
      cerr<<"Unexpected reading error: "<<file_lst<<endl;
      return;
   }
   int read_len;
   string file;
   while(ifs_lst>>read_len>>file) {
      const char* path = getenv("LMAT_DIR");
      if(path ) {
         file = string(path) + string("/") + file;
      } 
      cout<<"load: "<<read_len<<" "<<file<<endl;
      read_len_vec.push_back(read_len);


      ifstream ifs_pre(file.c_str());
      if( !ifs_pre) {
         cerr<<"Unexpected reading error: "<<file<<endl;
         continue;
      }
      ifs_pre.close();
      igzstream ifs(file.c_str());

      const int buff_size=20004;
      char buff[buff_size];
      while(ifs.getline(buff,buff_size)) {
         istringstream istrm(buff);
         int32_t taxid;
         float min_val=0,max_val=0,mean_val=0,stdev=0;
         istrm>>taxid>>min_val>>max_val>>mean_val>>stdev;
         const float cutoff = max_val;
         //cout<<"debug: crap"<<read_len<<" "<<taxid<<" "<<cutoff<<endl;
         rand_hits_all[read_len][taxid]=cutoff;
      }
   } 
   sort(read_len_vec.begin(),read_len_vec.end());
   int i; 

   for (i = 1; i < (signed)read_len_vec.size(); i++) {

     read_len_avgs.push_back((read_len_vec[i-1] + read_len_vec[i]) / 2);

   }
     

}

static float log_odds_score(float label_prob, float random_prob) {
      //if( label_prob >= 1 ) label_prob = 0.9999;
      //if( random_prob >= 1 ) random_prob = 0.9999;
      //const float denom = (1-label_prob)*random_prob; // must be > zero
      //const float numer =  label_prob*(1-random_prob); // must be > zero
      //if(verbose)  cout<<"log_odds_calc: "<<numer<<" "<<denom<<endl;
      const float numer = label_prob;
      const float denom = random_prob <= 0 ? 0.00001 : random_prob;
      const float log_odds = log( numer / denom );
      return log_odds;
}

pair<ufpair_t,match_t>
construct_labels(const TaxTree<TID_T>& tax_tree, const vector<label_info_t>& label_vec, const list<TID_T>& taxid_lst, const hmap_t& tax2idx, const hmap_t& idx2taxid, ofstream& ofs, size_t mer_len, const ScoreOptions& sopt) {
   const unsigned num_tax_ids = taxid_lst.size();
   vector<bool> any_kmer_match(label_vec.size(),false) ;
   unsigned cnt_fnd_kmers=0;
   vector<vector<TID_T> > label_matrix(label_vec.size());
   uint16_t cand_kmer_cnt = 0;
   uint16_t debug_bad_cand_kmer_cnt = 0;
   for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
      if(label_vec[pos].first >= 0 ) ++cand_kmer_cnt;
      else {
         if(verbose) cout<<"Debug bad kmers pos: "<<pos<<" "<<label_vec[pos].first<<endl;
         ++debug_bad_cand_kmer_cnt;
      }
      label_matrix[pos].resize(num_tax_ids,0);
      tax_data_t::const_iterator it = label_vec[pos].second.begin();
      const tax_data_t::const_iterator is = label_vec[pos].second.end();
      unsigned any=0;
      for(; it != is; ++it) {
         const TID_T tax_id = (*it).first;
         const uint16_t has_cnt = (*it).second;
         if( tax2idx.find(tax_id) == tax2idx.end() ) {
            cout<<"HOUSTON WE HAVE A PROBLEM: "<<tax_id<<endl;
         }
         const unsigned idx = (*(tax2idx.find(tax_id))).second;
         label_matrix[pos][idx] = has_cnt;
         ++any;
         if(verbose) cout<<"check: "<<pos<<" "<<idx<<" "<<tax_id<<" "<<has_cnt<<endl;
      }   
      any_kmer_match[pos] = !label_vec[pos].second.empty();
      if( any_kmer_match[pos] ) {
         ++cnt_fnd_kmers;
      }
   }
      vector<ufpair_t> rank_label(num_tax_ids,make_pair(0,0)); 
      unsigned sig_hits = 0;
      float log_sum=0.0;
      ufmap_t all_cand_set;
      bool hasHuman=false;
      for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         float found_genome_cnt = 0;
         const TID_T taxid = (*(idx2taxid.find(tax_idx))).second;
         if( isHuman(taxid) ) hasHuman=true;
         if(verbose) cout<<"col: "<<tax_idx<<" "<<taxid;
         bool noMatch=false;
         for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
            if( label_matrix[pos][tax_idx] > 0 ) {
               found_genome_cnt += 1;
               if(verbose) cout<<" pos="<<pos<<" "<<label_vec[pos].first<<" "<<label_matrix[pos][tax_idx]<<" taxid="<<taxid<<" "<<found_genome_cnt<<endl;
            }
         }
         float label_prob = 0.0;
         if( !noMatch ) {
            label_prob = (float)found_genome_cnt / (float) cand_kmer_cnt;
         }
         if(verbose) cout<<" "<<label_prob<<endl;
         float log_odds=label_prob;
         float random_prob = 0.0001;
         if( !sopt._comp_rand_hits ) {
            const  int cand_kmer_cnt_match = getReadLen(cand_kmer_cnt);
            u_ufmap_t::const_iterator mtch = sopt._rand_hits.find(cand_kmer_cnt_match);
            if( mtch != sopt._rand_hits.end() ) {
               const ufmap_t& rand_hits = (*mtch).second;
               const TID_T stid = taxid;
               if( rand_hits.find(stid) != rand_hits.end()) {
                  const float val = (*(rand_hits.find(stid))).second;
                  random_prob = val+0.0001;
               } else {
                  cerr<<"ERROR, ALL TAXIDS MUST HAVE NULL MODELS: "<<taxid<<" "<<cand_kmer_cnt_match<<" "<<cand_kmer_cnt<<endl;
                  random_prob = 1.0; // basically try to ignore these tax ids for now.
               } 
            }
            log_odds = log_odds_score(label_prob,random_prob);
            if(verbose) cout<<" cand kmer_cnt "<<cand_kmer_cnt_match<<" in="<<cand_kmer_cnt<<" "<<debug_bad_cand_kmer_cnt<<" "<<label_matrix.size()<<endl;
         }           
         rank_label[tax_idx] = make_pair(taxid,log_odds);
         all_cand_set.insert( rank_label[tax_idx] );
         log_sum += log_odds;
         sig_hits++;
         if(verbose) cout<<"Sig match chec: "<<taxid<<" "<<rank_label[tax_idx].second<<" "<<label_prob<<" "<<random_prob<<" "<<endl;
      }
      list<ufpair_t> valid_cand;
      pair<ufpair_t,match_t> res = make_pair(make_pair(0,0),eNoMatch);   
      string matchType=match2str(res.second);
      const float log_avg = sig_hits > 0 ? log_sum / (float)sig_hits : 0;
      float log_std=0, max_score = 0;
      bool first=true;
      for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         if( rank_label[tax_idx].second > 0 ) { 
            if( rank_label[tax_idx].second > max_score || first ) {
               max_score = rank_label[tax_idx].second;
               first=false;
            }
         }
         const float val = log_avg - rank_label[tax_idx].second;
         log_std += (val*val);
      }
      float stdev1 = sig_hits > 1 ? sqrt(log_std/(sig_hits-1)) : 0;
      if( sig_hits > 0 ) {
         if( hasHuman ) {
            for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
               //float found_genome_cnt = 0;
               const TID_T taxid = (*(idx2taxid.find(tax_idx))).second;
               if( isHuman(taxid) ) {
                  rank_label[tax_idx].second += (sopt._diff_thresh2*stdev1);
               }
            } 
         }
         TCmp tcmp(sopt._imap);
         sort(rank_label.begin(),rank_label.end(),tcmp);
         ofs<<log_avg<<" "<<stdev1<<" "<<cand_kmer_cnt<<"\t";
         stdev1 *= sopt._diff_thresh;
         
         res = findReadLabelVer2(rank_label,stdev1,tax_tree,tax2idx,valid_cand,sopt._imap,all_cand_set);    
         if(sopt._prn_all) {
            bool prn=false;
            for(signed i = rank_label.size()-1; i >= 0; --i) {
               if( rank_label[i].second >= 0  || verbose ) {
                  ofs<<" "<<rank_label[i].first<<" "<<rank_label[i].second;
                  prn=true;
               }
            }
            if(!prn) {
               ofs<<"-1 -1";
            }
            ofs<<"\t";
         }
         matchType=match2str(res.second);
      }
    
      ufpair_t best_guess = make_pair(0,0);
      if( res.second == eDirectMatch) {
         best_guess = res.first;
         ofs<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }  else if( res.second == eMultiMatch || res.second == ePartialMultiMatch ) {
         if( !sopt._prn_all) {
            list<ufpair_t>::const_iterator it = valid_cand.begin();
            const list<ufpair_t>::const_iterator is = valid_cand.end();
            for( ; it != is; ++it) {
               ofs<<" "<<(*it).first<<" "<<(*it).second;
            }
            if(valid_cand.empty() ) {
               ofs<<"-1 -1";
            }
            ofs<<"\t";
         }
         best_guess = res.first;
         ofs<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }
      else if( res.second == eNoMatch ) {
         ofs<<-1<<" "<<-1<<" "<<matchType;
      } else {
         cerr<<"Unexpected match type"<<endl;
      }
      ofs<<endl;
      //return make_pair(best_guess.first,make_pair(res.second,res.first));
      return make_pair(best_guess,res.second);
}

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

int retrieve_kmer_labels(INDEXDB<DBTID_T>* table, const char* str, const int slen, const kmer_t klen, 
                          vector<label_info_t>& label_vec, list<TID_T>& taxid_lst, hmap_t& tax2idx, hmap_t& idx2tax, 
                          const hmap_t& dmap, const TaxTree<TID_T>& tax_tree, uint16_t max_count) 
   {
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
           const int pos = j-klen+1;
           /* kmer_lookup(kmer); do k-mer lookup here... */
           label_vec[pos].first = 0; // marks the position as having a valid k-mer (for case where n-masked reads are used)
           if(verbose) cout<<"debug valid: "<<j<<" "<<pos<<endl;
           if( no_dups.find(kmer_id) != no_dups.end() ) continue;
           no_dups.insert(kmer_id);

           TaxNodeStat<DBTID_T> *h = new TaxNodeStat<DBTID_T>(*table);
           if(verbose) cout<<"lookup kmer at posit:"<<j<<endl;

#if (DBTID_T == uint32_t)
           h->begin(kmer_id, tid_rank_map,  max_count, tid_map_is_strain_species);
#else
           h->begin(kmer_id, &conv_map);
#endif
           unsigned dcnt = 0, mtch = 0;
           list<TID_T> obs_tids;
           while( h->next() ) {
              TID_T tid = h->taxid();
              if( tid == 20999999 || badGenomes(tid)) continue;
              uint16_t ng = h->taxidCount();
              //start: hack to drop kmers with too many tids
	      /* 
             if (ng > max_count && max_count != 0) {
                 if(verbose) cout<<"Too many tids "<<ng<<" max="<<max_count<<endl;
                 tid = 1;
                 label_vec[pos].first = 1;
                 const uint16_t pr_cnt = 1;
                 obs_tids.push_back(tid);
                 label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
                 if( tax2idx.find(tid) == tax2idx.end() ) {
                    const unsigned idx = taxid_lst.size();
                    tax2idx[tid] = idx;
                    idx2tax[idx] = tid;
                    taxid_lst.push_back(tid);
                 }
                 break;
              } 
	      */
 //end: hack to drop kmers with too many tids
               if(dcnt==0) {
                  if( ng <= 0 ) {
                     cout<<"Warning unexpected value for ng: "<<ng<<endl;
                     ng=1;
                  }
                  label_vec[pos].first = ng;
               }
               //collect the unique set of tax ids
               //const uint16_t pr_cnt = h->present();
               const uint16_t pr_cnt = 1;
               obs_tids.push_back(tid);
               label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
               if( tax2idx.find(tid) == tax2idx.end() ) {
                  const unsigned idx = taxid_lst.size();
                  tax2idx[tid] = idx;
                  idx2tax[idx] = tid;
                  taxid_lst.push_back(tid);
               }
               // Note: ng should not change with multiple calls to next here
               if(verbose) {
                  if( dcnt == 0 ) cout<<"debug kmer "<<kmer_id<<" gc="<<ng;
                  cout<<" ["<<pos<<" "<<tid<<" "<<pr_cnt<<"]" ;
               }
               dcnt++;
               ++mtch;
           }
           vector<TID_T> obs_tids_vec(obs_tids.size());
           list<TID_T>::const_iterator it = obs_tids.begin();
           const list<TID_T>::const_iterator is = obs_tids.end();
           for(unsigned i =0; it != is; ++it, ++i) {
              obs_tids_vec[i] = *it;
           }
           CmpDepth1 cd(dmap);
           sort(obs_tids_vec.begin(),obs_tids_vec.end(),cd);
           int last_depth = -1;
           for(unsigned i = 0; i < obs_tids_vec.size(); ++i) {
             const TID_T tid=obs_tids_vec[i];
             const int depth = (*dmap.find(tid)).second;
             if( depth == 0 ) {
               break;
             }
             if( last_depth == depth || last_depth == -1 ) {
                  vector<TID_T> path;
                  TaxTree<TID_T>& tax_tree_tmp = const_cast<TaxTree<TID_T>&>(tax_tree);
                  tax_tree_tmp.getPathToRoot(tid,path);
                  for(unsigned p = 0; p < path.size(); ++p) {
                     const TID_T ptid = path[p];
                     //if( verbose ) cout<<"lineage tids added: "<<ptid<<" for "<<tid<<" pos="<<pos<<" "<<p<<endl;
                     label_vec[pos].second.insert( make_pair(ptid,1) );
                     if( tax2idx.find(ptid) == tax2idx.end() ) {
                        const unsigned idx = taxid_lst.size();
                        tax2idx[ptid] = idx;
                        idx2tax[idx] = ptid;
                        if( verbose ) cout<<"lineage tids added: "<<ptid<<" for "<<tid<<" pos="<<pos<<endl;
                        taxid_lst.push_back(ptid);
                     }
                  }
             } else {
               break;
             }
          }
          if( verbose ) cout<<"num taxids: "<<taxid_lst.size();
          if(dcnt > 0 && verbose ) cout<<" end k-mer lookup"<<endl;
          if(mtch == 0 && verbose ) cout<<" no k-mer matches "<<endl;
          delete h;
       }
   }
   return valid_kmers;
}

void proc_line(const TaxTree<TID_T>& tax_tree, int ri_len, string &line, int k_size, INDEXDB<DBTID_T> *table, ofstream &ofs, float threshold, const ScoreOptions& sopt, uint16_t max_count
               ,map<TID_T,int>& track_taxids, map<nomatch_t,int>& track_nomatch, map<TID_T,float>& track_tscores, float min_label_score, int min_kmer) {

     if(ri_len < 0 || ri_len > (signed)line.length()) {
         cout<<"unexpected ri_len value: "<<ri_len<<endl;
         return;
     } else if( ri_len < k_size ) {
         ofs<<"-1 -1 -1"<<"\t-1 -1\t"<<ri_len<<" "<<k_size<<" ReadTooShort"<<endl;
        if( track_nomatch.find(eReadTooShort) != track_nomatch.end()) {
            track_nomatch[eReadTooShort] += 1;
        } else {
            track_nomatch[eReadTooShort] = 1;
        }
     } else {
        vector<label_info_t> label_vec(ri_len-k_size+1,make_pair(-1,tax_data_t()));
        list<TID_T> taxid_lst; 
        hmap_t tax2idx, idx2tax;
        const int valid_kmers = retrieve_kmer_labels(table, line.c_str(), ri_len, k_size,label_vec,taxid_lst,tax2idx,idx2tax, sopt._imap, tax_tree, max_count);
        if( !taxid_lst.empty() ) {  
           pair< ufpair_t, match_t> mtch = construct_labels(tax_tree,label_vec,taxid_lst,tax2idx,idx2tax,ofs,k_size,sopt); 
           if(mtch.second == eNoMatch ) {
              if( track_nomatch.find(eNoDbHits) != track_nomatch.end()) {
                  track_nomatch[eNoDbHits] += 1;
              } else {
                  track_nomatch[eNoDbHits] = 1;
              }
           } else if( mtch.first.second >= min_label_score && valid_kmers >= min_kmer) {
               if( track_taxids.find(mtch.first.first) == track_taxids.end()) {
                  track_taxids[mtch.first.first] = 1;
                  track_tscores[mtch.first.first] = mtch.first.second;
               } else {
                  track_taxids[mtch.first.first] += 1;
                  track_tscores[mtch.first.first] += mtch.first.second;
               }
           }
           
        } else {
           ofs<<"-1 -1 "<<valid_kmers<<"\t-1 -1\t"<<ri_len<<" "<<k_size<<" NoDbHits"<<endl;
           if( track_nomatch.find(eNoDbHits) != track_nomatch.end()) {
               track_nomatch[eNoDbHits] += 1;
           } else {
               track_nomatch[eNoDbHits] = 1;
           } 
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

void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:restore persistent db] [-h:max_tid_count\n";
  cout << "[-r <rank/tid-map-file>] [-h <tid-cutoff>] [-w <with-strain-species-map> (affects -r option)]\n";
  cout << "note: -r makes -k unnecessary\n"; 
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1;
   int n_threads = 0;

   bool restore = true;

   float threshold = 0.0, min_score = 0.0;
   int min_kmer = 35;

   string rank_ids, kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, gid_to_tid_fn, depth_file, rand_hits_file, rank_table_file, id_bit_conv_fn;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   size_t mmap_size = 0;
   bool fastq=false; 
   uint16_t max_count = ~0;
   bool prn_read = true;

   while ((c = getopt(argc, argv, "u:ah:n:j:b:ye:wmpk:c:v:k:i:d:l:t:s:r:o:x:f:g:z:q:")) != -1) {
      switch(c) {
      case 'f':
	id_bit_conv_fn = optarg;
	break;
      case 'j':
         min_kmer = atoi(optarg);
         break;
      case 'u':
         rank_ids = optarg;
         break;
      case 'x':
         min_score = atof(optarg);
         break;
      case 'a':
         prn_read=false;
         break;
      case 'w':
	      tid_map_is_strain_species = true;
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
      case 'n':
         rand_hits_file=optarg;
         break;
#if 0
      case 'j':
         verbose = true;
         break;
#endif
      case 'b':
         sopt._diff_thresh = atof(optarg);
         break;
      case 'l':
         sopt._diff_thresh2 = atof(optarg);
         break;
      case 'y':
         verbose = true;
         break;
      case 'e':
         depth_file = optarg;
         break;
      case 'q':
         fastq=true;
         break;
      case 'p':
         sopt._prn_all = true;
         break;   
      case 'r':
	      rank_table_file = optarg;
        break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
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
      case 'g':
         gid_to_tid_fn = optarg;
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

   if (depth_file == "") cout << "depth_file\n";
   if (ofbase == "") cout << "ofbase\n";
   if (n_threads == 0) cout << "n_threads\n";
   if (kmer_db_fn == "") cout << "kmer_db_fn\n";
   if (query_fn == "") cout << "query_fn\n";

   if (depth_file == "" ||  ofbase == "" || n_threads == 0 || kmer_db_fn == "" || query_fn == "")  {
     cout<<ofbase<<" "<<n_threads<<" "<<kmer_db_fn<<" "<<query_fn<<" "<<depth_file<<endl; 
     usage(argv[0]);
     return -1;

   }
   if (!restore && k_size == -1)  {
     cout << "missing kmer size!\n";
     usage(argv[0]);
     return -1;

   }
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


   cout << "Start kmer DB load\n";
   INDEXDB<DBTID_T> *taxtable;


#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDB<DBTID_T>*, std::size_t> ret = mfile.find<INDEXDB<DBTID_T>>("KmerDB");

     taxtable = ret.first;
     k_size = taxtable->get_kmer_length();
     cout << "k size:  " << k_size  <<   endl ;
     taxtable->conv_ptrs();
#else


#if WITH_PJMALLOC == 1

   if (restore) {
      
     perm(&taxtable, sizeof(taxtable));
     mopen(kmer_db_fn.c_str(), "r", mmap_size);
     //mopen(kmer_db_fn.c_str(), "r+", mmap_size);

     if (k_size < 1)
       k_size = taxtable->get_kmer_length();

     cout << "num kmers: " << taxtable->size() << " - " << k_size  <<   endl ;


   } else 
    
#endif

{


#if (USE_SORTED_DB == 0)
     taxtable = new INDEXDB<DBTID_T>;

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
     taxtable->ingest();
#endif


   }

#endif
   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;

   sopt._comp_rand_hits = (rand_hits_file.length() == 0);
   if( !sopt._comp_rand_hits ) {
      loadRandHits(rand_hits_file,sopt._rand_hits);
   }
   if( k_size <= 0 ) {
      cerr<<"Unable to read database, k_size="<<k_size<<endl;
      return -1;
   }
   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   size_t *arr = split_file(n_threads, ifs);
   ifs.close();

   bool finished;

   size_t pos = 0 ;

   ofstream ofs;

   if (rank_table_file.length() > 0) {
     if(max_count == ~0) {
       cout << "Need to set -h <tid-cutoff> to use rank file map!\n";
     } else {
       FILE * rmfp = fopen(rank_table_file.c_str(), "r");
     
       uint32_t src, dest;
       
       while (fscanf(rmfp,"%d%d", &src, &dest) > 0) {
	      tid_rank_map[src] = dest;
       }
       
       fclose(rmfp);
       
     }
	 
   }



   cout<<"Read taxonomy tree: "<<tax_tree_fn<<endl;
   TaxTree<TID_T> tax_tree(tax_tree_fn.c_str());

   cout<<"Read taxonomy depth: "<<depth_file<<endl;
   ifstream ifs1(depth_file.c_str());
   if(!ifs1) {
	cerr<<"unable to open: "<<depth_file<<endl;
	return -1;
   }
   TID_T taxid,depth;
   while(ifs1>>taxid>>depth) {
	sopt._imap[taxid] = depth;
   }
   

   StopWatch clock;
   clock.start();
   vector<map<TID_T,float> > track_tscoreall(n_threads); 
   vector<map<TID_T,int> > track_matchall(n_threads); 
   vector<map<nomatch_t,int> > track_nomatchall(n_threads); 
   size_t read_count = 0; 

#pragma omp parallel shared(arr, k_size, query_fn, ofbase, taxtable, tax_tree, sopt, prn_read,track_matchall,track_nomatchall,track_tscoreall,min_score,min_kmer)  private(ifs, finished, pos, ofs, ofname, line, read_count)
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
     outs << omp_get_thread_num();
     ofname += outs.str();
     ofname += ".out" ;

     ofs.open(ofname.c_str());

     ifs.seekg(arr[omp_get_thread_num()]);
     string read_buff, hdr_buff, save_hdr;
     while (!finished)   {
       getline(ifs, line);

       pos = ifs.tellg();

       if ((pos >= arr[1+omp_get_thread_num()]) || ((signed)pos == -1)) {
	      finished = true;
       } 
       if (line[0] == '>' || (fastq && line[0] == '@') ) {
         save_hdr=hdr_buff;
         // skip the ">"
         hdr_buff=line.substr(1,line.length()-1);
         //if(fastq) readOne=true;
       }
       if( finished ) {
         save_hdr=hdr_buff;
       }
       if (line[0] != '>' && line.length() > 1 && !fastq) {
          read_buff += line;
       } 
       if( fastq && line[0] != '@' && line[0] != '+' && line[0] != '-' ) {
          read_buff += line;
       }
       if( ((line[0] == '>' || finished) || (fastq && (line[0] == '+' || line[0] == '-'))) && read_buff.length() > 0 ) {
          if(save_hdr[0] == '\0') {
            ostringstream ostrm;
            ostrm<<"unknown_hdr:"<<read_count; 
            save_hdr=ostrm.str();
          }
          ofs<<save_hdr<<"\t";
          if( prn_read ) {
            ofs<<read_buff<<"\t";
          } else {
            ofs<<"X"<<"\t";
          } 
          int thread = omp_get_thread_num();
          map<TID_T,int>& track_match = track_matchall[thread]; 
          map<nomatch_t,int>& track_nomatch = track_nomatchall[thread]; 
          map<TID_T,float>& track_tscore = track_tscoreall[thread]; 
	       proc_line(tax_tree, read_buff.length(), read_buff, k_size, taxtable, ofs, threshold,sopt, max_count, track_match, track_nomatch, track_tscore,min_score, min_kmer);
          read_buff="";
	       read_count ++;
          if(fastq) getline(ifs, line); // skip quality values for now
       }
     }
     ofs.close();
   }

   map<TID_T,int>  merge_count;
   map<TID_T,float>  merge_score;
   for(unsigned thread = 0; thread < track_tscoreall.size();  ++thread) {
      map<TID_T,float>& gt = track_tscoreall[thread];
      map<TID_T,float>::const_iterator it = gt.begin();
      const map<TID_T,float>::const_iterator is = gt.end();
      for(; it != is; ++it) {
         TID_T tid = (*it).first;
         float score = (*it).second;
         if( merge_score.find(tid) == merge_score.end() ) {
            merge_score.insert( make_pair(tid,score));
         } else {
            merge_score[tid] += score;
         }
      }
      map<TID_T,int>& gtcnt = track_matchall[thread];
      map<TID_T,int>::const_iterator it1 = gtcnt.begin();
      const map<TID_T,int>::const_iterator is1 = gtcnt.end();
      for(; it1 != is1; ++it1) {
         TID_T tid = (*it1).first;
         int cnt = (*it1).second;
         if( merge_count.find(tid) == merge_count.end() ) {
            merge_count.insert( make_pair(tid,cnt));
         } else {
            merge_count[tid] += cnt;
         }
      }
   }
   set<TID_T> cand_tid;
   vector< pair<TID_T,float> > sort_val(merge_score.size());
   map<TID_T,float>::const_iterator it = merge_score.begin();
   map<TID_T,float>::const_iterator is = merge_score.end();
   for(unsigned i = 0; it != is; ++it, ++i) {
      sort_val[i] = (*it);
      const TID_T tid = sort_val[i].first;
      if( cand_tid.find( tid ) ==  cand_tid.end()) {
         //cout<<"confirm: "<<tid<<endl;
         cand_tid.insert(tid);
      }
   }
   const unsigned buff_size = 200000;
   char buff[buff_size];
   map<TID_T,string> save_id;
   if(verbose) cout<<"rank read: "<<rank_ids<<endl;
   ifstream tax_strm(rank_ids.c_str());
   while(tax_strm.getline(buff,buff_size)) {
      string proc(buff);
      char* val = strtok(buff,"=,");
      //cout<<"major: "<<buff<<" "<<val<<endl;
      while( val != NULL ) {
         if( strcmp(val,"taxid")==0 ) {
            val = strtok(NULL,"=,");
            istringstream istrm(val);
            TID_T cid;
            istrm>>cid;
            //cout<<"found taxid: "<<val<<endl;
            if( cand_tid.find(cid) != cand_tid.end() ) {
               size_t pos = proc.rfind('\t');
               string id=proc.substr(pos+1,proc.length()-pos);
               //cout<<"rankcheck: "<<val<<" "<<cid<<" "<<pos<<" ["<<id<<"] ["<<proc<<"]"<<endl;
               save_id.insert( make_pair(cid,id) );
            }
            break;
         }
         val = strtok(NULL,"=,");
      }
   }
   ostringstream ostrm;
   ostrm<<ofbase<<"."<<min_score<<"."<<min_kmer<<".fastsummary";
   ofstream sum_ofs(ostrm.str().c_str());
   if( !sum_ofs ) {
      cout<<"Could not open for writing "<<ostrm.str()<<endl;
      return -1;
   }
   sort(sort_val.begin(),sort_val.end(),SimpleCmp());
   for(unsigned i= 0; i < sort_val.size(); ++i) {
      const TID_T tid = sort_val[i].first;
      assert( merge_count.find( tid ) != merge_count.end());
      const int cnt = (*merge_count.find(tid)).second;
      const float wght_cnt = sort_val[i].second;
      string str_id = save_id[tid];
      sum_ofs<<wght_cnt<<"\t"<<cnt<<"\t"<<tid<<"\t"<<str_id<<endl;
   }

   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
