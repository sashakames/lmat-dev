/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory. 
 * Written by LMAT development team. LLNL-CODE-605872 All rights reserved.
 * This file is part of LMAT Please read COPYRIGHT file in root directory.
 * 
 * */

#include <unistd.h>
#include <fstream>
#include <math.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <omp.h>
#include "kencode.hpp"
#include "all_headers.hpp"
#include "TaxNodeStat.hpp"

#define MMAP_SIZE 0

//#define TID_T uint16_t
#define TID_T uint32_t

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
//typedef pair<uint16_t,tax_data_t> label_info_t;
typedef pair<int16_t,tax_data_t> label_info_t;
typedef map<TID_T,TID_T> hmap_t;
typedef map<TID_T,float> ufmap_t;
typedef map<uint16_t, ufmap_t> u_ufmap_t;

vector <int> read_len_vec;
vector <int> read_len_avgs;


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
struct ClosestCmp {
   bool operator()(const int & x, const int & y) { return x > y; }
};
*/

// vec is assumed to be sorted
int closest(int value)
{
//cout << "read_len_avgs.size: " << read_len_avgs.size() << " read_len_vec: " << read_len_vec.size() << endl;
  size_t i;
  for (i =0; i<read_len_avgs.size(); i++) {
    if (value <= read_len_avgs[i]) 
      return read_len_vec[i];
  }

  if (read_len_vec.size() > 0)
    return (read_len_vec[i]);
  else
    return 0;
  
  return -1;
}



static int getReadLen(int rl)
{

  int len = closest(rl);

  if (len > 0)
    return len;

  return 80;

}


typedef std::pair<TaxNode<TID_T>*,uint16_t> tncpair_t;

static bool isAncestor(const TaxTree<TID_T> & tax_tree, TID_T prev_taxid /*ancestor */, TID_T curr_taxid /* descendant */ ) {
   bool isOk=false;
   vector<TID_T> path;
   TaxTree<TID_T> & tax_tree_tmp = const_cast<TaxTree<TID_T> &>(tax_tree);
   tax_tree_tmp.getPathToRoot(curr_taxid,path);
   for(unsigned p = 0; p < path.size(); ++p) {
      if(path[p] == prev_taxid) {
         isOk=true;
         break;
      }
   }
   return isOk;
}


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
   //defualt:
      //cerr<<"Unexpected error: "<<match<<endl;
      //break;
   }
   return str;
}

static bool addToCandLineage(const ufpair_t cand, list<ufpair_t>& lineage, const hmap_t& dmap, const TaxTree<TID_T> & tax_tree) {
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

static bool cmpCompLineage(ufpair_t cand, const vector<ufpair_t>& lineage, set<TID_T>& no_good, float diff_thresh, const TaxTree<TID_T> & tax_tree) {
   const float undef = -10000;
   bool keep_going=true;
   for(unsigned i = 0; i < lineage.size(); ++i) {
      if( isAncestor(tax_tree, lineage[i].first, cand.first )) {
         break;
      } 
      // half to swap == test, because stdev can be 0 so first test must be > not >=
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

static pair<ufpair_t,match_t> findReadLabelVer2(const vector<ufpair_t>& rank_label, float diff_thresh, const TaxTree<TID_T> & tax_tree, const hmap_t& taxid2idx, list<ufpair_t>& cand_lin, const hmap_t& dmap, const ufmap_t& all_cand_set) {
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
      TaxTree<TID_T> & tax_tree_tmp = const_cast<TaxTree<TID_T> &>(tax_tree);
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
      // fill in lineage up to root here
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
   return make_pair(taxid_call,match);
}

void
fill_in_labels(const TaxTree<TID_T> & taxtree, vector<TID_T>& row, const tax_data_t& taxids, const hmap_t& idx2taxid, const hmap_t& taxid2idx) {
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const TID_T col_tax_val = (*mtch).second;
      const TaxTree<TID_T>::const_iterator col_tax_val_it = taxtree.find(col_tax_val); 
      if(col_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
         continue;
      }
      TaxNode<TID_T> * col_tax_val_node = (*col_tax_val_it).second;
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
            TaxNode<TID_T> * tax_val_node = (*tax_val_it).second;
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
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(true), _strict_kmer_match(false), _prn_all(true), _imap(imap), _diff_thresh(1.0) {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh;
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
      cout<<"load: "<<read_len<<" "<<file<<endl;
      read_len_vec.push_back(read_len);


      ifstream ifs(file.c_str());
      if( !ifs) {
         cerr<<"Unexpected reading error: "<<file<<endl;
         continue;
      }
      const int buff_size=20004;
      char buff[buff_size];
      while(ifs.getline(buff,buff_size)) {
         istringstream istrm(buff);
         int32_t taxid;
         float min_val=0,max_val=0,mean_val=0,stdev=0;
         istrm>>taxid>>min_val>>max_val>>mean_val>>stdev;
         //const float cutoff = mean_val+(3.0*stdev); // can play around with how cutoff is set
         const float cutoff = max_val;
         rand_hits_all[read_len][taxid]=cutoff;
      }
   } 
   
   for (size_t i = 1; i < read_len_vec.size(); i++) {
     read_len_avgs.push_back((read_len_vec[i-1] + read_len_vec[i]) / 2);
   }
     

}

static float log_odds_score(float label_prob, float random_prob) {
      if( label_prob >= 1 ) label_prob = 0.9999;
      if( random_prob >= 1 ) random_prob = 0.9999;
      const float denom = (1-label_prob)*random_prob; // must be > zero
      const float numer =  label_prob*(1-random_prob); // must be > zero
      if(verbose)  cout<<"log_odds_calc: "<<numer<<" "<<denom<<endl;
      const float log_odds = log( numer / denom );
      return log_odds;
}

void 
construct_labels(const TaxTree<TID_T> & tax_tree, const vector<label_info_t>& label_vec, const list<TID_T>& taxid_lst, const hmap_t& tax2idx, const hmap_t& idx2taxid, ofstream& ofs, size_t mer_len, const ScoreOptions& sopt, const vector<bool> & is_forward) {
   const unsigned num_tax_ids = taxid_lst.size();
   vector<bool> any_kmer_match(label_vec.size(),false) ;
   unsigned cnt_fnd_kmers=0, coverage = 0;
   int last_hit_pos = -1;
   vector<vector<TID_T> > label_matrix(label_vec.size());
   uint16_t cand_kmer_cnt = 0;
   for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
      if(label_vec[pos].first >= 0 ) ++cand_kmer_cnt;
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
            //fill_in_labels(tax_tree,label_matrix[pos],label_vec[pos].second,idx2taxid,tax2idx);
            ++cnt_fnd_kmers;
            if(last_hit_pos == -1 ) {
               coverage += mer_len;
            } else if( pos < (last_hit_pos + mer_len) ) {
               const int overlap = (last_hit_pos + mer_len) - pos;
               const int add_coverage = mer_len - overlap;
               coverage += add_coverage;
            } else { // pos >= (last_hit_pos+mer_len 
               coverage += mer_len;
            }
            last_hit_pos=pos;
         }
   }
      const unsigned read_len = label_matrix.size()+ mer_len - 1;
      const float pcnt_cov = (float)coverage /(float)(read_len);
      vector<ufpair_t> rank_label(num_tax_ids,make_pair(0,0)); 
      unsigned sig_hits = 0;
      float log_sum=0.0;
      ufmap_t all_cand_set;
      for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         float tot_genome_cnt = 0, found_genome_cnt = 0;
         const TID_T taxid = (*(idx2taxid.find(tax_idx))).second;
         if(verbose) cout<<"col: "<<tax_idx<<" "<<taxid;
         bool noMatch=false;
         unsigned add_coverage = 0;
         signed tax_last_hit_pos = -1;
         for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
            if( any_kmer_match[pos] ) {
               if( sopt._strict_kmer_match && label_matrix[pos][tax_idx] == 0 ) {
                  noMatch=true;
                  break;
               }
               tot_genome_cnt += 1; // equal weight;
               if( label_matrix[pos][tax_idx] > 0 ) {
                     found_genome_cnt += 1;
               }
               int local_add=0;
               if( label_matrix[pos][tax_idx] > 0 ) {
                  int overlap = -1;
                  if(tax_last_hit_pos == -1 ) {
                     local_add = mer_len;
                  } else if( pos < (tax_last_hit_pos + mer_len) ) {
                     overlap = (tax_last_hit_pos + mer_len) - pos;
                     if((int)mer_len >= overlap) {
                        local_add = (mer_len - overlap);
                     }
                  } else { // pos >= (last_hit_pos+mer_len 
                     local_add = mer_len;
                  }
                  add_coverage += local_add;
                  //cout<<"ddebug: "<<local_add<<" "<<add_coverage<<" "<<overlap<<" "<<mer_len<<" "<<tax_last_hit_pos<<endl;
                  tax_last_hit_pos=pos;
               }
               if(verbose) cout<<" pos="<<pos<<" "<<label_vec[pos].first<<" "<<label_matrix[pos][tax_idx]<<" taxid="<<taxid<<" "<<found_genome_cnt<<" "<<tot_genome_cnt<<" "<<local_add<<" "<<add_coverage<<endl;
            }
         }
         float label_prob = 0.0;
         if( !noMatch ) {
            if( sopt._equal_kmer_vote) {
               label_prob = (float)found_genome_cnt / (float) label_matrix.size();
            } else {
               label_prob = (float)add_coverage / (float)read_len;
            }
         }
         if(verbose) cout<<" "<<label_prob<<endl;
         float log_odds=label_prob;
         float random_prob = 0.0001;
         const  int cand_kmer_cnt_match = getReadLen(cand_kmer_cnt);
         if( !sopt._comp_rand_hits ) {
            u_ufmap_t::const_iterator mtch = sopt._rand_hits.find(cand_kmer_cnt_match);
            if( mtch != sopt._rand_hits.end() ) {
               const ufmap_t& rand_hits = (*mtch).second;
               if( rand_hits.find(taxid) != rand_hits.end()) {
                  random_prob = (*(rand_hits.find(taxid))).second+0.0001;
               }
            }
            //cout<<"Input: to los: "<<label_prob<<" "<<random_prob<<endl;
            log_odds = log_odds_score(label_prob,random_prob);
         }           
         rank_label[tax_idx] = make_pair(taxid,log_odds);
         all_cand_set.insert( rank_label[tax_idx] );
         log_sum += log_odds;
         sig_hits++;
         if(verbose) cout<<"Sig match chec: "<<taxid<<" "<<rank_label[tax_idx].second<<" "<<log_odds<<" "<<label_prob<<" "<<random_prob<<endl;
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
            const float val = log_avg - rank_label[tax_idx].second;
            log_std += (val*val);
         }
      }
      float stdev = sig_hits > 1 ? sqrt(log_std/(sig_hits-1)) : 0;
      if( sig_hits > 0 ) {
         TCmp tcmp(sopt._imap);
         sort(rank_label.begin(),rank_label.end(),tcmp);
         ofs<<"read_label "<< cnt_fnd_kmers<<" "<<pcnt_cov<<" "<<max_score<<" "<<log_avg<<" "<<stdev<<" ";
         stdev *= sopt._diff_thresh;
         res = findReadLabelVer2(rank_label,stdev,tax_tree,tax2idx,valid_cand,sopt._imap,all_cand_set);    
         if(sopt._prn_all) {
            for(signed i = rank_label.size()-1; i >= 0; --i) {
               if( rank_label[i].second >= 0  || verbose ) {
                  ofs<<" "<<rank_label[i].first<<" "<<rank_label[i].second;
               }
            }
         }
         matchType=match2str(res.second);
      }
      if( res.second == eDirectMatch) {
         const ufpair_t best_guess = res.first;
         ofs<<" "<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }  else if( res.second == eMultiMatch || res.second == ePartialMultiMatch ) {
         if( !sopt._prn_all) {
            list<ufpair_t>::const_iterator it = valid_cand.begin();
            const list<ufpair_t>::const_iterator is = valid_cand.end();
            for( ; it != is; ++it) {
               ofs<<" "<<(*it).first<<" "<<(*it).second;
            }
         }
         const ufpair_t best_guess = res.first;
         ofs<<" "<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }
      else if( eNoMatch ) {
         ofs<<" "<<-1<<" "<<-1<<" "<<matchType;
      } else {
         cerr<<"Unexpected match type"<<endl;
      }

      ofs<<endl;
      
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

   void retrieve_kmer_labels(INDEXDB<TID_T>* table, const string& read, const size_t mer_len, vector<label_info_t>& label_vec, list<TID_T>& taxid_lst, hmap_t& tax2idx, hmap_t& idx2tax, bool revStrnd, vector<bool>& is_forward, const hmap_t& dmap, const TaxTree<TID_T> & tax_tree, int16_t max_count) 
   {

     //  bool verbose=true;
     kencode_c ken(mer_len);
     size_t read_len = read.length();
     uint64_t kmer_id = 0x0;
     bool bad_last_kmer=true;
     unsigned next_valid_pos = 0;
     
      
     for(size_t j = 0; j < read.size(); ++j) {
       if(verbose) cout<<"retrieve kmer at posit:"<<j<<endl;
       if(!isBase(read[j])) {
           bad_last_kmer=true;
           next_valid_pos = (j + mer_len);
           continue;
       }
       if( bad_last_kmer && (j >= mer_len - 1) && j >= next_valid_pos ) {
           const string nkmer=read.substr((j-mer_len)+1,mer_len);
           kmer_id = ken.kencode(nkmer);
           if(verbose) cout<<"create new kmer at posit:"<<j<<endl;
           bad_last_kmer=false;
       } else if( !bad_last_kmer) {
           kmer_id = ken.kencode(read[j]);
           bad_last_kmer=false;
           if(verbose) cout<<"update kmer at posit:"<<j<<endl;
       }
       size_t pos = j - mer_len + 1;
       if (revStrnd) {
          pos = read_len - j - 1;
       }
       if(!bad_last_kmer /*&& label_vec[pos].first == 0*/) {
          label_vec[pos].first = 0; // marks the position as having a valid k-mer (for case where n-masked reads need to be identified)
          TaxNodeStat<TID_T> *h = new TaxNodeStat<TID_T>(*table);
          if(verbose) cout<<"lookup kmer at posit:"<<j<<endl;
          h->begin(kmer_id);
          unsigned dcnt = 0, mtch = 0;
          list<TID_T> obs_tids;
          while( h->next() ) {
            TID_T tid = h->taxid();
            if( tid == 20999999) continue;
            uint16_t ng = h->taxidCount();
            //start: hack to drop kmers with too many tids
            if (ng > max_count && max_count != 0) {
              if (add_root_on_kmer_drop ) { //start: assign root
                tid = 1;
                label_vec[pos].first = 1; 
                const uint16_t pr_cnt = 1; 
                obs_tids.push_back(tid);
                label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
                is_forward[pos] = !revStrnd;
                if( tax2idx.find(tid) == tax2idx.end() ) {
                  const unsigned idx = taxid_lst.size();
                  //cout<<"debug0 tax2idx: "<<tid<<" "<<idx<<endl;
                  tax2idx[tid] = idx;
                  idx2tax[idx] = tid;
                  taxid_lst.push_back(tid);
                }
              } //end: assign root
              break;
            }  //end: hack to drop kmers with too many tids
            if(dcnt==0) {
               assert( ng > 0 );
               if( ng <= 0 ) {
                  cout<<"Warning unexpected value for ng: "<<ng<<endl;
                  ng=1;
               }
               label_vec[pos].first = ng;
            }
            //collect the unique set of tax ids
            const uint16_t pr_cnt = 1;
            //const uint16_t pr_cnt = h->present();
               obs_tids.push_back(tid);
               label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
               is_forward[pos] = !revStrnd;
               if( tax2idx.find(tid) == tax2idx.end() ) {
                  const unsigned idx = taxid_lst.size();
                  //cout<<"debug1 tax2idx: "<<tid<<" "<<idx<<endl;
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
                  TaxTree<TID_T> & tax_tree_tmp = const_cast<TaxTree<TID_T> &>(tax_tree);
                  tax_tree_tmp.getPathToRoot(tid,path);
                  for(unsigned p = 0; p < path.size(); ++p) {
                     const TID_T ptid = path[p];
                     label_vec[pos].second.insert( make_pair(ptid,1) );
                     is_forward[pos] = !revStrnd;
                     if( tax2idx.find(ptid) == tax2idx.end() ) {
                        const unsigned idx = taxid_lst.size();
                        tax2idx[ptid] = idx;
                        idx2tax[idx] = ptid;
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
   }

   void proc_line(const TaxTree<TID_T> & tax_tree, int ri_len, string &line, int k_size, INDEXDB<TID_T> *table, ofstream &ofs, const ScoreOptions& sopt, int16_t max_count) {

     //int thread_num = omp_get_thread_num();
     if(ri_len < 0 || ri_len > (signed)line.length()) {
         cout<<"unexpected ri_len value: "<<ri_len<<endl;
         return;
      }

       for(unsigned chk = 0; (chk+ri_len) <= line.length(); chk += (ri_len+1)) {
          string rval;
          if( ri_len < (signed)line.length() ) {
            rval = line.substr(chk,ri_len);
          } else {
            rval = line;
          }
          const unsigned read_len = rval.length();
          if( static_cast<signed>(read_len) < ri_len ) {
            break;
          }  
          if( (signed)read_len < k_size ) {
            if( chk == 0 ) /// read too short
               ofs<<"read_label 0 -1 ReadTooShort "<<read_len<<" "<<k_size<<endl; 
            break;
          } 
          vector<label_info_t> label_vec(read_len-k_size+1,make_pair(-1,tax_data_t()));
          vector<bool> is_forward(read_len-k_size+1,false);
          list<TID_T> taxid_lst; 
          hmap_t tax2idx, idx2tax;

          retrieve_kmer_labels(table, rval, k_size, label_vec, taxid_lst, tax2idx, idx2tax, false, is_forward, sopt._imap, tax_tree, max_count);

          string rev_cmp(rval.size(),'\0');
          assert(rval.size()>0);
          for (int j=(rval.size()-1); j >= 0; --j) {
            const int bidx = (rval.size()- 1) - j;
            rev_cmp[j] = getComp(rval[bidx]);
          }
          retrieve_kmer_labels(table, rev_cmp, k_size, label_vec, taxid_lst, tax2idx, idx2tax, true, is_forward, sopt._imap, tax_tree, max_count);

          if( !taxid_lst.empty() ) {  
            construct_labels(tax_tree,label_vec,taxid_lst,tax2idx,idx2tax,ofs,k_size,sopt,is_forward); 
         } else {
            ofs<<"read_label 0 -1 NoDbHits"<<endl; 
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
  cout << "Usage:\n"
  << "  -d <string>  - database filename       [required]\n"
  << "  -o <string>  - output base filename    [required]\n"
  << "  -i <string>  - query filename          [required]\n"
  << "  -c <string>  - taxonomy filename       [required]\n"
  << "  -e <string>  - taxonomy depth filename [required]\n"
  << "  -k           - set k-mer length        [required]\n"
  << "  -m <string>  - mapping for 16 - 32 bit tids[optional]\n"
  << "  -n <string>  - null model filelist     [optional]\n"
  << "  -r <string>  - if given, <string> is a file containing a list of\n"
  << "                 tax_histo files to load; if absent, DB will be restored\n"
  << "                 from memory mapped file; if given, do not also give '-d'\n"
  << "                 [optional; default: restore mmapped file]\n"
  << "  -t <int>     - thread count            [optional; default: 1]\n"
  << "  -b <float>   - scoring differential between selected tax id and competing alternatives\n"
  << "                 (number of standard deviations difference) [optional; default: 1]\n"
  << "  -h <int>     - ignore k-mers if the map to more than this number of taxonomy nodes;\n"
  << "                 [optional; default: 50]\n"
  << "  -p           - print out all competing labels with score 0 or greater [optional]\n"
  << "\n"
  << "example invocation:\n"
  << execname << " -d marker_lib.18 -o output.marker.k=18 -i my_queries.fa -c ../runtime_inputs/kpath_taxonomy.dat -e ../runtime_inputs/depth_for_kpath_taxonomy.dat  -t 32\n"

  << "\nnotes: -n option is not necessary when using our supplied reference DB (available\n"
  <<   "       for download from gdo-bioinformatics.ucllnl.org/pub/lmat)\n";
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1, ri_len = -1;
   int n_threads = 1;
   bool ascii = false;
   omp_set_num_threads(1);


   string kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, depth_file, rand_hits_file;
   string map_file_32_16;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   size_t mmap_size = 0;
   
   int max_reads =0;
   int16_t max_count = 50;
   int count = 0;
   string tax_histo_list;

   while ((c = getopt(argc, argv, "r:h:n:jb:ye:m:pk:c:v:k:i:d:l:t:s:r o:x:f:g:z:q:a ")) != -1) {
      switch(c) {
      case 'm':
         map_file_32_16 =optarg;
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
      case 'b':
         sopt._diff_thresh = atof(optarg);
         break;
      case 'e':
         ++count;
         depth_file = optarg;
         break;
      case 'q':
	max_reads = atoi(optarg);
	break;
      case 'p':
         sopt._prn_all = true;
         break;   
      case 'z':
	verbose = true;
	break;
      case 'a':
	ascii = true;
	break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'l':
         ri_len = atoi(optarg);
         break;
      case 'c':
         ++count;
         tax_tree_fn = optarg;
	 break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'i':
         ++count;
         query_fn = optarg;
         break;
      case 'd':
         ++count;
         kmer_db_fn = optarg;
         break;
      case 'o':
         ++count;
         ofbase = optarg;
         break;
      case 'r':
         ++count;
         tax_histo_list = optarg;
         break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }


   if (count != 5) {
     cout << " count: " << count << endl;
     usage(argv[0]);
     return -1;
   }

   INDEXDB<TID_T> *taxtable;
   if (!tax_histo_list.size()) {

     cout << "starting restore for " << kmer_db_fn << endl;

#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<INDEXDB<TID_T>*, std::size_t> ret = mfile.find<INDEXDB<TID_T> >("KmerDB");

     taxtable = ret.first;
     taxtable->conv_ptrs();

#else



     perm(&taxtable, sizeof(taxtable));
     mopen(kmer_db_fn.c_str(), "r", mmap_size);
     //     mopen(kmer_db_fn.c_str(), "r+", mmap_size);
#endif

     if (k_size < 1) {
       k_size = taxtable->get_kmer_length();
     }  

     cout << "num kmers: " << taxtable->size() << " - " << k_size  <<   endl ;

     } else {

       cout << "loading tax histo files\n";

       string line;
       ifstream in(tax_histo_list.c_str());
       assert(in);

#if USE_SORTED_DB == 1

#if USE_BOOST == 1
       cout << "Recompile with BOOST=0!" << endl;
       exit(0);
#else
       taxtable = new INDEXDB<TID_T>(5000000000, 256000000000);
#endif
       while (true) {

         line = "";
         getline(in, line);
         if (! line.size()) break;
         cout << ">>> "<< line <<endl;
	 my_map mm;
         taxtable->add_data(line.c_str(), 0, true, mm, 0 );
       }


#else
       taxtable = new INDEXDB<TID_T>;

       while (true) {
         line = "";
         getline(in, line);
         if (! line.size()) break;
         cout << ">>> "<< line <<endl;
         taxtable->registerFile(line.c_str());
       }
       taxtable->ingest();
#endif

     }


   sopt._comp_rand_hits = (rand_hits_file.length() == 0);
   if( !sopt._comp_rand_hits ) {
      loadRandHits(rand_hits_file,sopt._rand_hits);
   }

   assert(k_size > 0 );
   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   size_t *arr = split_file(n_threads, ifs);
   ifs.close();

   bool finished;

   size_t pos = 0 ;

   ofstream ofs;

   cout<<"Read taxonomy tree: "<<tax_tree_fn<<endl;
   TaxTree<TID_T>  tax_tree(tax_tree_fn.c_str());

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
   

  
   size_t read_count = 0; 

#pragma omp parallel shared(arr, k_size, query_fn, ofbase, taxtable, tax_tree, max_reads, sopt, ri_len)  private(ifs, finished, pos, ofs, ofname, line, read_count)

   {


     StopWatch clock;
     clock.start();


     read_count = 0;

     

     finished = false;
     int use_len = -1; 

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

     string read_buff;
int ct = 0;
     while (!finished)   {
++ct;
//if (ct == 100) break;


       getline(ifs, line);

       pos = ifs.tellg();

       if ((pos >= arr[1+omp_get_thread_num()]) || ((signed)pos == -1)) {
	      finished = true;
       } 

       //#pragma omp critical    
       //cout << omp_get_thread_num()  << " - pos:" << ifs.tellg() << "\n"; 

       if (line[0] != '>' && line.length() > 1) {
          read_buff += line;
       } 
       if( (line[0] == '>' || finished) && read_buff.length() > 0 ) {
          if( ri_len == -1 ) {
            // if not specified use first read as the fixed length
            use_len = read_buff.length();
          } else {
            use_len = ri_len;
          } 
	       //#pragma omp critical
          ofs<<"read: "<<read_buff<< " pos: " << ifs.tellg() <<  endl;
       //split read into multiples of "ri_len" (read interval length)
	  proc_line(tax_tree, use_len, read_buff, k_size, taxtable, ofs, sopt, max_count);
          read_buff="";
	  read_count ++;
	  if ((signed)read_count == max_reads)
	   finished = true;

       }
     }
     ofs.close();

     cout << "threadtime: " << clock.stop() << endl;
   }



   return 0; 
}
